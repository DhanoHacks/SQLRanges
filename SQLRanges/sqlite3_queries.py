import pandas as pd
import ray
import sqlite3
import pyranges.methods.merge, pyranges.methods.intersection, pyranges.methods.coverage, pyranges.methods.subtraction

def count_exons(sql_table_name, conn):
    return pd.read_sql_query(f"SELECT gene_id, COUNT(*) as exon_count FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id", conn)

def exon_length(sql_table_name, conn):
    return pd.read_sql_query(f"SELECT gene_id, SUM(End - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id", conn)

def highest_transcripts(sql_table_name, conn):
    if sql_table_name == "arabidopsis":
        return pd.read_sql_query(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature IN (\'mRNA\', \'snoRNA\', \'lnc_RNA\', \'pseudogenic_transcript\') GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)
    else:
        return pd.read_sql_query(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)
    
def get_connection(sql_db_name):
    return sqlite3.connect(sql_db_name)
    
@ray.remote
def merge_exon_intervals_single(sql_table_name, sql_db_name, chrom_strand):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name)
    exon_intervals_sql = pd.read_sql_query(f"SELECT Start, End FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    conn.close()
    exon_intervals_sql = pyranges.methods.merge._merge(exon_intervals_sql, chromosome=chrom, count=None, strand=strand)
    return exon_intervals_sql

def merge_exon_intervals(sql_table_name, sql_db_name, chrom_strand_tup):
    futures = [merge_exon_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand) for chrom_strand in chrom_strand_tup]
    exon_intervals_sql = ray.get(futures)
    exon_intervals_sql = pd.concat(exon_intervals_sql).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals_sql

@ray.remote
def get_overlapping_genes_single(sql_table_name, sql_db_name, chrom_strand, other_genes):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name)
    self_genes_sql = pd.read_sql_query(f"SELECT * FROM {sql_table_name} WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    conn.close()

    other_genes_chrom_strand = other_genes[
        (other_genes["Chromosome"] == chrom) & (other_genes["Strand"] == strand)
    ]
    overlapping_genes_sql = pyranges.methods.intersection._overlap(self_genes_sql, other_genes_chrom_strand, how="first")
    return overlapping_genes_sql

def get_overlapping_genes(sql_table_name, sql_db_name, chrom_strand_tup, other_genes):
    futures = [get_overlapping_genes_single.remote(sql_table_name, sql_db_name, chrom_strand, other_genes) for chrom_strand in chrom_strand_tup]
    overlapping_genes_sql_list = ray.get(futures)
    return pd.concat(overlapping_genes_sql_list)

@ray.remote
def get_subtracted_exons_single(sql_table_name, sql_db_name, chrom_strand, other_cdf):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name)
    self_genes_sql = pd.read_sql_query(f"SELECT * FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    conn.close()
    other_cdf_clusters = other_cdf.merge(strand=True)
    other_chrom_clusters = other_cdf_clusters[(other_cdf_clusters.Chromosome == chrom) & (other_cdf_clusters.Strand == strand)].df
    # add __num__ column to self_genes_sql which counts how many intervals in self_genes_sql overlap with other_genes
    self_genes_sql = pyranges.methods.coverage._number_overlapping(self_genes_sql, other_chrom_clusters, strandedness="same", keep_nonoverlapping=True, overlap_col="__num__")
    subtracted_intervals_sql = pyranges.methods.subtraction._subtraction(self_genes_sql, other_chrom_clusters, strandedness="same")
    if subtracted_intervals_sql is not None:
        return subtracted_intervals_sql.drop(columns=["__num__"])
    else:
        return None

def get_subtracted_exons(sql_table_name, sql_db_name, chrom_strand_tup, other_cdf):
    futures = [get_subtracted_exons_single.remote(sql_table_name, sql_db_name, chrom_strand, other_cdf) for chrom_strand in chrom_strand_tup]
    subtracted_intervals_sql = ray.get(futures)
    subtracted_intervals_sql = pd.concat(subtracted_intervals_sql) #.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return subtracted_intervals_sql