import pandas as pd
import ray
import duckdb
import pyranges.methods.merge, pyranges.methods.intersection, pyranges.methods.coverage, pyranges.methods.subtraction

def count_exons(sql_table_name, conn):
    return conn.execute(f"SELECT gene_id, COUNT(*) as exon_count FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id").fetchdf()

def exon_length(sql_table_name, conn):
    return conn.execute(f"SELECT gene_id, SUM(\"End\" - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id").fetchdf()

def highest_transcripts(sql_table_name, conn):
    if sql_table_name == "arabidopsis":
        return conn.execute(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature IN (\'mRNA\', \'snoRNA\', \'lnc_RNA\', \'pseudogenic_transcript\') GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1").fetchdf()
    else:
        return conn.execute(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1").fetchdf()

def get_connection(sql_db_name):
    return duckdb.connect(sql_db_name, read_only=True)

@ray.remote
def  merge_exon_intervals_single(sql_table_name, sql_db_name, chrom_strand):
    chrom, strand = chrom_strand
    conn = duckdb.connect(sql_db_name, read_only=True)
    exon_intervals_duckdb = conn.execute(f"SELECT \"Start\", \"End\" FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'").fetchdf()
    conn.close()
    exon_intervals_duckdb = pyranges.methods.merge._merge(exon_intervals_duckdb, chromosome=chrom, count=None, strand=strand)
    return exon_intervals_duckdb

def merge_exon_intervals(sql_table_name, sql_db_name, chrom_strand_tup):
    futures = [merge_exon_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand) for chrom_strand in chrom_strand_tup]
    exon_intervals_duckdb = ray.get(futures)
    exon_intervals_duckdb = pd.concat(exon_intervals_duckdb).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals_duckdb

@ray.remote
def get_overlapping_genes_single(sql_table_name, sql_db_name, chrom_strand, other_genes):
    chrom, strand = chrom_strand
    conn = duckdb.connect(sql_db_name, read_only=True)
    self_genes_duckdb = conn.execute(f"SELECT * FROM {sql_table_name} WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'").fetchdf()
    conn.close()

    other_genes_chrom_strand = other_genes[
        (other_genes["Chromosome"] == chrom) & (other_genes["Strand"] == strand)
    ]
    overlapping_genes_duckdb = pyranges.methods.intersection._overlap(self_genes_duckdb, other_genes_chrom_strand, how="first")
    return overlapping_genes_duckdb

def get_overlapping_genes(sql_table_name, sql_db_name, chrom_strand_tup, other_genes):
    futures = [get_overlapping_genes_single.remote(sql_table_name, sql_db_name, chrom_strand, other_genes) for chrom_strand in chrom_strand_tup]
    overlapping_genes_duckdb_list = ray.get(futures)
    return pd.concat(overlapping_genes_duckdb_list)

@ray.remote
def get_subtracted_exons_single(sql_table_name, sql_db_name, chrom_strand, other_cdf):
    chrom, strand = chrom_strand
    conn = duckdb.connect(sql_db_name, read_only=True)
    self_genes_duckdb = conn.execute(f"SELECT * FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'").fetchdf()
    conn.close()
    other_cdf_clusters = other_cdf.merge(strand=True)
    other_chrom_clusters = other_cdf_clusters[(other_cdf_clusters.Chromosome == chrom) & (other_cdf_clusters.Strand == strand)].df
    # add __num__ column to self_genes_sql which counts how many intervals in self_genes_sql overlap with other_genes
    self_genes_duckdb = pyranges.methods.coverage._number_overlapping(self_genes_duckdb, other_chrom_clusters, strandedness="same", keep_nonoverlapping=True, overlap_col="__num__")
    subtracted_intervals_duckdb = pyranges.methods.subtraction._subtraction(self_genes_duckdb, other_chrom_clusters, strandedness="same")
    if subtracted_intervals_duckdb is not None:
        return subtracted_intervals_duckdb.drop(columns=["__num__"])
    else:
        return None

def get_subtracted_exons(sql_table_name, sql_db_name, chrom_strand_tup, other_cdf):
    futures = [get_subtracted_exons_single.remote(sql_table_name, sql_db_name, chrom_strand, other_cdf) for chrom_strand in chrom_strand_tup]
    subtracted_intervals_duckdb = ray.get(futures)
    subtracted_intervals_duckdb = pd.concat(subtracted_intervals_duckdb) #.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return subtracted_intervals_duckdb