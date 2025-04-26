import pandas as pd
import ray
import duckdb
import sqlite3
import pyranges.methods.merge, pyranges.methods.intersection, pyranges.methods.coverage, pyranges.methods.subtraction

def get_connection(sql_db_name, backend="duckdb"):
    if backend == "duckdb":
        return duckdb.connect(sql_db_name, read_only=True)
    else:
        return sqlite3.connect(sql_db_name)

def query_db(query, conn, backend="duckdb"):
    if backend == "duckdb":
        return conn.execute(query).fetchdf()
    else:
        return pd.read_sql_query(query, conn)
    

def count_exons(sql_table_name, conn, backend="duckdb"):
    return query_db(f"SELECT gene_id, COUNT(*) as exon_count FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id", conn, backend)

def exon_length(sql_table_name, conn, backend="duckdb"):
    if backend == "duckdb":
        query = f"SELECT gene_id, SUM(\"End\" - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id"
    else:
        query = f"SELECT gene_id, SUM(End - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id"
    return query_db(query, conn, backend)

def highest_transcripts(sql_table_name, conn, backend="duckdb"):
    if sql_table_name == "arabidopsis":
        return query_db(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature IN (\'mRNA\', \'snoRNA\', \'lnc_RNA\', \'pseudogenic_transcript\') GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn, backend)
    else:
        return query_db(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn, backend)

@ray.remote
def  merge_exon_intervals_single(sql_table_name, sql_db_name, chrom_strand, backend="duckdb"):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name, backend)
    if backend == "duckdb":
        query = f"SELECT \"Start\", \"End\" FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'"
    else:
        query = f"SELECT Start, End FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'"
    exon_intervals = query_db(query, conn, backend)
    conn.close()
    exon_intervals = pyranges.methods.merge._merge(exon_intervals, chromosome=chrom, count=None, strand=strand)
    return exon_intervals

def merge_exon_intervals(sql_table_name, sql_db_name, chrom_strand_tup, backend="duckdb"):
    futures = [merge_exon_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand, backend=backend) for chrom_strand in chrom_strand_tup]
    exon_intervals = ray.get(futures)
    exon_intervals = pd.concat(exon_intervals).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals

@ray.remote
def get_overlapping_genes_single(sql_table_name, sql_db_name, chrom_strand, other_genes, backend="duckdb"):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name, backend)
    self_genes = query_db(f"SELECT * FROM {sql_table_name} WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn, backend)
    conn.close()

    other_genes_chrom_strand = other_genes[
        (other_genes["Chromosome"] == chrom) & (other_genes["Strand"] == strand)
    ]
    overlapping_genes = pyranges.methods.intersection._overlap(self_genes, other_genes_chrom_strand, how="first")
    return overlapping_genes

def get_overlapping_genes(sql_table_name, sql_db_name, chrom_strand_tup, other_genes, backend="duckdb"):
    futures = [get_overlapping_genes_single.remote(sql_table_name, sql_db_name, chrom_strand, other_genes, backend=backend) for chrom_strand in chrom_strand_tup]
    overlapping_genes_list = ray.get(futures)
    return pd.concat(overlapping_genes_list)

@ray.remote
def get_subtracted_exons_single(sql_table_name, sql_db_name, chrom_strand, other_cdf, backend="duckdb"):
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name, backend)
    self_genes = query_db(f"SELECT * FROM {sql_table_name} WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn, backend)
    conn.close()
    other_cdf_clusters = other_cdf.merge(strand=True)
    other_chrom_clusters = other_cdf_clusters[(other_cdf_clusters.Chromosome == chrom) & (other_cdf_clusters.Strand == strand)].df
    # add __num__ column to self_genes_sql which counts how many intervals in self_genes_sql overlap with other_genes
    self_genes = pyranges.methods.coverage._number_overlapping(self_genes, other_chrom_clusters, strandedness="same", keep_nonoverlapping=True, overlap_col="__num__")
    subtracted_intervals = pyranges.methods.subtraction._subtraction(self_genes, other_chrom_clusters, strandedness="same")
    if subtracted_intervals is not None:
        return subtracted_intervals.drop(columns=["__num__"])
    else:
        return None

def get_subtracted_exons(sql_table_name, sql_db_name, chrom_strand_tup, other_cdf, backend="duckdb"):
    futures = [get_subtracted_exons_single.remote(sql_table_name, sql_db_name, chrom_strand, other_cdf, backend=backend) for chrom_strand in chrom_strand_tup]
    subtracted_intervals = ray.get(futures)
    subtracted_intervals = pd.concat(subtracted_intervals)
    return subtracted_intervals