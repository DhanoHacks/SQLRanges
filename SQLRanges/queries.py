import pandas as pd
import ray
import duckdb
import sqlite3
import pyranges.methods.merge, pyranges.methods.intersection, pyranges.methods.coverage, pyranges.methods.subtraction

def get_connection(sql_db_name: str, backend: str = "duckdb") -> sqlite3.Connection | duckdb.DuckDBPyConnection:
    """Get a connection to the database.

    Args:
        sql_db_name (str): Name of the database file.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        sqlite3.Connection | duckdb.DuckDBPyConnection: Database connection object.
    """
    if backend == "duckdb":
        return duckdb.connect(sql_db_name, read_only=True)
    else:
        return sqlite3.connect(sql_db_name)

def query_db(query: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, backend: str = "duckdb") -> pd.DataFrame:
    """Execute a SQL query and return the result as a pandas DataFrame.

    Args:
        query (str): SQL query to execute.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: Result of the query as a pandas DataFrame.
    """
    if backend == "duckdb":
        return conn.execute(query).fetchdf()
    else:
        return pd.read_sql_query(query, conn)
    

def count_intervals(sql_table_name: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, group_by: str = "gene_id", feature_filter: None | str = None, return_col_name: str = "count", backend: str = "duckdb") -> pd.DataFrame:
    """Count the number of intervals in the database, grouped by a specified column.

    Args:
        sql_table_name (str): Name of the SQL table.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        group_by (str, optional): Column to group by. Defaults to "gene_id".
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        return_col_name (str, optional): Column name for the count result. Defaults to "count".
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the grouped counts.
    """
    if feature_filter is None:
        return query_db(f"SELECT {group_by}, COUNT(*) as {return_col_name} FROM {sql_table_name} WHERE Feature = \'{feature_filter}\' GROUP BY {group_by}", conn, backend)
    else:
        return query_db(f"SELECT {group_by}, COUNT(*) as {return_col_name} FROM {sql_table_name} GROUP BY {group_by}", conn, backend)

def exon_length(sql_table_name: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, backend: str = "duckdb") -> pd.DataFrame:
    """Calculate the total length of exons for each gene in the database.

    Args:
        sql_table_name (str): Name of the SQL table.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the gene IDs and their corresponding total exon lengths.
    """
    if backend == "duckdb":
        query = f"SELECT gene_id, SUM(\"End\" - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id"
    else:
        query = f"SELECT gene_id, SUM(End - Start) as total_exon_length FROM {sql_table_name} WHERE Feature = \'exon\' GROUP BY gene_id"
    return query_db(query, conn, backend)

def highest_transcripts(sql_table_name: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, backend: str = "duckdb") -> pd.DataFrame:
    """Identify the chromosome with the highest number of transcripts in the database.

    Args:
        sql_table_name (str): Name of the SQL table.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the chromosome name and its corresponding transcript count.
    """
    if sql_table_name == "arabidopsis":
        return query_db(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature IN (\'mRNA\', \'snoRNA\', \'lnc_RNA\', \'pseudogenic_transcript\') GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn, backend)
    else:
        return query_db(f"SELECT Chromosome, COUNT(*) as transcript_count FROM {sql_table_name} WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn, backend)

@ray.remote
def  merge_exon_intervals_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, backend="duckdb") -> pd.DataFrame:
    """Merge overlapping exon intervals for a specific chromosome and strand.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the merged exon intervals.
    """
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

def merge_exon_intervals(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, backend="duckdb") -> pd.DataFrame:
    """Merge overlapping exon intervals for multiple chromosomes and strands.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the merged exon intervals for all specified chromosomes and strands.
    """
    futures = [merge_exon_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand, backend=backend) for chrom_strand in chrom_strand_tup]
    exon_intervals = ray.get(futures)
    exon_intervals = pd.concat(exon_intervals).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals

@ray.remote
def get_overlapping_genes_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, other_genes: pd.DataFrame, backend="duckdb") -> pd.DataFrame:
    """Get overlapping genes for a specific chromosome and strand.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        other_genes (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the overlapping genes for the specified chromosome and strand.
    """
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name, backend)
    self_genes = query_db(f"SELECT * FROM {sql_table_name} WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn, backend)
    conn.close()

    other_genes_chrom_strand = other_genes[
        (other_genes["Chromosome"] == chrom) & (other_genes["Strand"] == strand)
    ]
    overlapping_genes = pyranges.methods.intersection._overlap(self_genes, other_genes_chrom_strand, how="first")
    return overlapping_genes

def get_overlapping_genes(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, other_genes: pd.DataFrame, backend="duckdb") -> pd.DataFrame:
    """Get overlapping genes for multiple chromosomes and strands.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        other_genes (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
        backend (str, optional): Database backend to use. Defaults to "duckdb".
        The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.

    Returns:
        pd.DataFrame: A DataFrame containing the overlapping genes for all specified chromosomes and strands.
    """
    futures = [get_overlapping_genes_single.remote(sql_table_name, sql_db_name, chrom_strand, other_genes, backend=backend) for chrom_strand in chrom_strand_tup]
    overlapping_genes_list = ray.get(futures)
    return pd.concat(overlapping_genes_list)

@ray.remote
def get_subtracted_exons_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, other_cdf: pd.DataFrame, backend="duckdb") -> pd.DataFrame:
    """Subtract overlapping exons for a specific chromosome and strand.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        other_cdf (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the subtracted exon intervals for the specified chromosome and strand.
    """
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

def get_subtracted_exons(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, other_cdf: pd.DataFrame, backend="duckdb") -> pd.DataFrame:
    """Subtract overlapping exons for multiple chromosomes and strands.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        other_cdf (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
        The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the subtracted exon intervals for all specified chromosomes and strands.
    """
    futures = [get_subtracted_exons_single.remote(sql_table_name, sql_db_name, chrom_strand, other_cdf, backend=backend) for chrom_strand in chrom_strand_tup]
    subtracted_intervals = ray.get(futures)
    subtracted_intervals = pd.concat(subtracted_intervals)
    return subtracted_intervals