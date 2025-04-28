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
    """Count the number of intervals in the database, grouped by a specified column. The function can also optionaly filter the intervals based on a specific feature.

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
        feature_clause = ""
    else:
        feature_clause = f" WHERE Feature = '{feature_filter}'"
    query = f"SELECT \"{group_by}\", COUNT(*) as \"{return_col_name}\" FROM \"{sql_table_name}\"{feature_clause} GROUP BY \"{group_by}\""
    return query_db(query, conn, backend)

def total_length(sql_table_name: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, group_by: str = "gene_id", feature_filter: None | str = None, return_col_name: str = "total_length", backend: str = "duckdb") -> pd.DataFrame:
    """Calculate the total length of intervals in the database, grouped by a specified column. The function can also optionaly filter the intervals based on a specific feature.

    Args:
        sql_table_name (str): Name of the SQL table.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        group_by (str, optional): Column to group by. Defaults to "gene_id".
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        return_col_name (str, optional): Column name for the total length result. Defaults to "total_length".
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the total length of intervals grouped by the specified column.
    """
    if feature_filter is None:
        feature_clause = ""
    else:
        feature_clause = f" WHERE Feature = '{feature_filter}'"
    query = f"SELECT \"{group_by}\", SUM(\"End\" - \"Start\") as \"{return_col_name}\" FROM \"{sql_table_name}\"{feature_clause} GROUP BY \"{group_by}\""
    return query_db(query, conn, backend)

@ray.remote
def  merge_intervals_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Merge intervals for a specific chromosome and strand.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the merged intervals for the specified chromosome and strand.
    """
    chrom, strand = chrom_strand
    conn = get_connection(sql_db_name, backend)
    if feature_filter is None:
        feature_clause = ""
    else:
        feature_clause = f" Feature = '{feature_filter}' AND"
    query = f"SELECT \"Start\", \"End\" FROM \"{sql_table_name}\" WHERE{feature_clause} Chromosome = '{chrom}' AND Strand = '{strand}'"
    merged_intervals = query_db(query, conn, backend)
    conn.close()
    merged_intervals = pyranges.methods.merge._merge(merged_intervals, chromosome=chrom, count=None, strand=strand)
    return merged_intervals

def merge_intervals(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Merge intervals for multiple chromosomes and strands.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the merged intervals for all specified chromosomes and strands.
    """
    futures = [merge_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand, feature_filter=feature_filter, backend=backend) for chrom_strand in chrom_strand_tup]
    merged_intervals = ray.get(futures)
    merged_intervals = pd.concat(merged_intervals).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return merged_intervals

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