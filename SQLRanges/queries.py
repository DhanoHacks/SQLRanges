import pandas as pd
import ray
import duckdb
import sqlite3
import tempfile
import os
import pyranges.methods.merge, pyranges.methods.intersection, pyranges.methods.coverage, pyranges.methods.subtraction
from .utils import get_connection, query_db, get_intervals, to_db

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
    self_intervals = get_intervals(sql_table_name, sql_db_name, chrom, strand, feature_filter=feature_filter, backend=backend)
    
    merged_intervals = pyranges.methods.merge._merge(self_intervals, chromosome=chrom, count=None, strand=strand)
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
    merged_intervals = [df for df in merged_intervals if df is not None and not df.empty]
    if not merged_intervals:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand"])
    merged_intervals = pd.concat(merged_intervals)
    return merged_intervals

@ray.remote
def overlapping_intervals_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, other_intervals: pd.DataFrame, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Find overlapping intervals between the database and a set of other intervals for a specific chromosome and strand. The function can also optionaly filter the intervals based on a specific feature.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        other_intervals (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand' (and 'Feature' if feature_filter is set).
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the overlapping intervals for the specified chromosome and strand.
    """
    chrom, strand = chrom_strand
    self_intervals = get_intervals(sql_table_name, sql_db_name, chrom, strand, feature_filter=feature_filter, backend=backend)

    if feature_filter is None:
        other_intervals_chrom_strand = other_intervals[
            (other_intervals["Chromosome"] == chrom) & (other_intervals["Strand"] == strand)
        ]
    else:
        other_intervals_chrom_strand = other_intervals[
            (other_intervals["Chromosome"] == chrom) & (other_intervals["Strand"] == strand) & (other_intervals["Feature"] == feature_filter)
        ]
    overlapping_intervals = pyranges.methods.intersection._overlap(self_intervals, other_intervals_chrom_strand, how="first")
    return overlapping_intervals

def overlapping_intervals(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, other_intervals: pd.DataFrame, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Find overlapping intervals between the database and a set of other intervals for multiple chromosomes and strands. The function can also optionaly filter the intervals based on a specific feature.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        other_intervals (pd.DataFrame): DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand' (and 'Feature' if feature_filter is set).
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the overlapping intervals for all specified chromosomes and strands.
    """
    futures = [overlapping_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand, other_intervals, feature_filter=feature_filter, backend=backend) for chrom_strand in chrom_strand_tup]
    overlapping_intervals_list = ray.get(futures)
    return pd.concat(overlapping_intervals_list)

@ray.remote
def subtract_intervals_single(sql_table_name: str, sql_db_name: str, chrom_strand: tuple, other_intervals: pd.DataFrame, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Subtract a set of other intervals from the database intervals for a specific chromosome and strand. The function can also optionaly filter the intervals based on a specific feature.
    This function is designed to be run in parallel using Ray.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand (tuple): Tuple containing chromosome and strand information.
        other_intervals (pd.DataFrame): DataFrame containing the gene intervals to subtract from the database.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand' (and 'Feature' if feature_filter is set).
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the subtracted intervals for the specified chromosome and strand.
    """
    chrom, strand = chrom_strand
    self_intervals = get_intervals(sql_table_name, sql_db_name, chrom, strand, feature_filter=feature_filter, backend=backend)
    
    temp_path = tempfile.NamedTemporaryFile(delete=True).name
    other_intervals = other_intervals[(other_intervals.Chromosome == chrom) & (other_intervals.Strand == strand)]
    to_db(temp_path, "temp", other_intervals, backend=backend)
    other_intervals_merged = merge_intervals("temp", temp_path, [(chrom, strand)], feature_filter=feature_filter, backend=backend)
    # delete temp file
    os.remove(temp_path)
    
    # add __num__ column to self_genes_sql which counts how many intervals in self_genes_sql overlap with other_intervals
    self_intervals = pyranges.methods.coverage._number_overlapping(self_intervals, other_intervals_merged, strandedness="same", keep_nonoverlapping=True, overlap_col="__num__")
    subtracted_intervals = pyranges.methods.subtraction._subtraction(self_intervals, other_intervals_merged, strandedness="same")
    if subtracted_intervals is not None:
        return subtracted_intervals.drop(columns=["__num__"])
    else:
        return None

def subtract_intervals(sql_table_name: str, sql_db_name: str, chrom_strand_tup: list, other_intervals: pd.DataFrame, feature_filter: None | str = None, backend="duckdb") -> pd.DataFrame:
    """Subtract a set of other intervals from the database intervals for multiple chromosomes and strands. The function can also optionaly filter the intervals based on a specific feature.

    Args:
        sql_table_name (str): Name of the SQL table.
        sql_db_name (str): Name of the SQL database.
        chrom_strand_tup (list): List of tuples containing chromosome and strand information.
        other_intervals (pd.DataFrame): DataFrame containing the gene intervals to subtract from the database.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand' (and 'Feature' if feature_filter is set).
        feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
        backend (str, optional): Database backend to use. Defaults to "duckdb".

    Returns:
        pd.DataFrame: A DataFrame containing the subtracted intervals for all specified chromosomes and strands.
    """
    futures = [subtract_intervals_single.remote(sql_table_name, sql_db_name, chrom_strand, other_intervals, feature_filter=feature_filter, backend=backend) for chrom_strand in chrom_strand_tup]
    subtracted_intervals = ray.get(futures)
    subtracted_intervals = pd.concat(subtracted_intervals)
    return subtracted_intervals