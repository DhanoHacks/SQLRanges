from . import utils
from . import queries
import os
import pandas as pd
import ray
import pyranges as pr

class sqlranges:
    """
    This class provides methods to perform various operations on genomic data stored in a SQL database.
    
    It supports both SQLite and DuckDB backends. The input file must be in GTF or GFF3 format.
    """
    def __init__(self, input: str | pd.DataFrame, table_name: str, db_name: str, backend: str = "sqlite3", file_format: str = "gtf"):
        """
        Initialize the sqlranges class.

        Args:
            input (str | pandas.DataFrame): Path to the input file (GTF or GFF3 or SQLite3 or Duckdb) or a pandas DataFrame containing genomic data.
            table_name (str): Name of the table to be created in the database.
            db_name (str): Name of the database file (e.g., 'database.db').
            backend (str, optional): Database backend to use. Either 'sqlite3' or 'duckdb'. Defaults to 'sqlite3'.
            file_format (str, optional): Format of the input file. Either 'gtf' or 'gff3' or 'sqlite3' or 'duckdb'. Defaults to 'gtf'.
        """
        assert backend in ['sqlite3', 'duckdb'], "Backend must be either 'sqlite3' or 'duckdb'"
        assert file_format in ['gtf', 'gff3', 'sqlite3', 'duckdb'], "File format must be either 'gtf', 'gff3', 'sqlite3' or 'duckdb'"
        self.db_name = db_name
        self.table_name = table_name
        self.backend = backend
        if file_format in ['gtf', 'gff3']:
            utils.to_db(self.db_name, self.table_name, input, format=file_format, backend=self.backend)
        else:
            assert file_format == backend, "If file_format is 'sqlite3' or 'duckdb', the backend must be the same"
        self.chrom_strand_tup = utils.get_chrom_strand_tup(self.table_name, self.db_name, backend=self.backend)
        self.conn = queries.get_connection(self.db_name, backend=self.backend)
        
    def to_pandas(self) -> pd.DataFrame:
        """Convert the SQL table to a pandas DataFrame.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the genomic data from the SQL table.
        """
        return utils.to_pandas(self.conn, self.table_name, backend=self.backend)
    
    def to_pyranges(self) -> pr.PyRanges:
        """Convert the SQL table to a PyRanges object.

        Returns:
            pyranges.PyRanges: A PyRanges object containing the genomic data from the SQL table.
        """
        return utils.to_pyranges(self.conn, self.table_name, backend=self.backend)

    def to_gtf(self, output_path: str, comments: list[str] = []) -> None:
        """Export the SQL table to a GTF file.

        Args:
            output_path (str): Path to the output GTF file.
            comments (list[str], optional): List of comments to be added to the top of the GTF file. Defaults to [].
        """
        utils.to_gtf(self.conn, self.table_name, output_path, comments=comments, backend=self.backend)
        
    def to_gff3(self, output_path: str, comments: list[str] = []) -> None:
        """Export the SQL table to a GFF3 file.

        Args:
            output_path (str): Path to the output GFF3 file.
            comments (list[str], optional): List of comments to be added to the top of the GFF3 file. Defaults to [].
        """
        utils.to_gff3(self.conn, self.table_name, output_path, comments=comments, backend=self.backend)

    def query_sql(self, sql: str) -> pd.DataFrame:
        """Execute a SQL query on the database.

        Args:
            sql (str): The SQL query to be executed.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the results of the SQL query.
        """
        return queries.query_db(sql, self.conn, backend=self.backend)
        
    def count_intervals(self, group_by: str = "gene_id", feature_filter: None | str = None, return_col_name: str = "count") -> pd.DataFrame:
        """Count the number of intervals in the database, grouped by a specified column. The function can also optionally filter the intervals based on a specific feature.

        Args:
            group_by (str, optional): Column to group by. Defaults to "gene_id".
            feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
            return_col_name (str, optional): Column name for the count result. Defaults to "count".

        Returns:
            pd.DataFrame: A DataFrame containing the grouped counts.
        """
        return queries.count_intervals(self.table_name, self.conn, group_by=group_by, feature_filter=feature_filter, return_col_name=return_col_name, backend=self.backend)
    
    def total_length(self, group_by: str = "gene_id", feature_filter: None | str = None, return_col_name: str = "total_length") -> pd.DataFrame:
        """Calculate the total length of intervals in the database, grouped by a specified column. The function can also optionally filter the intervals based on a specific feature.

        Args:
            group_by (str, optional): Column to group by. Defaults to "gene_id".
            feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.
            return_col_name (str, optional): Column name for the total length result. Defaults to "total_length".

        Returns:
            pd.DataFrame: A DataFrame containing the grouped total lengths.
        """
        return queries.total_length(self.table_name, self.conn, group_by=group_by, feature_filter=feature_filter, return_col_name=return_col_name, backend=self.backend)
    
    def merge_intervals(self, feature_filter: None | str = None) -> pd.DataFrame:
        """Merge overlapping intervals in the database. The function can also optionally filter the intervals based on a specific feature.

        Args:
            feature_filter (None | str, optional): Filter for specific features. If None, no filter is applied. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the merged intervals.
        """
        return queries.merge_intervals(self.table_name, self.db_name, self.chrom_strand_tup, feature_filter=feature_filter, backend=self.backend)
    
    def overlapping_intervals(self, other_intervals: "sqlranges", feature_filter: None | str = None, other_feature_filter: None | str = None) -> pd.DataFrame:
        """Find overlapping intervals between the database and a set of other intervals. The function can also optionally filter the intervals based on a specific feature.

        Args:
            other_intervals (sqlranges): A sqlranges object containing the other intervals to find overlaps with.
                The other intervals should have the columns 'Chromosome', 'Start', 'End', and 'Strand'.
            feature_filter (None | str, optional): Filter for specific features on the database intervals. If None, no filter is applied. Defaults to None.
            other_feature_filter (None | str, optional): Filter for specific features on the other intervals. If None, no filter is applied. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the overlapping intervals.
        """
        return queries.overlapping_intervals(self.table_name, self.db_name, self.chrom_strand_tup, other_intervals.table_name, other_intervals.db_name, feature_filter=feature_filter, other_feature_filter=other_feature_filter, backend=self.backend, other_backend=other_intervals.backend)
    
    def subtract_intervals(self, other_intervals: "sqlranges", feature_filter: None | str = None, other_feature_filter: None | str = None) -> pd.DataFrame:
        """Subtract a set of other intervals from the database intervals. The function can also optionally filter the intervals based on a specific feature.

        Args:
            other_intervals (sqlranges): A sqlranges object containing the other intervals to subtract.
                The other intervals should have the columns 'Chromosome', 'Start', 'End', and 'Strand'.
            feature_filter (None | str, optional): Filter for specific features on the database intervals. If None, no filter is applied. Defaults to None.
            other_feature_filter (None | str, optional): Filter for specific features on the other intervals. If None, no filter is applied. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the intervals after subtraction.
        """
        return queries.subtract_intervals(self.table_name, self.db_name, self.chrom_strand_tup, other_intervals.table_name, other_intervals.db_name, feature_filter=feature_filter, other_feature_filter=other_feature_filter, backend=self.backend, other_backend=other_intervals.backend)
    
__all__ = ["sqlranges"]