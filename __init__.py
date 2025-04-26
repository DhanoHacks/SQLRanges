from . import duckdb_queries
from . import duckdb_utils
from . import sqlite3_queries
from . import sqlite3_utils
import os
import pandas as pd
import ray
import pyranges as pr

class SQLRanges:
    """
    This class provides methods to perform various operations on genomic data stored in a SQL database.
    
    It supports both SQLite and DuckDB backends. The input file must be in GTF or GFF3 format.
    """
    def __init__(self, input_file, table_name, db_name, backend='sqlite3', file_format="gtf"):
        """
        Initialize the SQLRanges class.

        Args:
            input_file (str): Path to the input file.
            table_name (str): Name of the table to be created in the database.
            db_name (str): Name of the database file (e.g., 'database.db').
            backend (str, optional): Database backend to use. Either 'sqlite3' or 'duckdb'. Defaults to 'sqlite3'.
            file_format (str, optional): Format of the input file. Either 'gtf' or 'gff3'. Defaults to 'gtf'.
        """
        assert backend in ['sqlite3', 'duckdb'], "Backend must be either 'sqlite3' or 'duckdb'"
        if backend == 'sqlite3':
            self.utils_backend = sqlite3_utils
            self.query_backend = sqlite3_queries
        else:
            self.utils_backend = duckdb_utils
            self.query_backend = duckdb_queries
        self.db_name = db_name
        self.table_name = table_name
        self.utils_backend.to_db(self.db_name, self.table_name, input_file, format=file_format)
        self.chrom_strand_tup = self.utils_backend.get_chrom_strand_tup(self.table_name, self.db_name)
        self.conn = self.query_backend.get_connection(self.db_name)
    
    def to_pyranges(self) -> pr.PyRanges:
        """Convert the SQL table to a PyRanges object.

        Returns:
            pyranges.PyRanges: A PyRanges object containing the genomic data from the SQL table.
        """
        return self.utils_backend.to_pyranges(self.conn, self.table_name)
        
    def count_exons(self) -> pd.DataFrame:
        """Count the number of exons for each gene in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the gene IDs and their corresponding exon counts.
        """
        return self.query_backend.count_exons(self.table_name, self.conn)
    
    def exon_length(self) -> pd.DataFrame:
        """Calculate the total length of exons for each gene in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the gene IDs and their corresponding total exon lengths.
        """
        return self.query_backend.exon_length(self.table_name, self.conn)
    
    def highest_transcripts(self) -> pd.DataFrame:
        """Identify the chromosome with the highest number of transcripts in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the chromosome name and its corresponding transcript count.
        """
        return self.query_backend.highest_transcripts(self.table_name, self.conn)
    
    def merge_exon_intervals(self) -> pd.DataFrame:
        """Merge overlapping exon intervals.

        Returns:
            pandas.DataFrame: A DataFrame containing the merged exon intervals.
        """
        return self.query_backend.merge_exon_intervals(self.table_name, self.db_name, self.chrom_strand_tup)
    
    def get_overlapping_genes(self, other_genes):
        """Get overlapping genes with the provided gene intervals.

        Args:
            other_genes (pandas.DataFrame): A DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.

        Returns:
            pandas.DataFrame: A DataFrame containing the overlapping genes.
        """
        return self.query_backend.get_overlapping_genes(self.table_name, self.db_name, self.chrom_strand_tup, other_genes)
    
    def get_subtracted_exons(self, other_cdf):
        """Remove a set of repetitive intervals from the exon features.

        Args:
            other_cdf (pyranges.PyRanges): A PyRanges object containing the intervals to be subtracted.
            The PyRanges object should have columns 'Chromosome', 'Start', 'End', and 'Strand'.

        Returns:
            pandas.DataFrame: A DataFrame containing the subtracted exons.
        """
        return self.query_backend.get_subtracted_exons(self.table_name, self.db_name, self.chrom_strand_tup, other_cdf)