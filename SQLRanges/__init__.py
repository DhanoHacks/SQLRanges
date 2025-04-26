from . import utils
from . import queries
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
        self.db_name = db_name
        self.table_name = table_name
        self.backend = backend
        utils.to_db(self.db_name, self.table_name, input_file, format=file_format, backend=self.backend)
        self.chrom_strand_tup = utils.get_chrom_strand_tup(self.table_name, self.db_name, backend=self.backend)
        self.conn = queries.get_connection(self.db_name, backend=self.backend)
    
    def to_pyranges(self) -> pr.PyRanges:
        """Convert the SQL table to a PyRanges object.

        Returns:
            pyranges.PyRanges: A PyRanges object containing the genomic data from the SQL table.
        """
        return utils.to_pyranges(self.conn, self.table_name, backend=self.backend)
        
    def count_exons(self) -> pd.DataFrame:
        """Count the number of exons for each gene in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the gene IDs and their corresponding exon counts.
        """
        return queries.count_exons(self.table_name, self.conn, backend=self.backend)
    
    def exon_length(self) -> pd.DataFrame:
        """Calculate the total length of exons for each gene in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the gene IDs and their corresponding total exon lengths.
        """
        return queries.exon_length(self.table_name, self.conn, backend=self.backend)
    
    def highest_transcripts(self) -> pd.DataFrame:
        """Identify the chromosome with the highest number of transcripts in the database.

        Returns:
            pandas.DataFrame: A DataFrame containing the chromosome name and its corresponding transcript count.
        """
        return queries.highest_transcripts(self.table_name, self.conn, backend=self.backend)
    
    def merge_exon_intervals(self) -> pd.DataFrame:
        """Merge overlapping exon intervals.

        Returns:
            pandas.DataFrame: A DataFrame containing the merged exon intervals.
        """
        return queries.merge_exon_intervals(self.table_name, self.db_name, self.chrom_strand_tup, backend=self.backend)
    
    def get_overlapping_genes(self, other_genes):
        """Get overlapping genes with the provided gene intervals.

        Args:
            other_genes (pandas.DataFrame): A DataFrame containing the gene intervals to check for overlaps with.
            The DataFrame should have columns 'Chromosome', 'Start', 'End', and 'Strand'.

        Returns:
            pandas.DataFrame: A DataFrame containing the overlapping genes.
        """
        return queries.get_overlapping_genes(self.table_name, self.db_name, self.chrom_strand_tup, other_genes, backend=self.backend)
    
    def get_subtracted_exons(self, other_cdf):
        """Remove a set of repetitive intervals from the exon features.

        Args:
            other_cdf (pyranges.PyRanges): A PyRanges object containing the intervals to be subtracted.
            The PyRanges object should have columns 'Chromosome', 'Start', 'End', and 'Strand'.

        Returns:
            pandas.DataFrame: A DataFrame containing the subtracted exons.
        """
        return queries.get_subtracted_exons(self.table_name, self.db_name, self.chrom_strand_tup, other_cdf, backend=self.backend)
    
__all__ = ["SQLRanges"]