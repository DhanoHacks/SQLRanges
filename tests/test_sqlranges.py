import pytest
import pandas as pd
from sqlranges import sqlranges
import gzip
import shutil
import urllib.request
import os
import warnings
import pathlib
warnings.filterwarnings("ignore")

GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz"

@pytest.fixture(scope="session")
def downloaded_gtf(tmp_path_factory: pytest.TempPathFactory):
    """Downloads the GTF file if it doesn't exist and returns its path.
    This fixture is used to ensure that the GTF file is available for testing.

    Args:
        tmp_path_factory (pytest.TempPathFactory): Fixture for creating temporary paths.
    """
    tmp_dir = tmp_path_factory.mktemp("data")
    gz_path = tmp_dir / "gencode.vM36.annotation.gtf.gz"
    gtf_path = tmp_dir / "gencode.vM36.annotation.gtf"

    if not gz_path.exists():
        urllib.request.urlretrieve(GTF_URL, gz_path)

    with gzip.open(gz_path, "rb") as f_in:
        with open(gtf_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

@pytest.fixture(params=["duckdb", "sqlite3"])
def sqlr(downloaded_gtf: str, tmp_path: pathlib.Path, request: pytest.FixtureRequest):
    """Fixture to create an instance of sqlranges for testing.
    This fixture sets up the database connection and provides a clean instance
    of sqlranges for each test. It uses the downloaded GTF file and the specified
    database backend (either duckdb or sqlite3).
    The database is created in a temporary path for each test run.

    Args:
        downloaded_gtf (str): Path to the downloaded GTF file.
        tmp_path (pathlib.Path): Temporary path for creating the database.
        request (pytest.FixtureRequest): Fixture request object to access the parameterized value.

    Returns:
        sqlranges: An instance of sqlranges connected to the specified database.
    """
    db_path = tmp_path / f"test-{request.param}.db"
    return sqlranges(str(downloaded_gtf), table_name="mouse", db_name=str(db_path), backend=request.param, file_format="gtf")

def test_count_exons(sqlr: sqlranges):
    """Test the count of exons in the database.
    This test checks the number of exons in the database by grouping
    by gene_id and filtering for the "exon" feature. The result is compared
    to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
    """
    df = sqlr.count_intervals(group_by="gene_id", feature_filter="exon", return_col_name="exon_count")
    df["exon_count"] = df["exon_count"].astype(int)
    df = df.sort_values("gene_id").reset_index(drop=True)
    expected = pd.read_csv("tests/expected_outputs/exon_counts_mouse.csv")
    assert df.equals(expected)

def test_total_exon_length(sqlr: sqlranges):
    """Test the total length of exons in the database.
    This test calculates the total length of exons by grouping
    by gene_id and filtering for the "exon" feature. The result is compared
    to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
    """
    df = sqlr.total_length(group_by="gene_id", feature_filter="exon", return_col_name="total_exon_length")
    df["total_exon_length"] = df["total_exon_length"].astype(int)
    df = df.sort_values("gene_id").reset_index(drop=True)
    expected = pd.read_csv("tests/expected_outputs/total_exon_length_mouse.csv")
    assert df.equals(expected)

def test_highest_transcripts(sqlr : sqlranges):
    """Test the highest number of transcripts in a gene.
    This test counts the number of transcripts in each gene by grouping
    by gene_id and filtering for the "transcript" feature. The result is compared
    to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
    """
    df = sqlr.count_intervals(group_by="Chromosome", feature_filter="transcript")
    chrom = df["Chromosome"].iloc[df["count"].idxmax()]
    assert chrom == "chr2"

def test_merge_exon_intervals(sqlr : sqlranges):
    """Test the merging of overlapping exon intervals.
    This test merges overlapping exon intervals in the database by filtering
    for the "exon" feature. The result is compared to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
    """
    df = sqlr.merge_intervals(feature_filter="exon").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    expected = pd.read_csv("tests/expected_outputs/merged_exon_intervals_mouse.csv")
    assert df.equals(expected)

def test_overlapping_genes(sqlr: sqlranges, tmp_path: pathlib.Path):
    """Test the overlapping genes in the database.
    This test checks for overlapping genes in the database by filtering
    for the "gene" feature. The result is compared to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
        tmp_path (pathlib.Path): Temporary path for creating the output file.
    """
    other = pd.DataFrame({"Chromosome": ["chr1"], "Start": [3000000], "End": [4000000], "Strand": ["+"], "Feature": ["gene"]})
    sqlr.overlapping_intervals(other, feature_filter="gene").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True).to_csv(tmp_path / "overlapping_genes.csv", index=False)
    diff = os.popen(f"diff {tmp_path / 'overlapping_genes.csv'} tests/expected_outputs/overlapping_genes_mouse.csv").read()
    assert diff == "", "Overlapping genes do not match expected output."

def test_subtracted_exons(sqlr: sqlranges, tmp_path: pathlib.Path):
    """Test the subtraction of exon intervals from the database.
    This test subtracts exon intervals from the database by filtering
    for the "exon" feature. The result is compared to an expected output file.

    Args:
        sqlr (sqlranges): An instance of sqlranges connected to the database.
        tmp_path (pathlib.Path): Temporary path for creating the output file.
    """
    other_genes = sqlr.query_sql("""
        SELECT Chromosome, MIN(Start)+1000000 AS Start, MAX("End")-1000000 AS "End", 
               Strand, 'exon' AS Feature 
        FROM mouse WHERE Feature = 'exon' 
        GROUP BY Chromosome, Strand
    """)
    out_path = tmp_path / "subtracted.csv"
    sqlr.subtract_intervals(other_genes, feature_filter="exon").to_csv(out_path, index=False)

    # Sort both files for comparison
    sorted_test = tmp_path / "sorted_test.csv"
    sorted_expected = tmp_path / "sorted_expected.csv"
    os.system(f"sort {out_path} > {sorted_test}")
    os.system(f"sort tests/expected_outputs/subtracted_exons_mouse.csv > {sorted_expected}")

    diff = os.popen(f"diff {sorted_test} {sorted_expected}").read()
    assert diff == "", "Subtracted exons do not match expected output."
