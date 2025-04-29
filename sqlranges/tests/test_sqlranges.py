import pytest
import pandas as pd
from sqlranges import sqlranges
import gzip
import shutil
import urllib.request
import os

GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz"

@pytest.fixture(scope="session")
def downloaded_gtf(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("data")
    gz_path = tmp_dir / "gencode.vM36.annotation.gtf.gz"
    gtf_path = tmp_dir / "gencode.vM36.annotation.gtf"

    if not gz_path.exists():
        urllib.request.urlretrieve(GTF_URL, gz_path)

    with gzip.open(gz_path, "rb") as f_in:
        with open(gtf_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return gtf_path

@pytest.fixture(params=["duckdb", "sqlite3"])
def sqlr(downloaded_gtf, tmp_path, request):
    db_path = tmp_path / f"test-{request.param}.db"
    return sqlranges(str(downloaded_gtf), table_name="mouse", db_name=str(db_path), backend=request.param, file_format="gtf")

def test_count_exons(sqlr):
    df = sqlr.count_intervals(group_by="gene_id", feature_filter="exon", return_col_name="exon_count")
    df["exon_count"] = df["exon_count"].astype(int)
    expected = pd.read_csv("expected_outputs/exon_counts_mouse.csv")
    assert df.equals(expected)

def test_total_exon_length(sqlr):
    df = sqlr.total_length(group_by="gene_id", feature_filter="exon", return_col_name="total_exon_length")
    df["total_exon_length"] = df["total_exon_length"].astype(int)
    expected = pd.read_csv("expected_outputs/total_exon_length_mouse.csv")
    assert df.equals(expected)

def test_highest_transcripts(sqlr):
    df = sqlr.count_intervals(group_by="Chromosome", feature_filter="transcript")
    chrom = df["Chromosome"].iloc[df["count"].idxmax()]
    assert chrom == "chr2"

def test_merge_exon_intervals(sqlr):
    df = sqlr.merge_intervals(feature_filter="exon").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    expected = pd.read_csv("expected_outputs/merged_exon_intervals_mouse.csv")
    assert df.equals(expected)

def test_overlapping_genes(sqlr):
    other = pd.DataFrame({"Chromosome": ["chr1"], "Start": [3000000], "End": [4000000], "Strand": ["+"], "Feature": ["gene"]})
    df = sqlr.overlapping_intervals(other, feature_filter="gene").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    expected = pd.read_csv("expected_outputs/overlapping_genes_mouse.csv")
    assert df.equals(expected)

def test_subtracted_exons(sqlr, tmp_path):
    other_genes = sqlr.query_sql("""
        SELECT Chromosome, MIN(Start)+1000000 AS Start, MAX("End")-1000000 AS "End", 
               Strand, ANY_VALUE(Feature) AS Feature 
        FROM mouse WHERE Feature = 'exon' 
        GROUP BY Chromosome, Strand
    """)
    out_path = tmp_path / "subtracted.csv"
    sqlr.subtract_intervals(other_genes, feature_filter="exon").to_csv(out_path, index=False)

    # Sort both files for comparison
    sorted_test = tmp_path / "sorted_test.csv"
    sorted_expected = tmp_path / "sorted_expected.csv"
    os.system(f"sort {out_path} > {sorted_test}")
    os.system(f"sort expected_outputs/subtracted_exons_mouse.csv > {sorted_expected}")

    diff = os.popen(f"diff {sorted_test} {sorted_expected}").read()
    assert diff == "", "Subtracted exons do not match expected output."
