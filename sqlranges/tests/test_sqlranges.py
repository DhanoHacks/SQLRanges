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

    # Download the file
    urllib.request.urlretrieve(GTF_URL, gz_path)

    # Decompress
    with gzip.open(gz_path, "rb") as f_in:
        with open(gtf_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return gtf_path

def run_all_tests(downloaded_gtf, tmp_path):
    def test_count_exons(sqlr: sqlranges):
        df = sqlr.count_intervals(group_by="gene_id", feature_filter="exon", return_col_name="exon_count")
        df["exon_count"] = df["exon_count"].astype(int)
        expected_output = pd.read_csv("expected_outputs/exon_counts_mouse.csv")
        assert df.equals(expected_output), "Exon counts do not match expected output."
    
    def test_total_exon_length(sqlr: sqlranges):
        df = sqlr.total_length(group_by="gene_id", feature_filter="exon", return_col_name="total_exon_length")
        df["total_exon_length"] = df["total_exon_length"].astype(int)
        expected_output = pd.read_csv("expected_outputs/total_exon_length_mouse.csv")
        assert df.equals(expected_output), "Total exon lengths do not match expected output."
        
    def test_highest_transcripts(sqlr: sqlranges):
        df = sqlr.count_intervals(group_by="Chromosome", feature_filter="transcript")
        chrom = df["Chromosome"].iloc[df["count"].idxmax()]
        expected_chrom = "chr2"
        assert chrom == expected_chrom, f"Expected highest transcript chromosome to be {expected_chrom}, but got {chrom}."
    
    def test_merge_exon_intervals(sqlr: sqlranges):
        df = sqlr.merge_intervals(feature_filter="exon").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
        expected_output = pd.read_csv("expected_outputs/merged_exon_intervals_mouse.csv")
        assert df.equals(expected_output), "Merged exon intervals do not match expected output."

    def test_overlapping_genes(sqlr: sqlranges):
        other_genes = pd.DataFrame({"Chromosome": ["chr1"], "Start": [3000000], "End": [4000000], "Strand": ["+"], "Feature": ["gene"]})
        df = sqlr.overlapping_intervals(other_genes, feature_filter="gene").sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
        expected_output = pd.read_csv("expected_outputs/overlapping_genes_mouse.csv")
        assert df.equals(expected_output), "Overlapping genes do not match expected output."
        
    def test_subtracted_exons(sqlr: sqlranges, tmp_path: str):
        other_genes = sqlr.query_sql("SELECT Chromosome, MIN(Start)+1000000 AS Start, MAX(\"End\")-1000000 AS \"End\", Strand, ANY_VALUE(Feature) AS Feature FROM mouse WHERE Feature = 'exon' GROUP BY Chromosome, Strand")
        df = sqlr.subtract_intervals(other_genes, feature_filter="exon").to_csv(tmp_path / "subtracted_exons_sqlr.csv", index=False)
        # cat subtracted_exons_sqlr.csv, sort, and save to subtracted_exons_sqlr_sorted.csv
        os.system(f"cat {tmp_path}/subtracted_exons_sqlr.csv | sort > {tmp_path}/subtracted_exons_sorted.csv")
        os.system(f"cat expected_outputs/subtracted_exons_mouse.csv | sort > {tmp_path}/subtracted_exons_expected_sorted.csv")
        # Compare the two sorted files
        diff = os.popen(f"diff {tmp_path}/subtracted_exons_sorted.csv {tmp_path}/subtracted_exons_expected_sorted.csv").read()
        assert diff == "", "Subtracted exons do not match expected output."
    
    db_path = tmp_path / "db-mouse.duckdb"
    sqlr = sqlranges(str(downloaded_gtf), table_name="mouse", db_name=str(db_path), backend="duckdb", file_format="gtf")
    test_count_exons(sqlr)
    test_total_exon_length(sqlr)
    test_highest_transcripts(sqlr)
    test_merge_exon_intervals(sqlr)
    test_overlapping_genes(sqlr)
    test_subtracted_exons(sqlr, tmp_path)
    
    db_path = tmp_path / "db-mouse.sqlite3"
    sqlr = sqlranges(str(downloaded_gtf), table_name="mouse", db_name=str(db_path), backend="sqlite3", file_format="gtf")
    test_count_exons(sqlr)
    test_total_exon_length(sqlr)
    test_highest_transcripts(sqlr)
    test_merge_exon_intervals(sqlr)
    test_overlapping_genes(sqlr)
    test_subtracted_exons(sqlr, tmp_path)

