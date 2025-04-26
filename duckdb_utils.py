import duckdb
import time
import ray
import pandas as pd
import pyranges as pr

def get_chrom_strand_tup(sql_table_name, db_name):
    conn = duckdb.connect(db_name, read_only=True)
    chrom_strand_tup = conn.execute(f"SELECT DISTINCT Chromosome, Strand FROM {sql_table_name}").fetchdf()
    conn.close()
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    return chrom_strand_tup
    
def process_line(line, format):
    """
    Process a single GTF file line and return an SQL INSERT statement.
    
    The function splits the line on tab characters,
    processes the first 8 fields and subtracts 1 from the Start field.
    Then it parses the attribute field (the ninth column) by splitting on ';'
    and further splitting each attribute on the first space.
    It finally builds an INSERT statement inserting the collected key/value pairs.
    """
    # dont process comment lines i.e lines starting with '#'
    if line.startswith('#'):
        return None
    line = line.rstrip("\n")
    # Strip the trailing newline and split by tab:
    fields = line.strip().split('\t')
    if len(fields) < 9:
        raise ValueError("Error: Line has less than 9 tab-separated fields.")
    
    # Prepare the first columns in order:
    keys_list = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"]
    data = {}
    for i, key in enumerate(keys_list):
        data[key] = fields[i]
    
    # Adjust the Start field (subtract 1)
    try:
        data["Start"] = int(data["Start"]) - 1
        data["End"] = int(data["End"])
    except Exception as e:
        raise ValueError(f"Error processing 'Start' value: {data['Start']}") from e
    
    # Process the attribute field (9th field)
    attributes = fields[8]
    for attribute in attributes.split(';'):
        attribute = attribute.lstrip()  # Remove leading whitespace
        if not attribute:
            continue  # Skip empty tokens
        # Split attribute into key and value at first space
        parts = attribute.split(' ' if format == "gtf" else "=", 1)
        if len(parts) == 2:
            key, value = parts
            # Remove any quotation marks and extra whitespace from the value
            value = value.replace('"', '').strip()
            data[key] = value
    return data

@ray.remote
def process_batch(lines_batch, format):
    # Process each line in the batch
    line_dicts = []
    for line in lines_batch:
        # You can call your original process_line function here.
        # Or, directly implement the processing logic.
        line_dict = process_line(line, format)
        if line_dict:
            line_dicts.append(line_dict)
    return pd.DataFrame(line_dicts)

def to_db(duckdb_db_name, duckdb_table_name, input_file, chunk_size=4000000, format="gtf"):
    """
    Main function to process a GTF file and write its data into an CSV file.
    
    Parameters:
        csv_name    : Name (path) of the CSV file to write to.
        input_file : GTF file to process.
        batch_size : Number of bytes/characters to read from the file at once. Default is 1000000.
        chunk_size : Number of lines to process in each batch. Default is 8000.
    """
    assert format in ["gtf", "gff3"], "Format must be either 'gtf' or 'gff3'."

    with open(input_file, "r") as f:
        # Process the lines in batches using Ray
        start_time = time.time()
        futures = []
        while True:
            # Each batch processes chunk_size bytes at once.
            lines_batch = f.readlines(chunk_size)
            if not lines_batch:
                break
            futures.append(process_batch.remote(lines_batch, format))
        elapsed = time.time() - start_time
        print(f"Submitted {len(futures)} tasks to ray in {elapsed*1000:.0f}ms")
        
        # Wait for all futures to complete and collect the results.
        start_time = time.time()
        dfs = ray.get(futures)
        elapsed = time.time() - start_time
        print(f"Processed all lines in {elapsed*1000:.0f}ms")
    # Concatenate all line_dicts into a single DataFrame
    start_time = time.time()
    df = pd.concat(dfs, ignore_index=True)
    elapsed = time.time() - start_time
    print(f"Created DataFrame with {len(df)} rows in {elapsed*1000:.0f}ms")
    
    duckdb_conn = duckdb.connect(duckdb_db_name)
    duckdb_conn.execute(f"DROP TABLE IF EXISTS {duckdb_table_name}")
    duckdb_conn.execute(f"CREATE TABLE {duckdb_table_name} AS SELECT * FROM df")
    duckdb_conn.commit()
    duckdb_conn.close()

def to_pyranges(conn, table_name):
    return pr.PyRanges(conn.execute(f"SELECT * FROM {table_name}").fetchdf())