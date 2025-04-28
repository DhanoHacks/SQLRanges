import sqlite3
import duckdb
import time
import ray
import pandas as pd
import pyranges as pr

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

def query_db(query: str, conn: sqlite3.Connection | duckdb.DuckDBPyConnection, backend: str = "duckdb", return_df=True) -> pd.DataFrame | None:
    """Execute a SQL query and return the result as a pandas DataFrame (or nothing).

    Args:
        query (str): SQL query to execute.
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        backend (str, optional): Database backend to use. Defaults to "duckdb".
        return_df (bool, optional): Whether to return the result (as a DataFrame) or not. Defaults to True.

    Returns:
        pd.DataFrame | None: Result of the query as a pandas DataFrame, or None if return_df is False.
    """
    if return_df:
        if backend == "duckdb":
            return conn.execute(query).fetchdf()
        else:
            return pd.read_sql_query(query, conn)
    else:
        return conn.execute(query)

def get_chrom_strand_tup(sql_table_name: str, db_name: str, backend: str = "duckdb") -> list:
    """Get a list of unique Chromosome and Strand tuples from the database.

    Args:
        sql_table_name (str): Name of the SQL table to query.
        db_name (str): Name of the database file.
        backend (str, optional): Database backend to use, either "sqlite3" or "duckdb". Defaults to "duckdb".

    Returns:
        list: A list of tuples containing unique Chromosome and Strand values.
    """
    conn = get_connection(db_name, backend=backend)
    chrom_strand_tup = query_db(f"SELECT DISTINCT Chromosome, Strand FROM \"{sql_table_name}\"", conn, backend=backend)
    conn.close()
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    return chrom_strand_tup
    
def process_line(line: str, format: str = "gtf") -> dict:
    """
    Process a single line of GTF or GFF3 file and return a dictionary of attributes.

    Args:
        line (str): A line from the GTF or GFF3 file.
        format (str, optional): The format of the file, either "gtf" or "gff3". Defaults to "gtf".

    Raises:
        ValueError: If the line has less than 9 tab-separated fields.
        ValueError: If there is an error processing the 'Start' or 'End' value.

    Returns:
        dict: A dictionary containing the attributes of the line.
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
        raise ValueError(f"Error processing 'Start' or 'End' value: {data['Start']}, {data['End']}") from e
    
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
def process_batch(lines_batch: list, format: str = "gtf") -> pd.DataFrame:
    """
    Process a batch of lines from a GTF or GFF3 file and return a DataFrame.
    This function is designed to be run in parallel using Ray.

    Args:
        lines_batch (list): A batch of lines from the GTF or GFF3 file.
        format (str): The format of the file, either "gtf" or "gff3". Defaults to "gtf".

    Returns:
        pd.DataFrame: A DataFrame containing the processed lines.
    """
    # Process each line in the batch
    line_dicts = []
    for line in lines_batch:
        line_dict = process_line(line, format)
        if line_dict:
            line_dicts.append(line_dict)
    return pd.DataFrame(line_dicts)

# def to_db(sql_db_name, sql_table_name, input_file, chunk_size=4000000, format="gtf", backend="duckdb"):
def to_db(sql_db_name: str, sql_table_name: str, input: str | pd.DataFrame, chunk_size: int = 4000000, format: str = "gtf", backend: str = "duckdb"):
    """
    Convert a GTF or GFF3 file to a SQL database table.

    Args:
        sql_db_name (str): Name of the SQL database file (e.g., 'database.db').
        sql_table_name (str): Name of the SQL table to create.
        input (str | pandas.DataFrame): Path to the input file (GTF or GFF3) or a pandas DataFrame containing genomic data.
        chunk_size (int, optional): Size of the chunks (in bytes) to read from the file. Defaults to 4000000 bytes.
        format (str, optional): Format of the input file, either "gtf" or "gff3". Defaults to "gtf".
        backend (str, optional): Database backend to use, either "sqlite3" or "duckdb". Defaults to "duckdb".
    """
    assert format in ["gtf", "gff3"], "Format must be either 'gtf' or 'gff3'."
    assert backend in ["sqlite3", "duckdb"], "Backend must be either 'sqlite3' or 'duckdb'."
    
    # If input is a DataFrame, convert it to a list of lines
    if isinstance(input, pd.DataFrame):
        if backend == "duckdb":
            conn = duckdb.connect(sql_db_name)
            conn.execute(f"DROP TABLE IF EXISTS {sql_table_name}")
            conn.execute(f"CREATE TABLE {sql_table_name} AS SELECT * FROM input")
            conn.commit()
            conn.close()
        else:
            conn = sqlite3.connect(sql_db_name)
            conn.execute(f"DROP TABLE IF EXISTS {sql_table_name}")
            input.to_sql(sql_table_name, conn, if_exists="replace", index=False)
            conn.commit()
            conn.close()
        return
    
    # if ray is not initialized, initialize it
    if not ray.is_initialized():
        ray.init()
    
    # If input is a file, read it in chunks
    with open(input, "r") as f:
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
    
    conn = get_connection(sql_db_name, backend=backend)
    query_db(f"DROP TABLE IF EXISTS \"{sql_table_name}\"", conn, backend=backend, return_df=False)
    if backend == "duckdb":
        conn.execute(f"CREATE TABLE \"{sql_table_name}\" AS SELECT * FROM df")
    else:
        df.to_sql(sql_table_name, conn, if_exists="replace", index=False)
    conn.commit()
    conn.close()

def to_pandas(conn: sqlite3.Connection | duckdb.DuckDBPyConnection, table_name: str, backend: str = "duckdb") -> pd.DataFrame:
    """Convert a SQL table to a pandas DataFrame.

    Args:
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        table_name (str): Name of the SQL table to query.
        backend (str, optional): Database backend to use, either "sqlite3" or "duckdb". Defaults to "duckdb".

    Returns:
        pd.DataFrame: A pandas DataFrame containing the genomic data from the SQL table.
    """
    return query_db(f"SELECT * FROM \"{table_name}\"", conn, backend=backend)

def to_pyranges(conn: sqlite3.Connection | duckdb.DuckDBPyConnection, table_name: str, backend: str = "duckdb") -> pr.PyRanges:
    """
    Convert a SQL table to a PyRanges object.

    Args:
        conn (sqlite3.Connection | duckdb.DuckDBPyConnection): Database connection object.
        table_name (str): Name of the SQL table to query.
        backend (str, optional): Database backend to use, either "sqlite3" or "duckdb". Defaults to "duckdb".

    Returns:
        pr.PyRanges: A PyRanges object containing the genomic data from the SQL table.
    """
    return pr.PyRanges(to_pandas(conn, table_name, backend=backend))