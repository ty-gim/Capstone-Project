import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml

def build_cms_query(cms_folder, cms_files, targets, controls, drug_mappings):
    """
    Builds an SQL query to extract and aggregate CMS Medicare Part D prescription data.

    Parameters:
    cms_folder (str): Path to the folder containing CMS files.
    cms_files (list): List of CMS file names to process.
    targets (list): List of target drug names.
    controls (list): List of control drug names.
    drug_mappings (dict): Dictionary mapping drug names to their brand/generic names.

    Returns:
    sql_query (str): SQL query that extracts/aggregates prescription data for the specified drugs.
    """
    file_paths = [str(Path(cms_folder) / f) for f in cms_files] # Convert file names to full paths
    queries = [] # Empty list to hold one SQL query for each file

    for fp in file_paths:
        # Look for a four-digit year pattern in the file name to extract the year
        year_groups = re.findall(r'(?<!\d)\d{4}(?!\d)', Path(fp).stem)
        year = year_groups[-1] if year_groups else 'Unknown'

        try:
            # Open the file with utf-8 encoding and ignore errors just to read the header line
            with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
                header_line = f.readline().strip().replace('\ufeff', '').replace('"', '').lower()
                header = header_line.split(',') # Split the header line into a list of columns
        except Exception as e:
            print(f"Error reading header for {fp}: {e}")
            continue

        # Search for columns that contain drug names and total claim counts
        # Government files often update formatting and naming conventions differently across years
        search_cols = [c for c in header if c in ['gnrc_name', 'generic_name', 'brnd_name', 'brand_name', 'drug_name']]
        clms_col = next((c for c in header if c in ['tot_clms', 'total_claim_count', 'tot_clm_cnt']), None)

        if not search_cols or not clms_col:
            print(f"Warning: Missing required columns in {fp}.")
            continue

        case_when = ""
        where_conditions = []

        # Identify the generic name column to check for combination products (e.g., Sitagliptin/Metformin)
        gnrc_col = next((c for c in header if c in ['gnrc_name', 'generic_name']), None)

        # Build a regex exclusion string: If it contains '/', '+', '&', or 'AND', it is a combination drug and must be excluded
        no_combo = f"AND NOT REGEXP_MATCHES(COALESCE({gnrc_col}, ''), '(?i)[/+&]|\\bAND\\b')" if gnrc_col else ""

        # Apply drug mappings to search for Brand AND Generic names using strict whole-word regex
        for d in targets + controls:
            terms = drug_mappings.get(d, [d])
            # \b ensures word boundaries. (?i) makes it case-insensitive.
            regex = f"(?i)\\b({'|'.join(re.escape(t) for t in terms)})\\b"

            # Create an OR condition for every column, explicitly excluding combination products
            drug_conds = "(" + " OR ".join([f"REGEXP_MATCHES(COALESCE({c}, ''), '{regex}')" for c in search_cols]) + f") {no_combo}"

            case_when += f"        WHEN {drug_conds} THEN '{d}'\n" # Map matched rows to the generic cohort name
            where_conditions.append(f"({drug_conds})") # Add to the WHERE clause to filter the dataset

        where_clause = " OR ".join(where_conditions) # Combine all WHERE conditions into a single OR clause

        # Build the query for this specific year's file
        q = f"""
        SELECT
            CASE
                {case_when}
            END AS generic_name,
            '{year}' AS year,
            TRY_CAST(REPLACE({clms_col}, ',', '') AS BIGINT) AS total_prescriptions
        FROM read_csv('{fp}', header=True, all_varchar=True, ignore_errors=True)
        WHERE ({where_clause})
        """
        queries.append(q)

    if not queries:
        raise SystemExit("No usable CMS files found (all skipped). Check the warnings above.")

    union_query = " UNION ALL ".join(queries) # Stack all yearly queries together

    # Final aggregation: Sum the total prescriptions grouped by drug and year
    sql_query = f"""
    SELECT generic_name, year, SUM(total_prescriptions) AS total_prescriptions
    FROM ({union_query})
    GROUP BY generic_name, year
    ORDER BY generic_name, year;
    """
    return sql_query

def main():
    """
    Parser for command line arguments and main execution function.
    """
    parser = argparse.ArgumentParser(description="CMS Volume Trend Pipeline")
    parser.add_argument('--cms_folder', type=str, required=True)
    parser.add_argument('--cms_files', type=str, nargs='+', required=True)
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)

    args = parser.parse_args()
    print(f"Running CMS Volume Analysis for {args.targets} vs {args.controls}...")

    with open(args.config_file, 'r') as f: # Snakemake passes the config file to this script
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    # Map target and control roles for downstream visualization
    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})

    con = duckdb.connect(database=':memory:') # Create an in-memory DuckDB connection
    query = build_cms_query(args.cms_folder, args.cms_files, args.targets, args.controls, drug_mappings)

    result_df = con.execute(query).df() # Execute the query and convert to a pandas DataFrame
    result_df['role'] = result_df['generic_name'].map(role_map)
    result_df.to_csv(args.output, index=False)
    print(f"CMS Volume data successfully saved to {args.output}")

if __name__ == "__main__":
    main()
