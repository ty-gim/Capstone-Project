import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml

def build_cms_query(cms_folder, cms_files, targets, controls, drug_mappings):
    file_paths = [str(Path(cms_folder) / f) for f in cms_files]
    queries = []

    for fp in file_paths:
        year_groups = re.findall(r'(?<!\d)\d{4}(?!\d)', Path(fp).stem)
        year = year_groups[-1] if year_groups else 'Unknown'

        try:
            with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
                header_line = f.readline().strip().replace('\ufeff', '').replace('"', '').lower()
                header = header_line.split(',')
        except Exception as e:
            print(f"Error reading header for {fp}: {e}")
            continue

        search_cols = [c for c in header if c in ['gnrc_name', 'generic_name', 'brnd_name', 'brand_name', 'drug_name']]
        clms_col = next((c for c in header if c in ['tot_clms', 'total_claim_count', 'tot_clm_cnt']), None)

        if not search_cols or not clms_col:
            print(f"Warning: Missing required columns in {fp}.")
            continue

        case_when = ""
        where_conditions = []

        # Apply drug mappings to search for Brand AND Generic names
        for d in targets + controls:
            terms = drug_mappings.get(d, [d])
            # Create an OR condition for every column and every mapped term
            drug_conds = " OR ".join([f"{c} ILIKE '%{t}%'" for c in search_cols for t in terms])
            case_when += f"        WHEN {drug_conds} THEN '{d}'\n"
            where_conditions.append(f"({drug_conds})")

        where_clause = " OR ".join(where_conditions)

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

    union_query = " UNION ALL ".join(queries)

    sql_query = f"""
    SELECT generic_name, year, SUM(total_prescriptions) AS total_prescriptions
    FROM ({union_query})
    GROUP BY generic_name, year
    ORDER BY generic_name, year;
    """
    return sql_query

def main():
    parser = argparse.ArgumentParser(description="CMS Volume Trend Pipeline")
    parser.add_argument('--cms_folder', type=str, required=True)
    parser.add_argument('--cms_files', type=str, nargs='+', required=True)
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)

    args = parser.parse_args()
    print(f"Running CMS Volume Analysis for {args.targets} vs {args.controls}...")

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})

    con = duckdb.connect(database=':memory:')
    query = build_cms_query(args.cms_folder, args.cms_files, args.targets, args.controls, drug_mappings)

    result_df = con.execute(query).df()
    result_df['role'] = result_df['generic_name'].map(role_map)
    result_df.to_csv(args.output, index=False)
    print(f"CMS Volume data successfully saved to {args.output}")

if __name__ == "__main__":
    main()
