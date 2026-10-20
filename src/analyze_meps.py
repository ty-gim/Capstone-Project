import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml
import sys

def build_regex_conditions(drugs, drug_mappings):
    """Builds strict RE2-compatible regex to prevent combination-drug contamination."""
    conditions = []
    for d in drugs:
        terms = drug_mappings.get(d, [d])
        terms_regex = "|".join([re.escape(t) for t in terms])
        regex = f"(?i)\\b({terms_regex})\\b"
        conditions.append((d, regex))
    return conditions

def process_meps_files(meps_folder, meps_files, targets, controls, drug_mappings, con):
    file_paths = [str(Path(meps_folder) / f) for f in meps_files]

    all_drugs = targets + controls
    regex_mapping = build_regex_conditions(all_drugs, drug_mappings)
    combo_exclusion = r"[/+&]|\bAND\b"

    case_when_clauses = []
    drug_filter_conds = []

    for drug_name, regex in regex_mapping:
        pos_match = f"REGEXP_MATCHES(drug_name, '{regex}')"
        neg_match = f"REGEXP_MATCHES(drug_name, '(?i){combo_exclusion}')"
        full_cond = f"({pos_match} AND NOT {neg_match})"
        case_when_clauses.append(f"WHEN {full_cond} THEN '{drug_name}'")
        drug_filter_conds.append(full_cond)

    drug_case_sql = "\n".join(case_when_clauses)
    drug_filter_sql = " OR ".join(drug_filter_conds)

    for fp in file_paths:
        print(f"Loading and sanitizing {fp}...")
        if not Path(fp).exists():
            raise FileNotFoundError(f"MEPS file not found: {fp}")

        try:
            # Read Stata file in chunks/iterator if possible, but pandas read_stata loads to RAM.
            # MEPS files are usually small enough for Pandas RAM.
            df = pd.read_stata(fp)
            df.columns = [c.lower() for c in df.columns]
        except Exception as e:
            raise RuntimeError(f"Failed to read Stata file {fp}: {e}")

        # 1. Identify required columns dynamically
        drug_col = 'rxdrgnam' if 'rxdrgnam' in df.columns else 'rxname' if 'rxname' in df.columns else None
        total_exp_col = next((c for c in df.columns if c.startswith('rxxp') and c.endswith('x')), None)
        weight_col = next((c for c in df.columns if c.startswith('perwt') and c.endswith('f')), None)

        # MEPS Survey Design Variables
        varstr_col = 'varstr' if 'varstr' in df.columns else None
        varpsu_col = 'varpsu' if 'varpsu' in df.columns else None

        # 2. Strict Validation (No silent fallbacks)
        if not drug_col or not total_exp_col:
            raise ValueError(f"Missing required drug or expenditure columns in {fp}.")
        if not weight_col:
            raise ValueError(f"CRITICAL: Missing survey weight column in {fp}. Cannot compute national estimates.")
        if not varstr_col or not varpsu_col:
            raise ValueError(f"CRITICAL: Missing survey design variables (VARSTR/VARPSU) in {fp}.")

        year_match = re.search(r'rxxp([0-9]{2})x', total_exp_col)
        year = f"20{year_match.group(1)}" if year_match else "Unknown"

        # 3. Sanitize and register to DuckDB
        df_subset = df[[drug_col, total_exp_col, weight_col, varstr_col, varpsu_col]].copy()
        df_subset.columns = ['drug_name', 'spending', 'weight', 'varstr', 'varpsu']

        # Drop rows with missing weights or negative spending (MEPS uses negative codes for missing)
        df_subset = df_subset[(df_subset['weight'] > 0) & (df_subset['spending'] >= 0)]

        con.register('current_meps', df_subset)

        # 4. Filter and format using DuckDB
        q = f"""
        SELECT
            CASE {drug_case_sql} ELSE 'OTHER' END AS generic_name,
            '{year}' AS year,
            spending,
            weight,
            varstr,
            varpsu
        FROM current_meps
        """
        con.execute(f"CREATE TEMP TABLE meps_{year} AS {q}")
        con.unregister('current_meps')

    # 5. Union all years into a single row-level dataset
    tables = con.execute("SELECT table_name FROM information_schema.tables WHERE table_name LIKE 'meps_%'").fetchall()
    if not tables:
        raise ValueError("No MEPS data was successfully processed.")

    union_query = " UNION ALL ".join([f"SELECT * FROM {t[0]}" for t in tables])

    final_query = f"""
    SELECT generic_name, CAST(year AS INTEGER) as year, spending, weight, varstr, varpsu
    FROM ({union_query})
    ORDER BY generic_name, year
    """
    return final_query

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meps_folder', type=str, required=True)
    parser.add_argument('--meps_files', type=str, nargs='+', required=True)
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)

    args = parser.parse_args()
    print("Running MEPS Row-Level Extraction Pipeline...")

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    con = duckdb.connect(database=':memory:')
    query = process_meps_files(args.meps_folder, args.meps_files, args.targets, args.controls, drug_mappings, con)

    result_df = con.execute(query).df()

    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})
    result_df['role'] = result_df['generic_name'].map(role_map)

    result_df.to_csv(args.output, index=False)
    print(f"MEPS row-level survey data successfully saved to {args.output}")

if __name__ == "__main__":
    main()
