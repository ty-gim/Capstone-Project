import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml

def build_and_run_meps_pipeline(meps_folder, meps_files, targets, controls, drug_mappings, con):
    file_paths = [str(Path(meps_folder) / f) for f in meps_files]

    for fp in file_paths:
        print(f"Loading {fp} into memory...")

        try:
            df = pd.read_stata(fp)
            df.columns = [c.lower() for c in df.columns]
        except Exception as e:
            print(f"Error reading {fp}: {e}")
            continue

        # CLAUDE FIX: Prioritize rxdrgnam over rxname to avoid combination-drug pollution
        if 'rxdrgnam' in df.columns:
            drug_col = 'rxdrgnam'
        elif 'rxname' in df.columns:
            drug_col = 'rxname'
        else:
            drug_col = None

        total_exp_col = next((c for c in df.columns if c.startswith('rxxp') and c.endswith('x')), None)
        weight_col = next((c for c in df.columns if c.startswith('perwt') and c.endswith('f')), None)

        if not drug_col or not total_exp_col:
            print(f"Warning: Missing required MEPS financial columns in {fp}.")
            continue

        if not weight_col:
            print(f"Warning: Missing weight column in {fp}. Defaulting to 1.0 (Unweighted).")
            df['default_weight'] = 1.0
            weight_col = 'default_weight'

        year_match = re.search(r'rxxp([0-9]{2})x', total_exp_col)
        year = f"20{year_match.group(1)}" if year_match else "Unknown"

        # 2. Sanitize Data
        df_subset = df[[drug_col, total_exp_col, weight_col]].copy()
        df_subset[drug_col] = df_subset[drug_col].astype(str)
        df_subset[total_exp_col] = pd.to_numeric(df_subset[total_exp_col], errors='coerce').fillna(0)
        df_subset[weight_col] = pd.to_numeric(df_subset[weight_col], errors='coerce').fillna(0)

        # 3. Build dynamic CASE WHEN using the config mappings for MULTIPLE targets and controls
        case_when = ""
        all_where_conds = []

        # Handle Targets
        for d in targets:
            terms = drug_mappings.get(d, [d])
            conds = " OR ".join([f"{drug_col} ILIKE '%{t}%'" for t in terms])
            case_when += f"        WHEN {conds} THEN '{d}'\n"
            all_where_conds.extend([f"{drug_col} ILIKE '%{t}%'" for t in terms])

        # Handle Controls
        for d in controls:
            terms = drug_mappings.get(d, [d])
            conds = " OR ".join([f"{drug_col} ILIKE '%{t}%'" for t in terms])
            case_when += f"        WHEN {conds} THEN '{d}'\n"
            all_where_conds.extend([f"{drug_col} ILIKE '%{t}%'" for t in terms])

        where_conditions = " OR ".join(all_where_conds)

        # 4. Query the subset
        q = f"""
        SELECT
            CASE
                {case_when}
                ELSE 'OTHER'
            END AS generic_name,
            '{year}' AS year,
            {total_exp_col} AS total_spending,
            {weight_col} AS person_weight
        FROM df_subset
        WHERE ({where_conditions})
        """

        con.execute(f"CREATE TEMP TABLE meps_{year} AS {q}")

    tables = con.execute("SELECT table_name FROM information_schema.tables WHERE table_name LIKE 'meps_%'").fetchall()
    if not tables:
        raise ValueError("No MEPS data was successfully processed.")

    union_query = " UNION ALL ".join([f"SELECT * FROM {t[0]}" for t in tables])

    # 5. Calculate True National Estimates using the Survey Weights
    final_query = f"""
    WITH CombinedData AS (
        {union_query}
    )
    SELECT
        generic_name,
        year,
        -- National Prescriptions = Sum of the weights
        SUM(person_weight) AS national_prescription_count,
        -- National Spending = Spending * Weight
        SUM(total_spending * person_weight) AS national_system_spending,
        -- Average Cost = National Spending / National Prescriptions
        CASE
            WHEN SUM(person_weight) = 0 THEN 0
            ELSE SUM(total_spending * person_weight) / SUM(person_weight)
        END AS avg_cost_per_prescription
    FROM CombinedData
    GROUP BY generic_name, year
    ORDER BY generic_name, year;
    """
    return final_query

def main():
    parser = argparse.ArgumentParser(description="MEPS Financial Pipeline")
    parser.add_argument('--meps_folder', type=str, required=True)
    parser.add_argument('--meps_files', type=str, nargs='+', required=True)
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)

    args = parser.parse_args()
    print(f"Running MEPS Financial Analysis for {args.targets} vs {args.controls}...")

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    # Map roles for R visualization
    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})

    con = duckdb.connect(database=':memory:')
    query = build_and_run_meps_pipeline(args.meps_folder, args.meps_files, args.targets, args.controls, drug_mappings, con)

    result_df = con.execute(query).df()
    result_df['role'] = result_df['generic_name'].map(role_map)
    result_df.to_csv(args.output, index=False)
    print(f"MEPS Financial data successfully saved to {args.output}")

if __name__ == "__main__":
    main()
