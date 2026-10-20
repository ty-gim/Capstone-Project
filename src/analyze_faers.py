import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml
import sys

def get_files_by_prefix(base_folder, quarters, prefix):
    """
    Strictly discovers files. If a quarter is missing a required file,
    it fails loudly rather than silently skipping.
    """
    files = []
    for q in quarters:
        path = Path(base_folder) / q
        if not path.exists():
            raise FileNotFoundError(f"Expected directory does not exist: {path}")

        # FAERS ASCII files are typically named like DEMO20Q1.txt
        matched_files = [p for p in path.glob('*') if p.is_file() and p.name.lower().startswith(prefix.lower()) and p.suffix.lower() == '.txt']

        if not matched_files:
            raise FileNotFoundError(f"Missing {prefix} file in {path}")
        if len(matched_files) > 1:
            raise ValueError(f"Multiple {prefix} files found in {path}. Expected exactly one.")

        files.append(str(matched_files[0].as_posix()))
    return files

def build_regex_conditions(drugs, drug_mappings):
    """
    Builds strict RE2-compatible regex conditions for the target drugs.
    """
    conditions = []
    for d in drugs:
        terms = drug_mappings.get(d, [d])
        # Escape terms and join with OR
        terms_regex = "|".join([re.escape(t) for t in terms])
        # \b ensures word boundaries. (?i) makes it case-insensitive.
        regex = f"(?i)\\b({terms_regex})\\b"
        conditions.append((d, regex))
    return conditions

def run_faers_pipeline(targets, controls, drug_mappings, faers_folder, faers_quarters, output_ml, output_visuals):
    # 1. Strict File Discovery
    print("Discovering and validating FAERS files...")
    demo_files = get_files_by_prefix(faers_folder, faers_quarters, 'DEMO')
    drug_files = get_files_by_prefix(faers_folder, faers_quarters, 'DRUG')
    reac_files = get_files_by_prefix(faers_folder, faers_quarters, 'REAC')
    indi_files = get_files_by_prefix(faers_folder, faers_quarters, 'INDI')

    # 2. Build Boolean SQL Logic for Drugs (RE2 Compatible)
    all_drugs = targets + controls
    regex_mapping = build_regex_conditions(all_drugs, drug_mappings)

    # Regex to detect combination drugs (e.g., SITAGLIPTIN/METFORMIN)
    combo_exclusion = r"[/+&]|\bAND\b"

    case_when_clauses = []
    drug_filter_conds = []

    for drug_name, regex in regex_mapping:
        # Positive match for the drug
        pos_match = f"(REGEXP_MATCHES(drugname, '{regex}') OR REGEXP_MATCHES(prod_ai, '{regex}'))"
        # Negative match for combination delimiters
        neg_match = f"(REGEXP_MATCHES(drugname, '(?i){combo_exclusion}') OR REGEXP_MATCHES(prod_ai, '(?i){combo_exclusion}'))"

        # Combine them: Must contain the drug, MUST NOT contain combination delimiters
        full_cond = f"({pos_match} AND NOT {neg_match})"

        case_when_clauses.append(f"WHEN {full_cond} THEN '{drug_name}'")
        drug_filter_conds.append(full_cond)

    drug_case_sql = "\n".join(case_when_clauses)
    drug_filter_sql = " OR ".join(drug_filter_conds)

    # 3. Initialize DuckDB (File-backed to prevent RAM blowouts on 20+ quarters)
    print("Initializing DuckDB out-of-core processing...")
    con = duckdb.connect('faers_processing.db')

    try:
        # STEP A: Canonical DEMO Deduplication
        print("Building Canonical DEMO table (Deduplicating case versions)...")
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE canonical_demo AS
            SELECT primaryid, caseid, occp_cod,
                   SUBSTRING(init_fda_dt, 1, 4) AS init_year
            FROM (
                SELECT primaryid, caseid, occp_cod, init_fda_dt, caseversion, fda_dt,
                       ROW_NUMBER() OVER (
                           PARTITION BY caseid
                           ORDER BY CAST(caseversion AS INT) DESC, fda_dt DESC
                       ) as rn
                FROM read_csv({demo_files}, delim='$', quote='', header=True, all_varchar=True, normalize_names=True)
            )
            WHERE rn = 1
        """)

        # STEP B: Filtered DRUG Table
        print("Filtering DRUG table...")
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE filtered_drug AS
            SELECT primaryid, drug_seq,
                   CASE {drug_case_sql} END AS drug_cohort
            FROM read_csv({drug_files}, delim='$', quote='', header=True, all_varchar=True, normalize_names=True)
            WHERE TRIM(role_cod) = 'PS'
              AND ({drug_filter_sql})
        """)

        # STEP C: Indication Resolution
        print("Resolving INDI table...")
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE resolved_indi AS
            SELECT primaryid, indi_drug_seq,
                   MAX(CASE WHEN REGEXP_MATCHES(indi_pt, '(?i)weight loss|obesity|weight management') THEN 1 ELSE 0 END) AS is_weight_loss,
                   MAX(CASE WHEN REGEXP_MATCHES(indi_pt, '(?i)diabet') THEN 1 ELSE 0 END) AS is_diabetes
            FROM read_csv({indi_files}, delim='$', quote='', header=True, all_varchar=True, normalize_names=True)
            GROUP BY primaryid, indi_drug_seq
        """)

        # STEP D: Reaction/Outcome Resolution
        print("Resolving REAC table...")
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE resolved_reac AS
            SELECT primaryid,
                   MAX(CASE WHEN pt IN ('Off label use', 'Product use in unapproved indication', 'Prescription drug used without a prescription') THEN 1 ELSE 0 END) AS is_off_label,
                   MAX(CASE WHEN pt IN ('Medication error', 'Suspected counterfeit product', 'Intentional product misuse', 'Wrong technique in product usage process', 'Underdose', 'Overdose') THEN 1 ELSE 0 END) AS is_med_error,
                   MAX(CASE WHEN pt IN ('Nausea', 'Vomiting', 'Diarrhoea') THEN 1 ELSE 0 END) AS is_gi_event
            FROM read_csv({reac_files}, delim='$', quote='', header=True, all_varchar=True, normalize_names=True)
            GROUP BY primaryid
        """)

        # STEP E: Final Master Join and Aggregation
        print("Executing Master Join and Aggregation...")
        final_query = """
            SELECT
                demo.init_year AS year,
                demo.occp_cod AS reporter_type,
                drug.drug_cohort AS drug,
                CASE
                    WHEN indi.is_weight_loss = 1 AND indi.is_diabetes = 1 THEN 'Dual Indication'
                    WHEN indi.is_weight_loss = 1 THEN 'Weight Loss'
                    WHEN indi.is_diabetes = 1 THEN 'Diabetes'
                    ELSE 'Other/Unknown'
                END AS indication,
                COUNT(DISTINCT demo.primaryid) AS total_cases,
                SUM(reac.is_off_label) AS off_label_cases,
                SUM(reac.is_med_error) AS med_error_cases,
                SUM(reac.is_gi_event) AS gi_cases
            FROM canonical_demo demo
            JOIN filtered_drug drug ON demo.primaryid = drug.primaryid
            LEFT JOIN resolved_indi indi ON drug.primaryid = indi.primaryid AND drug.drug_seq = indi.indi_drug_seq
            LEFT JOIN resolved_reac reac ON demo.primaryid = reac.primaryid
            WHERE demo.init_year IS NOT NULL
            GROUP BY 1, 2, 3, 4
            ORDER BY 1, 3, 4
        """

        df = con.execute(final_query).df()

    finally:
        con.close()
        Path('faers_processing.db').unlink(missing_ok=True)

    # 4. Map Roles and Save
    role_map = {d: 'target' for d in targets}
    role_map.update({d: 'control' for d in controls})
    df['role'] = df['drug'].map(role_map)

    print(f"Saving outputs to {output_ml} and {output_visuals}...")
    df.to_csv(output_ml, index=False)
    df.to_csv(output_visuals, index=False)
    print("FAERS Pipeline completed successfully.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--faers_folder', type=str, required=True)
    parser.add_argument('--faers_quarters', type=str, nargs='+', required=True)
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output_visuals', type=str, required=True)
    parser.add_argument('--output_ml', type=str, required=True)

    args = parser.parse_args()

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    run_faers_pipeline(
        targets=args.targets,
        controls=args.controls,
        drug_mappings=drug_mappings,
        faers_folder=args.faers_folder,
        faers_quarters=args.faers_quarters,
        output_ml=args.output_ml,
        output_visuals=args.output_visuals
    )

if __name__ == "__main__":
    main()
