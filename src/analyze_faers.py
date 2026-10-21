import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml

def build_faers_queries(targets, controls, drug_mappings, files, min_year, max_year):
    """
    Builds the ordered SQL steps that turn raw FAERS quarters into case-level counts.

    Parameters:
    targets (list): Target drug names.
    controls (list): Control drug names.
    drug_mappings (dict): Dictionary mapping generic drugs to brand names.
    files (dict): Dictionary of file paths for DEMO, DRUG, REAC, and INDI.
    min_year (int): Start year for the analysis window.
    max_year (int): End year for the analysis window.
    """
    # Standardize CSV reading options for DuckDB
    csv_opts = "delim='$', quote='', header=True, all_varchar=True, normalize_names=True"
    demo_src = f"read_csv({files['DEMO']}, {csv_opts})"
    drug_src = f"read_csv({files['DRUG']}, {csv_opts})"
    reac_src = f"read_csv({files['REAC']}, {csv_opts})"
    indi_src = f"read_csv({files['INDI']}, {csv_opts})"

    # Regex to detect combination products (e.g., drug A / drug B)
    combo_ai = r"[\\/+&;]|\bAND\b"
    combo_name = r"[A-Za-z]{4,}\s*[/+&]\s*[A-Za-z]{4,}|\bAND\b"

    drug_when, product_when, drug_filter = [], [], []
    for d in targets + controls:
        terms = drug_mappings.get(d, [d])
        regex = f"(?i)\\b({'|'.join(re.escape(t) for t in terms)})\\b"

        # Positive match for the drug, negative match for combination delimiters
        pos = f"(REGEXP_MATCHES(dn, '{regex}') OR REGEXP_MATCHES(ai, '{regex}'))"
        neg = f"(REGEXP_MATCHES(dn, '(?i){combo_name}') OR REGEXP_MATCHES(ai, '(?i){combo_ai}'))"

        drug_when.append(f"WHEN ({pos} AND NOT {neg}) THEN '{d}'")
        drug_filter.append(f"({pos} AND NOT {neg})")

        # Brand terms (everything after the generic) label the specific product (e.g., Ozempic)
        for t in terms[1:]:
            product_when.append(f"WHEN REGEXP_MATCHES(dn, '(?i)\\b{re.escape(t)}\\b') THEN '{t}'")

    product_sql = f"CASE {' '.join(product_when)} ELSE 'GENERIC/UNSPECIFIED' END" if product_when else "'GENERIC/UNSPECIFIED'"

    # MedDRA preferred terms behind each outcome flag
    off_label = ['Off label use', 'Product use in unapproved indication', 'Prescription drug used without a prescription']
    med_core = ['Medication error', 'Suspected counterfeit product', 'Intentional product misuse', 'Wrong technique in product usage process']
    med_all = med_core + ['Underdose', 'Overdose']
    gi = ['Nausea', 'Vomiting', 'Diarrhoea']
    pt_list = lambda pts: ", ".join(f"'{p.lower()}'" for p in pts)

    ddl = {}

    # STEP A: Canonical DEMO Deduplication
    # Keep the latest version of each case. Maps 'OT' (Other Health Professional) to Healthcare Professional.
    ddl['canonical_demo'] = f"""
        CREATE OR REPLACE TEMP TABLE canonical_demo AS
        SELECT primaryid, occp_cod,
               CASE WHEN occp_cod IN ('MD', 'PH', 'OT') THEN 'Healthcare professional'
                    WHEN occp_cod = 'CN' THEN 'Consumer'
                    WHEN occp_cod = 'LW' THEN 'Lawyer'
                    ELSE 'Other/Unknown' END AS reporter_group,
               init_year
        FROM (
            SELECT primaryid, TRIM(occp_cod) AS occp_cod, TRY_CAST(SUBSTRING(init_fda_dt, 1, 4) AS INT) AS init_year,
                   ROW_NUMBER() OVER (PARTITION BY caseid ORDER BY TRY_CAST(caseversion AS INT) DESC, fda_dt DESC, primaryid DESC) AS rn
            FROM {demo_src}
        )
        WHERE rn = 1 AND init_year BETWEEN {min_year} AND {max_year}
    """

    # STEP B: Filtered DRUG Table
    # Primary suspect drugs only; COALESCE prevents rows with empty prod_ai from being silently dropped
    ddl['filtered_drug'] = f"""
        CREATE OR REPLACE TEMP TABLE filtered_drug AS
        SELECT primaryid, drug_seq, dn AS drugname, ai AS prod_ai,
               CASE {' '.join(drug_when)} END AS drug,
               {product_sql} AS product
        FROM (
            SELECT primaryid, drug_seq, role_cod, COALESCE(drugname, '') AS dn, COALESCE(prod_ai, '') AS ai
            FROM {drug_src}
        )
        WHERE TRIM(role_cod) = 'PS' AND ({' OR '.join(drug_filter)})
    """

    ddl['indi_raw'] = f"""
        CREATE OR REPLACE TEMP TABLE indi_raw AS
        SELECT primaryid, indi_drug_seq, indi_pt FROM {indi_src}
    """

    # STEP C: Indication flags per drug line
    ddl['resolved_indi'] = """
        CREATE OR REPLACE TEMP TABLE resolved_indi AS
        SELECT primaryid, indi_drug_seq,
               MAX(CASE WHEN REGEXP_MATCHES(indi_pt, '(?i)weight (loss|control|management|decreased)|obes|overweight') THEN 1 ELSE 0 END) AS is_weight_loss,
               MAX(CASE WHEN REGEXP_MATCHES(indi_pt, '(?i)diabet') THEN 1 ELSE 0 END) AS is_diabetes
        FROM indi_raw
        GROUP BY primaryid, indi_drug_seq
    """

    # STEP D: Outcome flags per case
    ddl['resolved_reac'] = f"""
        CREATE OR REPLACE TEMP TABLE resolved_reac AS
        SELECT primaryid,
               MAX(CASE WHEN LOWER(pt) IN ({pt_list(off_label)}) THEN 1 ELSE 0 END) AS is_off_label,
               MAX(CASE WHEN LOWER(pt) IN ({pt_list(med_all)}) THEN 1 ELSE 0 END) AS is_med_error,
               MAX(CASE WHEN LOWER(pt) IN ({pt_list(med_core)}) THEN 1 ELSE 0 END) AS is_med_core,
               MAX(CASE WHEN LOWER(pt) IN ({pt_list(gi)}) THEN 1 ELSE 0 END) AS is_gi
        FROM {reac_src}
        GROUP BY primaryid
    """

    # Audit queries printed for manual inspection
    audit = {}
    audit['Top matched drugname / prod_ai strings per cohort'] = """
        SELECT * FROM (
            SELECT drug, drugname, prod_ai, COUNT(*) AS n,
                   ROW_NUMBER() OVER (PARTITION BY drug ORDER BY COUNT(*) DESC) AS rk
            FROM filtered_drug GROUP BY drug, drugname, prod_ai
        ) WHERE rk <= 8 ORDER BY drug, rk
    """
    audit['Top indication terms for target drugs'] = f"""
        SELECT i.indi_pt, COUNT(*) AS n
        FROM indi_raw i JOIN filtered_drug fd ON i.primaryid = fd.primaryid AND i.indi_drug_seq = fd.drug_seq
        WHERE fd.drug IN ({', '.join(repr(t) for t in targets)})
        GROUP BY i.indi_pt ORDER BY n DESC LIMIT 25
    """
    audit['Reporter codes (OCCP_COD) and the group each maps to'] = """
        SELECT occp_cod, reporter_group, COUNT(*) AS n FROM canonical_demo GROUP BY 1, 2 ORDER BY n DESC
    """
    audit['Share of target drug lines with no weight-loss or diabetes indication, by year and reporter'] = f"""
        SELECT demo.init_year AS year, demo.reporter_group, COUNT(*) AS drug_lines,
               ROUND(AVG(CASE WHEN COALESCE(indi.is_weight_loss, 0) + COALESCE(indi.is_diabetes, 0) = 0 THEN 1.0 ELSE 0.0 END), 3) AS unclassified_share
        FROM canonical_demo demo
        JOIN filtered_drug fd ON demo.primaryid = fd.primaryid
        LEFT JOIN resolved_indi indi ON fd.primaryid = indi.primaryid AND fd.drug_seq = indi.indi_drug_seq
        WHERE fd.drug IN ({', '.join(repr(t) for t in targets)})
        GROUP BY 1, 2 ORDER BY 1, 2
    """
    audit['Outcome PTs present in REAC (case-insensitive)'] = f"""
        SELECT LOWER(pt) AS pt, COUNT(*) AS n FROM {reac_src}
        WHERE LOWER(pt) IN ({pt_list(off_label + med_all + gi)})
        GROUP BY 1 ORDER BY n DESC
    """

    # STEP E: Final Master Join and Aggregation
    # One row per year x reporter group x drug x product x indication; outcomes counted as distinct cases to prevent fan-out bug
    final = """
        SELECT demo.init_year AS year,
               demo.reporter_group,
               fd.drug,
               fd.product,
               CASE WHEN indi.is_weight_loss = 1 AND indi.is_diabetes = 1 THEN 'Dual Indication'
                    WHEN indi.is_weight_loss = 1 THEN 'Weight Loss'
                    WHEN indi.is_diabetes = 1 THEN 'Diabetes'
                    ELSE 'Other/Unknown' END AS indication,
               COUNT(DISTINCT demo.primaryid) AS total_cases,
               COUNT(DISTINCT CASE WHEN reac.is_off_label = 1 THEN demo.primaryid END) AS off_label_cases,
               COUNT(DISTINCT CASE WHEN reac.is_med_error = 1 THEN demo.primaryid END) AS med_error_cases,
               COUNT(DISTINCT CASE WHEN reac.is_med_core = 1 THEN demo.primaryid END) AS med_error_core_cases,
               COUNT(DISTINCT CASE WHEN reac.is_gi = 1 THEN demo.primaryid END) AS gi_cases
        FROM canonical_demo demo
        JOIN filtered_drug fd ON demo.primaryid = fd.primaryid
        LEFT JOIN resolved_indi indi ON fd.primaryid = indi.primaryid AND fd.drug_seq = indi.indi_drug_seq
        LEFT JOIN resolved_reac reac ON demo.primaryid = reac.primaryid
        GROUP BY 1, 2, 3, 4, 5
        ORDER BY 1, 3, 5, 2, 4
    """
    return ddl, audit, final

def main():
    parser = argparse.ArgumentParser(description="FAERS Pharmacovigilance Pipeline")
    parser.add_argument('--targets', type=str, nargs='+', required=True)
    parser.add_argument('--controls', type=str, nargs='*', default=[])
    parser.add_argument('--faers_folder', type=str, required=True)
    parser.add_argument('--faers_quarters', type=str, nargs='+', required=True)
    parser.add_argument('--config_file', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    parser.add_argument('--db_path', type=str, default='faers_processing.db')

    args = parser.parse_args()

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    drug_mappings = config.get('drug_mappings', {})

    # Exactly one file per table per quarter, otherwise fail loudly
    files = {}
    for prefix in ['DEMO', 'DRUG', 'REAC', 'INDI']:
        files[prefix] = []
        for q in args.faers_quarters:
            path = Path(args.faers_folder) / q
            matched = [p for p in path.glob('*') if p.is_file() and p.name.lower().startswith(prefix.lower()) and p.suffix.lower() == '.txt']
            if len(matched) != 1:
                raise FileNotFoundError(f"Expected exactly one {prefix} file in {path}, found {len(matched)}")
            files[prefix].append(matched[0].as_posix())

    # Analysis window comes from the quarter folder names
    years = [int(re.search(r'(\d{4})', q).group(1)) for q in args.faers_quarters]
    min_year, max_year = min(years), max(years)

    ddl, audit, final = build_faers_queries(args.targets, args.controls, drug_mappings, files, min_year, max_year)

    # File-backed DuckDB so 20 quarters do not have to fit in RAM
    con = duckdb.connect(args.db_path)
    try:
        for name, sql in ddl.items():
            print(f"Building {name}...")
            con.execute(sql)
        for title, sql in audit.items():
            print(f"\n--- {title} ---")
            print(con.execute(sql).df().to_string(index=False))
        df = con.execute(final).df()
    finally:
        con.close()
        Path(args.db_path).unlink(missing_ok=True)

    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})
    df['role'] = df['drug'].map(role_map)

    # Post-execution assertion: Flag counts can never exceed total cases (verifies fan-out bug is fixed)
    for c in ['off_label_cases', 'med_error_cases', 'med_error_core_cases', 'gi_cases']:
        if (df[c] > df['total_cases']).any():
            raise ValueError(f"{c} exceeds total_cases in some rows")

    df.to_csv(args.output, index=False)
    print(f"\nFAERS case counts saved to {args.output} ({len(df)} rows, years {min_year}-{max_year})")

if __name__ == "__main__":
    main()
