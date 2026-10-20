import argparse
import duckdb
import pandas as pd
from pathlib import Path
import re
import yaml

def get_single_file(base_folder, quarter, prefix):
    path = Path(base_folder) / quarter
    for p in path.glob('*'):
        if p.is_file() and p.name.lower().startswith(prefix.lower()) and p.suffix.lower() == '.txt':
            return str(p.as_posix())
    return None

def build_faers_query(targets, controls, drug_mappings, faers_folder, faers_quarters):
    quarter_queries = []

    for q in faers_quarters:
        drug_file = get_single_file(faers_folder, q, 'DRUG')
        reac_file = get_single_file(faers_folder, q, 'REAC')
        indi_file = get_single_file(faers_folder, q, 'INDI')

        if not drug_file or not reac_file:
            continue

        quarter_match = re.search(r'([0-9]{2}Q[1-4])', drug_file.upper())
        quarter_str = quarter_match.group(1) if quarter_match else q

        headers = {}
        try:
            for key, fpath in [('drug', drug_file), ('reac', reac_file), ('indi', indi_file)]:
                if fpath:
                    with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                        headers[key] = f.readline().strip().replace('\ufeff', '').replace('"', '').lower().split('$')
        except Exception:
            continue

        id_cols = {k: next((c for c in h if c in ['primaryid', 'isr']), None) for k, h in headers.items()}
        search_cols = [c for c in headers['drug'] if c in ['prod_ai', 'drugname']]

        if not id_cols['drug'] or not id_cols['reac'] or not search_cols or 'role_cod' not in headers['drug'] or 'pt' not in headers['reac']:
            continue

        case_when = ""
        where_conditions = []

        for d in targets + controls:
            terms = drug_mappings.get(d, [d])
            drug_conds = " OR ".join([f"d.{c} ILIKE '%{t}%'" for c in search_cols for t in terms])
            case_when += f"\n            WHEN {drug_conds} THEN '{d}'"
            where_conditions.append(f"({drug_conds})")

        all_conds = " OR ".join(where_conditions)

        # CLAUDE FIX: Explicitly define Diabetes vs Weight Loss vs Unknown
        if indi_file and id_cols['indi'] and 'indi_pt' in headers['indi'] and 'indi_drug_seq' in headers['indi'] and 'drug_seq' in headers['drug']:
            indi_join = f"""
            LEFT JOIN (
                SELECT
                    {id_cols['indi']} AS join_id,
                    indi_drug_seq,
                    MAX(CASE WHEN indi_pt ILIKE '%weight%' OR indi_pt ILIKE '%obesity%' THEN 1 ELSE 0 END) AS is_weight_loss,
                    MAX(CASE WHEN indi_pt ILIKE '%diabet%' THEN 1 ELSE 0 END) AS is_diabetes
                FROM read_csv('{indi_file}', delim='$', quote='', header=True, ignore_errors=True, normalize_names=True, all_varchar=True)
                GROUP BY {id_cols['indi']}, indi_drug_seq
            ) i ON d.{id_cols['drug']} = i.join_id AND d.drug_seq = i.indi_drug_seq
            """
            indication_logic = """
                CASE
                    WHEN i.is_weight_loss = 1 THEN ' (Weight Loss)'
                    WHEN i.is_diabetes = 1 THEN ' (Diabetes)'
                    ELSE ' (Unknown/Other)'
                END
            """
        else:
            indi_join = ""
            indication_logic = "' (Unknown/Other)'"

        q_sql = f"""
        SELECT
            '{quarter_str}' AS quarter_year,
            CASE
                {case_when}
            END || {indication_logic} AS drug_and_indication,
            r.pt AS adverse_event,
            COUNT(DISTINCT d.{id_cols['drug']}) AS cases
        FROM read_csv('{drug_file}', delim='$', quote='', header=True, ignore_errors=True, normalize_names=True, all_varchar=True) AS d
        JOIN read_csv('{reac_file}', delim='$', quote='', header=True, ignore_errors=True, normalize_names=True, all_varchar=True) AS r
          ON d.{id_cols['drug']} = r.{id_cols['reac']}
        {indi_join}
        WHERE TRIM(d.role_cod) = 'PS'
          AND ({all_conds})
        GROUP BY 1, 2, 3
        """
        quarter_queries.append(q_sql)

    union_query = " UNION ALL ".join(quarter_queries)

    # Output the FULL dataset for ML
    full_query = f"""
    WITH AECounts AS ({union_query})
    SELECT quarter_year, drug_and_indication, adverse_event, SUM(cases) as cases
    FROM AECounts
    GROUP BY quarter_year, drug_and_indication, adverse_event
    """
    return full_query

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

    role_map = {d: 'target' for d in args.targets}
    role_map.update({d: 'control' for d in args.controls})

    con = duckdb.connect(database=':memory:')
    query = build_faers_query(args.targets, args.controls, drug_mappings, args.faers_folder, args.faers_quarters)

    # Get full data
    full_df = con.execute(query).df()
    full_df['drug'] = full_df['drug_and_indication'].str.rsplit(' (', n=1).str[0]
    full_df['role'] = full_df['drug'].map(role_map)

    # Save full data for ML
    full_df.to_csv(args.output_ml, index=False)

    # CLAUDE FIX: Deterministic Tie-Breaking for the Top 5 Visuals
    # Sort by cases descending, then adverse_event alphabetically ascending
    full_df = full_df.sort_values(['quarter_year', 'drug_and_indication', 'cases', 'adverse_event'], ascending=[True, True, False, True])
    full_df['ae_rank'] = full_df.groupby(['quarter_year', 'drug_and_indication']).cumcount() + 1

    top5_df = full_df[full_df['ae_rank'] <= 5]
    top5_df.to_csv(args.output_visuals, index=False)

if __name__ == "__main__":
    main()
