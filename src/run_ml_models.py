import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

def run_logistic_regression(faers_df):
    # Extract year as a continuous variable
    faers_df['year_num'] = faers_df['quarter_year'].str[:2].astype(int)
    # Center the year to 2020 (0) to fix the intercept issue Claude found
    faers_df['year_centered'] = faers_df['year_num'] - 20

    faers_df['indication'] = 'Unknown'
    faers_df.loc[faers_df['drug_and_indication'].str.contains('Weight Loss'), 'indication'] = 'Weight Loss'
    faers_df.loc[faers_df['drug_and_indication'].str.contains('Diabetes'), 'indication'] = 'Diabetes'

    model_df = faers_df[faers_df['indication'] != 'Unknown'].copy()

    # CLAUDE FIX: GLM requires a 2D array of [successes, failures] for binomial data
    model_df['failures'] = model_df['total_cases'] - model_df['misuse_cases']

    # We restrict to Semaglutide to test the interaction of the indication directly
    sema_df = model_df[model_df['drug'] == 'SEMAGLUTIDE'].copy()

    if len(sema_df) == 0:
        return None

    model = smf.glm(
        formula="misuse_cases + failures ~ C(indication, Treatment('Diabetes')) + year_centered",
        data=sema_df,
        family=sm.families.Binomial()
    ).fit()

    return model

def calculate_reporting_odds_ratio(faers_df):
    wl_df = faers_df[faers_df['drug_and_indication'] == 'SEMAGLUTIDE (Weight Loss)']
    a = wl_df['misuse_cases'].sum()
    b = wl_df['total_cases'].sum() - a

    dm_df = faers_df[faers_df['drug_and_indication'] == 'SEMAGLUTIDE (Diabetes)']
    c = dm_df['misuse_cases'].sum()
    d = dm_df['total_cases'].sum() - c

    # CLAUDE FIX: Haldane 0.5 correction for zero cells
    if a == 0 or b == 0 or c == 0 or d == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    ror = (a * d) / (b * c)

    se_ln_ror = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lower = np.exp(np.log(ror) - 1.96 * se_ln_ror)
    ci_upper = np.exp(np.log(ror) + 1.96 * se_ln_ror)

    return ror, ci_lower, ci_upper

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cms_file', required=True)
    parser.add_argument('--faers_file', required=True)
    parser.add_argument('--output_report', required=True)
    args = parser.parse_args()

    faers_df = pd.read_csv(args.faers_file)

    log_model = run_logistic_regression(faers_df)
    ror, ci_low, ci_high = calculate_reporting_odds_ratio(faers_df)

    report = []
    report.append("=========================================================")
    report.append("      CAPSTONE ML & STATISTICAL INSIGHTS REPORT          ")
    report.append("=========================================================\n")

    report.append("1. MACHINE LEARNING: LOGISTIC REGRESSION (RISK CLASSIFICATION)")
    report.append("---------------------------------------------------------")
    report.append("Model: Misuse Cases ~ Indication + Year (Restricted to Semaglutide)")
    if log_model:
        report.append("\nModel Coefficients (Log-Odds):")
        for index, row in log_model.summary2().tables[1].iterrows():
            coef = row['Coef.']
            pval = row['P>|z|']
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
            report.append(f"  {index}: {coef:.4f} (p={pval:.4e}) {sig}")
    else:
        report.append("Insufficient data to fit model.")

    report.append("\n2. PHARMACOVIGILANCE: REPORTING ODDS RATIO (ROR)")
    report.append("---------------------------------------------------------")
    report.append("Comparison: Semaglutide (Weight Loss) vs. Semaglutide (Diabetes)")
    report.append(f"Reporting Odds Ratio (ROR): {ror:.2f}")
    report.append(f"95% Confidence Interval: [{ci_low:.2f}, {ci_high:.2f}]")

    with open(args.output_report, 'w') as f:
        f.write("\n".join(report))

if __name__ == "__main__":
    main()
