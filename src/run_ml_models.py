import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

def run_logistic_regression(faers_df):
    # 1. Filter for Semaglutide and valid indications
    df = faers_df[(faers_df['drug'] == 'SEMAGLUTIDE') &
                  (faers_df['indication'].isin(['Weight Loss', 'Diabetes']))].copy()

    if len(df) == 0:
        return None

    # 2. Center year to 2020 for a meaningful intercept
    df['year'] = pd.to_numeric(df['year'], errors='coerce')
    df = df.dropna(subset=['year'])
    df['year_centered'] = df['year'] - 2020

    # 3. Define the outcome: Medication Errors (Avoids definitional circularity of 'Off-label')
    df['failures'] = df['total_cases'] - df['med_error_cases']

    # 4. Fit Grouped Binomial GLM
    # We include C(reporter_type) to adjust for stimulated reporting (media/lawyers)
    # We use cov_type='HC0' (Robust Standard Errors) to account for overdispersion/clustering
    model = smf.glm(
        formula="med_error_cases + failures ~ C(indication, Treatment('Diabetes')) + year_centered + C(reporter_type)",
        data=df,
        family=sm.families.Binomial()
    ).fit(cov_type='HC0')

    return model

def calculate_reporting_odds_ratio(faers_df):
    # Filter for Semaglutide and valid indications
    df = faers_df[(faers_df['drug'] == 'SEMAGLUTIDE') &
                  (faers_df['indication'].isin(['Weight Loss', 'Diabetes']))].copy()

    # Aggregate across years and reporter types for the crude ROR
    agg = df.groupby('indication')[['med_error_cases', 'total_cases']].sum()

    if 'Weight Loss' not in agg.index or 'Diabetes' not in agg.index:
        return np.nan, np.nan, np.nan

    a = agg.loc['Weight Loss', 'med_error_cases']
    b = agg.loc['Weight Loss', 'total_cases'] - a

    c = agg.loc['Diabetes', 'med_error_cases']
    d = agg.loc['Diabetes', 'total_cases'] - c

    # Haldane-Anscombe 0.5 correction for zero cells
    if a == 0 or b == 0 or c == 0 or d == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    ror = (a * d) / (b * c)
    se_ln_ror = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lower = np.exp(np.log(ror) - 1.96 * se_ln_ror)
    ci_upper = np.exp(np.log(ror) + 1.96 * se_ln_ror)

    return ror, ci_lower, ci_upper

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cms_file', required=True) # Kept to satisfy current Snakefile
    parser.add_argument('--faers_file', required=True)
    parser.add_argument('--output_report', required=True)
    args = parser.parse_args()

    faers_df = pd.read_csv(args.faers_file)

    log_model = run_logistic_regression(faers_df)
    ror, ci_low, ci_high = calculate_reporting_odds_ratio(faers_df)

    report = []
    report.append("=========================================================")
    report.append("      CAPSTONE STATISTICAL INSIGHTS REPORT               ")
    report.append("=========================================================\n")

    report.append("1. INFERENTIAL GROUPED-BINOMIAL GLM")
    report.append("---------------------------------------------------------")
    report.append("Model: med_error_cases ~ Indication + Year + Reporter Type")
    report.append("Note: Adjusted for stimulated reporting via reporter_type.")
    report.append("Note: Standard errors are robust (HC0) to account for overdispersion.\n")

    if log_model:
        report.append("Model Coefficients (Log-Odds):")
        # Extract robust standard errors and p-values
        for index, row in log_model.summary2().tables[1].iterrows():
            coef = row['Coef.']
            pval = row['P>|z|']
            # Convert log-odds to Odds Ratio for easier interpretation
            odds_ratio = np.exp(coef)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
            report.append(f"  {index}: Coef={coef:.4f} | OR={odds_ratio:.2f} | (p={pval:.4e}) {sig}")
    else:
        report.append("Insufficient data to fit model.")

    report.append("\n2. PHARMACOVIGILANCE: DISPROPORTIONALITY ANALYSIS")
    report.append("---------------------------------------------------------")
    report.append("Metric: Within-drug indication-contrast reporting odds ratio (ROR)")
    report.append("Comparison: Semaglutide (Weight Loss) vs. Semaglutide (Diabetes)")
    report.append("Outcome: Medication Errors & Counterfeits (Excludes 'Off-label' to prevent circularity)\n")

    if not np.isnan(ror):
        report.append(f"Reporting Odds Ratio (ROR): {ror:.2f}")
        report.append(f"95% Confidence Interval: [{ci_low:.2f}, {ci_high:.2f}]")
    else:
        report.append("Insufficient data to calculate ROR.")

    with open(args.output_report, 'w') as f:
        f.write("\n".join(report))

    print(f"Statistical report successfully generated at {args.output_report}")

if __name__ == "__main__":
    main()
