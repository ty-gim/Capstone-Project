import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

def run_logistic_regression(faers_df):
    """
    Fits a Logistic Regression model to predict the probability of an
    'Inappropriate Use / Access' adverse event based on Drug, Indication, and Time.
    """
    # Define the signal based on MedDRA Preferred Terms related to access/administration
    inappropriate_use_terms = [
        'Off label use',
        'Product use in unapproved indication',
        'Prescription drug used without a prescription',
        'Suspected counterfeit product',
        'Medication error',
        'Intentional product misuse'
    ]

    # Create the binary target variable
    faers_df['is_inappropriate_use'] = faers_df['adverse_event'].isin(inappropriate_use_terms).astype(int)

    # Extract year as a continuous variable for the model
    faers_df['year_num'] = faers_df['quarter_year'].str[:2].astype(int)

    # Extract indication cleanly
    faers_df['indication'] = 'Unknown'
    faers_df.loc[faers_df['drug_and_indication'].str.contains('Weight Loss'), 'indication'] = 'Weight Loss'
    faers_df.loc[faers_df['drug_and_indication'].str.contains('Diabetes'), 'indication'] = 'Diabetes'

    # Filter out Unknowns for a clean model
    model_df = faers_df[faers_df['indication'] != 'Unknown'].copy()

    # Fit the Logistic Regression (Binomial GLM)
    # Formula: is_inappropriate_use ~ drug + indication + year_num
    model = smf.glm(
        formula="is_inappropriate_use ~ C(drug, Treatment('METFORMIN')) + C(indication, Treatment('Diabetes')) + year_num",
        data=model_df,
        family=sm.families.Binomial(),
        freq_weights=model_df['cases']
    ).fit()

    return model

def calculate_reporting_odds_ratio(faers_df):
    """
    Calculates the Reporting Odds Ratio (ROR) with 95% Confidence Intervals
    for Inappropriate Use events: Semaglutide (Weight Loss) vs Semaglutide (Diabetes).
    """
    inappropriate_use_terms = [
        'Off label use', 'Product use in unapproved indication',
        'Prescription drug used without a prescription', 'Suspected counterfeit product'
    ]

    # Cohort 1: Semaglutide Weight Loss
    wl_df = faers_df[faers_df['drug_and_indication'] == 'SEMAGLUTIDE (Weight Loss)']
    a = wl_df[wl_df['adverse_event'].isin(inappropriate_use_terms)]['cases'].sum()
    b = wl_df[~wl_df['adverse_event'].isin(inappropriate_use_terms)]['cases'].sum()

    # Cohort 2: Semaglutide Diabetes
    dm_df = faers_df[faers_df['drug_and_indication'] == 'SEMAGLUTIDE (Diabetes)']
    c = dm_df[dm_df['adverse_event'].isin(inappropriate_use_terms)]['cases'].sum()
    d = dm_df[~dm_df['adverse_event'].isin(inappropriate_use_terms)]['cases'].sum()

    if b == 0 or c == 0 or d == 0:
        return None, None, None

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
    report.append("Model: is_inappropriate_use ~ Drug + Indication + Year")
    report.append("Objective: Predict the probability that an adverse event is related to inappropriate use or administration.")
    report.append("\nModel Coefficients (Log-Odds):")

    for index, row in log_model.summary2().tables[1].iterrows():
        coef = row['Coef.']
        pval = row['P>|z|']
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        report.append(f"  {index}: {coef:.4f} (p={pval:.4e}) {sig}")

    report.append("\nInterpretation:")
    report.append("Positive coefficients indicate a higher probability of inappropriate use events.")
    report.append("If 'Indication[T.Weight Loss]' is positive and significant, it suggests the weight-loss")
    report.append("indication is an independent predictor of these events, controlling for the drug itself.\n")

    report.append("2. PHARMACOVIGILANCE: REPORTING ODDS RATIO (ROR)")
    report.append("---------------------------------------------------------")
    report.append("Comparison: Semaglutide (Weight Loss) vs. Semaglutide (Diabetes)")
    report.append("Outcome: Inappropriate Use and Administration Events")
    if ror:
        report.append(f"Reporting Odds Ratio (ROR): {ror:.2f}")
        report.append(f"95% Confidence Interval: [{ci_low:.2f}, {ci_high:.2f}]")
        if ci_low > 1.0:
            report.append("Conclusion: STATISTICALLY SIGNIFICANT SIGNAL DETECTED.")
            report.append("Patients taking Semaglutide for weight loss have significantly higher odds")
            report.append("of reporting inappropriate use events than patients taking it for diabetes.")
    else:
        report.append("Insufficient data to calculate ROR.")

    with open(args.output_report, 'w') as f:
        f.write("\n".join(report))

if __name__ == "__main__":
    main()
