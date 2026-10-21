import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm, chi2

def main():
    """
    Executes statistical models to contrast reporting patterns between Weight Loss and Diabetes indications.
    Calculates Crude ROR, Quasi-Binomial GLM, Stratified GLMs, and an Interaction Model.
    """
    parser = argparse.ArgumentParser(description="Weight-loss vs diabetes reporting contrasts in FAERS")
    parser.add_argument('--faers_file', required=True)
    parser.add_argument('--targets', nargs='+', required=True)
    parser.add_argument('--comparators', nargs='*', default=[])
    parser.add_argument('--brand_scopes', nargs='*', default=[])
    parser.add_argument('--exclude_products', nargs='*', default=[])
    parser.add_argument('--output_csv', required=True)
    parser.add_argument('--output_report', required=True)
    args = parser.parse_args()

    raw = pd.read_csv(args.faers_file)
    raw['reporter_group'] = raw['reporter_group'].fillna('Other/Unknown')

    # Outcome columns and their labels
    outcomes = {
        'med_error_cases': 'Medication errors & counterfeits (primary)',
        'med_error_core_cases': 'Med errors excluding over/underdose (sensitivity)',
        'gi_cases': 'GI events (calibration)',
        'off_label_cases': 'Off-label use (manipulation check)',
    }
    cols = ['total_cases'] + list(outcomes)
    out_cols = ['target', 'scope', 'outcome', 'method', 'odds_ratio', 'ci_low', 'ci_high', 'p_value', 'events_wl', 'n_wl', 'events_dm', 'n_dm', 'note']

    # Only the two clean indication strata enter the contrast; Other/Unknown and Dual are audited in the report
    df = raw[raw['indication'].isin(['Weight Loss', 'Diabetes'])].copy()
    rows, report = [], []

    report.append("FAERS weight-loss vs diabetes reporting contrast")
    report.append("=" * 60)
    report.append("OR > 1 means weight-loss reports have higher odds of including the event than diabetes reports.")
    report.append("These are reporting odds, not risks; FAERS cannot estimate incidence.\n")

    for target in args.targets:
        # Data audit: indication mix by year for the target
        base_all = raw[raw['drug'] == target]
        mix = base_all.pivot_table(index='year', columns='indication', values='total_cases', aggfunc='sum', fill_value=0)
        mix['Unclassified share'] = (mix.get('Other/Unknown', 0) + mix.get('Dual Indication', 0)) / mix.sum(axis=1)
        report.append(f"{target}: cases by year and indication")
        report.append(mix.to_string(float_format=lambda x: f"{x:.2f}") + "\n")

        base = df[df['drug'] == target]
        scopes = {'All products': base}

        # Brand-only scopes isolate one product (e.g., Ozempic weight-loss vs Ozempic diabetes) to truly isolate off-label use
        for brand in args.brand_scopes:
            if (base['product'] == brand).any():
                scopes[f"{brand} only"] = base[base['product'] == brand]

        # Oral products (like Rybelsus) are dropped as a formulation sensitivity check
        if args.exclude_products:
            scopes[f"Excluding {'/'.join(args.exclude_products)}"] = base[~base['product'].isin(args.exclude_products)]

        for scope, sdf in scopes.items():
            tot_by = sdf.groupby('indication')['total_cases'].sum()
            if tot_by.get('Weight Loss', 0) < 30 or tot_by.get('Diabetes', 0) < 30:
                report.append(f"Scope skipped for {target} / {scope}: fewer than 30 cases in a stratum.")
                continue

            # Aggregate to the covariate level of the model
            cells = sdf.groupby(['year', 'reporter_group', 'indication'], as_index=False)[cols].sum()
            cells['wl'] = (cells['indication'] == 'Weight Loss').astype(int)
            cells['year_f'] = cells['year'].astype(str)

            # Reporter groups with at least 20 cases in both strata get their own stratified estimate
            gsz = cells.pivot_table(index='reporter_group', columns='wl', values='total_cases', aggfunc='sum', fill_value=0)
            groups = [g for g in gsz.index if gsz.loc[g].get(0, 0) >= 20 and gsz.loc[g].get(1, 0) >= 20]

            for col, label in outcomes.items():
                cells['fail'] = cells['total_cases'] - cells[col]
                strata = [('ALL', cells)] + [(g, cells[cells['reporter_group'] == g]) for g in groups]

                for name, cc in strata:
                    tot = cc.groupby('wl')[['total_cases', col]].sum()
                    if set(tot.index) != {0, 1}:
                        continue
                    ev1, n1 = tot.loc[1, col], tot.loc[1, 'total_cases']
                    ev0, n0 = tot.loc[0, col], tot.loc[0, 'total_cases']
                    counts = dict(events_wl=ev1, n_wl=n1, events_dm=ev0, n_dm=n0)

                    # --- MODEL 1: CRUDE ROR ---
                    # Crude reporting odds ratio with Haldane-Anscombe correction for zero cells
                    a, b, c, d = ev1, n1 - ev1, ev0, n0 - ev0
                    if min(a, b, c, d) == 0:
                        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
                    ror = (a * d) / (b * c)
                    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
                    crude = dict(target=target, scope=scope, outcome=label, odds_ratio=ror, ci_low=np.exp(np.log(ror) - 1.96 * se),
                                 ci_high=np.exp(np.log(ror) + 1.96 * se), p_value=2 * norm.sf(abs(np.log(ror) / se)), **counts)
                    if name == 'ALL':
                        rows.append(dict(method='Crude ROR', **crude))

                    # --- MODEL 2: ADJUSTED QUASI-BINOMIAL GLM ---
                    # Overall model adjusts for year and reporter; stratum models adjust for year only
                    method = 'Adjusted OR (year + reporter)' if name == 'ALL' else f"Stratum: {name}"
                    rhs = 'wl + C(year_f) + C(reporter_group)' if name == 'ALL' else 'wl + C(year_f)'
                    if min(ev1, ev0, n1 - ev1, n0 - ev0) == 0:
                        report.append(f"Adjusted GLM skipped for {target} / {scope} / {label} / {name}: zero events or non-events in a stratum.")
                        if name != 'ALL':
                            rows.append(dict(method=method, note='crude: adjusted model not estimable', **crude))
                        continue
                    try:
                        glm = smf.glm(f"{col} + fail ~ {rhs}", data=cc, family=sm.families.Binomial())
                        # Dispersion is floored at 1 so the CI is never artificially narrower than the plain binomial one
                        first = glm.fit()
                        phi = max(1.0, first.pearson_chi2 / first.df_resid) if first.df_resid > 0 else 1.0
                        fit = glm.fit(scale=phi)
                        ci = fit.conf_int().loc['wl']
                        rows.append(dict(target=target, scope=scope, outcome=label, method=method, odds_ratio=np.exp(fit.params['wl']),
                                         ci_low=np.exp(ci[0]), ci_high=np.exp(ci[1]), p_value=fit.pvalues['wl'], **counts))
                    except Exception as e:
                        report.append(f"GLM failed for {target} / {scope} / {label} / {name}: {e}")

                # --- HETEROGENEITY TEST (Likelihood Ratio Test) ---
                # Tests if the weight-loss effect differs significantly between reporter groups (Simpson's Paradox check)
                if len(groups) >= 2:
                    hc = cells[cells['reporter_group'].isin(groups)]
                    chk = hc.groupby(['reporter_group', 'wl'])[[col, 'fail']].sum()
                    if len(chk) == 2 * len(groups) and chk.min().min() > 0:
                        try:
                            full = smf.glm(f"{col} + fail ~ wl * C(reporter_group) + C(year_f)", data=hc, family=sm.families.Binomial()).fit()
                            red = smf.glm(f"{col} + fail ~ wl + C(reporter_group) + C(year_f)", data=hc, family=sm.families.Binomial()).fit()
                            phi = max(1.0, full.pearson_chi2 / full.df_resid) if full.df_resid > 0 else 1.0
                            p = chi2.sf((red.deviance - full.deviance) / phi, full.df_model - red.df_model)
                            rows.append(dict(target=target, scope=scope, outcome=label, method='Heterogeneity by reporter (LRT)', p_value=p, note=f"groups: {', '.join(groups)}"))
                        except Exception as e:
                            report.append(f"Heterogeneity test failed for {target} / {scope} / {label}: {e}")
                    else:
                        report.append(f"Heterogeneity test skipped for {target} / {scope} / {label}: a reporter-by-stratum cell has no events or no non-events.")

        # --- MODEL 3: SAME-CLASS INTERACTION MODEL ---
        # Tests if the weight-loss contrast is larger for the target than for comparator GLP-1s
        pool = df[df['drug'].isin([target] + args.comparators)]
        comp_wl = pool[(pool['drug'] != target) & (pool['indication'] == 'Weight Loss')]['total_cases'].sum()
        if args.comparators and comp_wl >= 30:
            cells = pool.groupby(['year', 'reporter_group', 'drug', 'indication'], as_index=False)[cols].sum()
            cells['wl'] = (cells['indication'] == 'Weight Loss').astype(int)
            cells['is_target'] = (cells['drug'] == target).astype(int)
            cells['year_f'] = cells['year'].astype(str)
            for col, label in outcomes.items():
                cells['fail'] = cells['total_cases'] - cells[col]
                # Every indication x drug-group cell needs both events and non-events
                grp = cells.groupby(['wl', 'is_target'])[[col, 'fail']].sum()
                if len(grp) < 4 or (grp.min().min() == 0):
                    report.append(f"Interaction GLM skipped for {target} / {label}: zero events or non-events in a cell.")
                    continue
                try:
                    glm = smf.glm(f"{col} + fail ~ wl * is_target + C(year_f) + C(reporter_group)", data=cells, family=sm.families.Binomial())
                    first = glm.fit()
                    phi = max(1.0, first.pearson_chi2 / first.df_resid) if first.df_resid > 0 else 1.0
                    fit = glm.fit(scale=phi)
                    ci = fit.conf_int()
                    for term, method in [('wl', 'Comparator GLP-1s only: WL vs DM'), ('wl:is_target', 'Interaction: extra WL effect for target')]:
                        rows.append(dict(target=target, scope='All products', outcome=label, method=method, odds_ratio=np.exp(fit.params[term]),
                                         ci_low=np.exp(ci.loc[term, 0]), ci_high=np.exp(ci.loc[term, 1]), p_value=fit.pvalues[term]))
                except Exception as e:
                    report.append(f"Interaction GLM failed for {target} / {label}: {e}")
        elif args.comparators:
            report.append(f"Interaction model skipped: only {int(comp_wl)} comparator weight-loss cases (need 30).\n")

    res = pd.DataFrame(rows).reindex(columns=out_cols)
    res['note'] = res['note'].fillna('')
    res.to_csv(args.output_csv, index=False)

    report.append("\nEstimates (odds ratios with 95% CI)")
    report.append("Crude ROR: pooled 2x2 table. Adjusted OR: quasi-binomial GLM (dispersion floored at 1) with year and reporter-group terms.")
    report.append("Stratum: same GLM within one reporter group, year only. Heterogeneity: likelihood-ratio test of weight-loss x reporter group.")
    report.append("Interaction: ratio of the weight-loss odds ratio for the target to that of the comparator GLP-1s.\n")
    if len(res):
        report.append(res.drop(columns=['target']).to_string(index=False, float_format=lambda x: f"{x:.3g}"))

        # Plain-language read of the primary adjusted estimate
        prim = res[(res['method'] == 'Adjusted OR (year + reporter)') & (res['scope'] == 'All products') & res['outcome'].str.contains('primary')]
        for _, r in prim.iterrows():
            if r['ci_low'] > 1:
                verdict = "higher"
            elif r['ci_high'] < 1:
                verdict = "lower"
            else:
                verdict = "not clearly different"
            report.append(f"\nPrimary outcome, {r['target']}: weight-loss reports have {verdict} odds of a medication-error report than diabetes reports "
                          f"(adjusted OR {r['odds_ratio']:.2f}, 95% CI {r['ci_low']:.2f}-{r['ci_high']:.2f}).")
    else:
        report.append("No estimates could be computed (check that both indication strata have cases).")

    with open(args.output_report, 'w') as f:
        f.write("\n".join(report))
    print(f"Statistics written to {args.output_csv} and {args.output_report}")

if __name__ == "__main__":
    main()
