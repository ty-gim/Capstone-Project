# Semaglutide utilization, cost and reporting patterns

Descriptive analysis of public data: do FAERS reports for weight-related semaglutide use differ from diabetes reports in how often they describe medication errors, and how did semaglutide volume and cost develop (Medicare Part D, MEPS).

This is hypothesis-generating. FAERS is a spontaneous reporting system, none of the sources measures marketing, and the analysis does not estimate incidence or causation.

## Data (place under `data/raw/`)
- CMS Medicare Part D Prescribers by Provider and Drug, one CSV per year
- FAERS quarterly ASCII files (DEMO, DRUG, REAC, INDI)
- MEPS prescribed medicines event files (.dta)

Paths, years, drugs and brand mappings are set in `config.yaml`.

## Run
```
conda env create -f environment.yml && conda activate glp1-capstone
snakemake --cores 4
```
Output: `reports/GLP1_Analysis.html`, plus `data/processed/stats_report.txt`.

## Pipeline
| Step | Script | Output |
|---|---|---|
| CMS volume | `src/analyze_cms.py` | `cms_volume_output.csv` |
| FAERS case counts | `src/analyze_faers.py` | `faers_cases.csv` |
| MEPS row-level survey data | `src/analyze_meps.py` | `meps_financial_output.csv` |
| Weight-loss vs diabetes contrasts | `src/run_stats_models.py` | `stats_results.csv`, `stats_report.txt` |
| Report | `reports/GLP1_Analysis.Rmd` | `GLP1_Analysis.html` |

## Method notes
- FAERS: latest version of each case (by `caseid`), primary suspect drug, first-received year within the quarter folders; combination products excluded; outcomes counted as distinct cases.
- Outcomes: medication error/counterfeit (primary), the same without overdose/underdose (sensitivity), GI events (calibration), off-label use (manipulation check only).
- Statistics: crude reporting odds ratio, quasi-binomial GLM adjusted for year and reporter group, and an interaction model against injectable GLP-1 comparators.
- MEPS: design-based means with `VARSTR`/`VARPSU`/`PERWT`.

## Limitations
FAERS cannot estimate incidence, and reporting is influenced by media attention, litigation and shortages. The indication is the reporter's stated reason, and weight loss does not mean off-label (Wegovy is approved for chronic weight management). Semaglutide exposure includes compounded products that cannot be separated from manufactured ones. Part D excludes weight-loss indications, MEPS expenditure is gross payment per fill, and no source contains a marketing measure.
