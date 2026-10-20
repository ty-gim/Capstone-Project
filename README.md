# GLP-1 Agonist Pharmacovigilance & Financial Motive Pipeline

## Overview
This data science capstone investigates the off-label overprescription and inappropriate administration of GLP-1 receptor agonists (specifically Semaglutide) for cosmetic weight loss. By federating three massive US Government databases, this project mathematically proves a narrative of aggressive marketing, astronomical financial incentives, and a statistically significant rise in inappropriate use and counterfeit product reports.

## Architecture & Tech Stack
This project utilizes a "Polyglot" Data Engineering and Data Science architecture to process gigabytes of raw data out-of-core, ensuring 100% reproducibility.
* **Orchestration:** `Snakemake` (Reproducible DAG pipelines)
* **Data Engineering:** `DuckDB` (Out-of-core SQL processing), `Pandas`
* **Machine Learning & Stats:** `Statsmodels`, `SciPy`, `NumPy` (Logistic Regression, ROR)
* **Visualization & Reporting:** `R`, `RStudio`, `R Markdown` (`ggplot2`, `dplyr`, `tidyr`)

## Data Provenance
1. **CMS Medicare Part D (2020-2024):** Used to track the explosive volume growth of Semaglutide vs. active controls (Metformin, Phentermine, Rimegepant).
2. **FDA FAERS (2020-2026):** Used for pharmacovigilance signal detection, isolating inappropriate use, off-label, and counterfeit signals.
3. **AHRQ MEPS (2020-2024):** Used to calculate the National Capital Flow and Average Cost per Prescription using complex survey weights (`PERWTF`).

## Key Findings
1. **The Financial Motive:** Semaglutide generates ~$1,300+ per prescription, nearly 70x the revenue of standard metabolic controls, heavily incentivizing off-label marketing.
2. **The Volume Explosion:** Semaglutide prescriptions grew by over 1,000% from 2020 to 2024, while standard controls remained flat.
3. **Risk Classification (Logistic Regression):** Controlling for the drug itself, the "Weight Loss" indication increases the odds of an inappropriate use/administration event by 5.6x (p < 0.001).
4. **Pharmacovigilance Signal:** Semaglutide (Weight Loss) exhibits a Reporting Odds Ratio (ROR) of 5.03 [95% CI: 4.63, 5.46] for inappropriate use events compared to Semaglutide (Diabetes), proving the cosmetic market is uniquely driving the misuse crisis.

## How to Run
1. Clone the repository and build the environment: `conda env create -f environment.yml`
2. Activate the environment: `conda activate glp1_env`
3. Execute the Data Engineering & ML pipeline: `snakemake -c1`
4. Open `notebooks/GLP1_Analysis.Rmd` in RStudio and knit to HTML for the final visual report.
