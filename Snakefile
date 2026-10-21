configfile: "config.yaml"

rule all:
    input:
        "reports/GLP1_Analysis.html"

rule analyze_cms:
    input:
        cms_dir = config["cms_folder"],
        config_file = "config.yaml"
    output:
        "data/processed/cms_volume_output.csv"
    benchmark:
        "data/processed/benchmarks/analyze_cms.tsv"
    params:
        targets = config["target_drugs"],
        controls = config["control_drugs"],
        files = config["cms_files"]
    shell:
        """
        python src/analyze_cms.py \
            --cms_folder {input.cms_dir} \
            --cms_files {params.files} \
            --targets {params.targets} \
            --controls {params.controls} \
            --config_file {input.config_file} \
            --output {output}
        """

rule analyze_faers:
    input:
        faers_dir = config["faers_folder"],
        config_file = "config.yaml"
    output:
        "data/processed/faers_cases.csv"
    benchmark:
        "data/processed/benchmarks/analyze_faers.tsv"
    log:
        "data/processed/logs/analyze_faers.log"
    params:
        targets = config["target_drugs"],
        controls = config["control_drugs"],
        quarters = config["faers_quarters"]
    shell:
        """
        python src/analyze_faers.py \
            --targets {params.targets} \
            --controls {params.controls} \
            --faers_folder {input.faers_dir} \
            --faers_quarters {params.quarters} \
            --config_file {input.config_file} \
            --output {output} > {log} 2>&1
        """

rule analyze_meps:
    input:
        meps_dir = config["meps_folder"],
        config_file = "config.yaml"
    output:
        "data/processed/meps_financial_output.csv"
    benchmark:
        "data/processed/benchmarks/analyze_meps.tsv"
    log:
        "data/processed/logs/analyze_meps.log"
    params:
        targets = config["target_drugs"],
        controls = config["control_drugs"],
        files = config["meps_files"]
    shell:
        """
        python src/analyze_meps.py \
            --meps_folder {input.meps_dir} \
            --meps_files {params.files} \
            --targets {params.targets} \
            --controls {params.controls} \
            --config_file {input.config_file} \
            --output {output} > {log} 2>&1
        """

rule run_stats_models:
    input:
        faers = "data/processed/faers_cases.csv"
    output:
        csv = "data/processed/stats_results.csv",
        report = "data/processed/stats_report.txt"
    benchmark:
        "data/processed/benchmarks/run_stats_models.tsv"
    params:
        targets = config["target_drugs"],
        comparators = config["glp1_comparators"],
        brand_scopes = config["stats"]["brand_scopes"],
        exclude_products = config["stats"]["exclude_products"]
    shell:
        """
        python src/run_stats_models.py \
            --faers_file {input.faers} \
            --targets {params.targets} \
            --comparators {params.comparators} \
            --brand_scopes {params.brand_scopes} \
            --exclude_products {params.exclude_products} \
            --output_csv {output.csv} \
            --output_report {output.report}
        """

rule render_report:
    input:
        cms = "data/processed/cms_volume_output.csv",
        faers = "data/processed/faers_cases.csv",
        meps = "data/processed/meps_financial_output.csv",
        stats = "data/processed/stats_results.csv",
        rmd = "reports/GLP1_Analysis.Rmd"
    output:
        "reports/GLP1_Analysis.html"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd}', output_file='GLP1_Analysis.html', output_dir='reports', quiet=TRUE)"
        """
