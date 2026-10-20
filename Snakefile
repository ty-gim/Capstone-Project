configfile: "config.yaml"

rule all:
    input:
        "data/processed/cms_volume_output.csv",
        "data/processed/faers_top5_visuals.csv",
        "data/processed/faers_full_ml_data.csv",
        "data/processed/meps_financial_output.csv",
        "data/processed/ml_insights_report.txt"

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
        visuals = "data/processed/faers_top5_visuals.csv",
        ml_data = "data/processed/faers_full_ml_data.csv"
    benchmark:
        "data/processed/benchmarks/analyze_faers.tsv"
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
            --output_visuals {output.visuals} \
            --output_ml {output.ml_data}
        """

rule analyze_meps:
    input:
        meps_dir = config["meps_folder"],
        config_file = "config.yaml"
    output:
        "data/processed/meps_financial_output.csv"
    benchmark:
        "data/processed/benchmarks/analyze_meps.tsv"
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
            --output {output}
        """

rule run_ml_models:
    input:
        cms = "data/processed/cms_volume_output.csv",
        faers_ml = "data/processed/faers_full_ml_data.csv"
    output:
        report = "data/processed/ml_insights_report.txt"
    benchmark:
        "data/processed/benchmarks/run_ml_models.tsv"
    shell:
        """
        python src/run_ml_models.py \
            --cms_file {input.cms} \
            --faers_file {input.faers_ml} \
            --output_report {output.report}
        """
