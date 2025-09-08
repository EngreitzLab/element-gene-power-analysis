# Snakemake rules to run power analysis using sceptre


# Run the power simulation with sceptre for each split
rule sceptre_power_analysis:
    input:
        gene_grna_group_pairs_split="results/{sample}/pair_splits/gene_grna_group_pairs_{split}.txt",
        final_sceptre_object="results/{sample}/differential_expression/final_sceptre_object.rds",
        grna_groups_table="results/{sample}/grna_groups_table.rds",
        perturb_sce="results/{sample}/perturb_sce.rds",
    output:
        power_analysis_output="results/{sample}/power_analysis/effect_size_{effect_size}/power_analysis_output_{split}.tsv",
    params:
        reps=config["num_replicates"],
    log:
        "results/logs/{sample}/sceptre_power_analysis_es_{effect_size}_split{split}.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="8G",
        time="1:00:00",
    script:
        "../R/sceptre_power_analysis.R"


# Create and split the discovery pairs file
# N_BATCHES = config["num_batches"]


rule split_target_response_pairs:
    input:
        gene_grna_group_pairs="results/{sample}/gene_grna_group_pairs.rds",
    output:
        splits=lambda wildcards: [
            f"results/{wildcards.sample}/pair_splits/gene_grna_group_pairs_{split}.txt"
            for split in range(1, config["num_batches"] + 1)
        ],
        # Default to 100 splits - this only works if every dataset has more than 100 genes, which is true (error will be thrown if not)
    params:
        batches=config["num_batches"],
    log:
        "results/{sample}/logs/split_target_response_pairs.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="24G",
        time="2:00:00",
    script:
        "../R/split_target_response_pairs.R"


# Run sceptre differential expression with "union"
rule sceptre_differential_expression:
    input:
        sceptre_diffex_input="results/{sample}/sceptre_diffex_input.rds",
    output:
        discovery_results="results/{sample}/results_run_discovery_analysis.rds",
        final_sceptre_object="results/{sample}/final_sceptre_object.rds",
    log:
        "results/logs/{sample}/sceptre_differential_expression.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="32G",
        time="12:00:00",
    script:
        "../R/sceptre_differential_expression.R"


# Create the sce object from sceptre object for simulations
rule create_sce:
    input:
        final_sceptre_object="results/{sample}/final_sceptre_object.rds",
    output:
        perturb_sce="results/{sample}/perturb_sce.rds",
    log:
        "results/logs/{sample}/create_sce.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="64G",
        time="4:00:00",
    script:
        "../R/create_sce_object.R"


# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis:
    input:
        splits=lambda wildcards: [
            f"results/{wildcards.sample}/power_analysis/effect_size_{wildcards.effect_size}/power_analysis_output_{split}.tsv"
            for split in range(1, config["num_batches"] + 1)
        ],
    output:
        combined_power_analysis_output="results/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
    log:
        "results/logs/{sample}/combine_sceptre_power_analysis_es_{effect_size}.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="32G",
        time="2:00:00",
    script:
        "../R/combine_sceptre_power_analysis.R"


# Compute the power from the power simulations
rule compute_power_from_simulations:
    input:
        combined_power_analysis_output="results/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
        discovery_results="results/{sample}/differential_expression/results_run_discovery_analysis.rds",
    output:
        power_analysis_results="results/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv",
    log:
        "results/logs/{sample}/compute_power_from_simulations_es_{effect_size}.log",
    conda:
        "../envs/power_analysis.yml"
    resources:
        mem="24G",
        time="1:00:00",
    script:
        "../R/compute_power_from_simulations.R"
