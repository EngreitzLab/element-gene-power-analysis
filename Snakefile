# config file containing parameters
configfile: "config/config.yml"

# Include all rules for processing validation datasets
include: "rules/sceptre_power_analysis.smk"

# Perform all analyses to output benchmarked datasets
rule all:
  input:
    expand(
      "results/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
      sample=config["samples"],
      effect_size=config["effect_sizes"]
    )
