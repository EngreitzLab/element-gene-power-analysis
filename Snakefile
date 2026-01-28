#Configuration file for the workflow
configfile: "config/config.yml"

# Include all rules from the separate file
include: "rules/sceptre_power_analysis.smk"

# Define the final target file(s)
rule all:
  input:
    expand(
     "results/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv",
     sample=config["samples"],
     effect_size=config["effect_sizes"]
    )

