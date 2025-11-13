# element-gene-power-analysis
“This repository implements a SCEPTRE-based power analysis that tests element–gene pairs in single-cell CRISPR screens, returning per-pair power, expected discoveries, and minimal sample sizes given effect sizes, detection rates, and α/FDR.”

The original code for the SCEPTRE method is available at https://github.com/EngreitzLab/DC_TAP_Paper but it was really tailored for the specific analysis in the paper. This repository generalizes the code and provides a more user-friendly interface for performing power analysis using SCEPTRE.

At the moment it should work but it is still a work in progress.

```sh
snakemake --use-conda all --conda-frontend conda --cores 1
```

Folder structure and file requirements:
 - The following input files needs to be places in `results/<sample>/`:
    - `final_sceptre_object.rds`
    - `grna_groups_table.rds`
    - `gene_grna_group_pairs.rds`



Example of `grna_groups_table.rds`
```R
                    grna_id          grna_target
1  240905-HPC-screen-guide-1 chr1:7885554-7885848
2 240905-HPC-screen-guide-10 chr1:7885554-7885848
3 240905-HPC-screen-guide-11 chr1:7885554-7885848
```


Example of `gene_grna_group_pairs.rds`
```R
   grna_target        response_id
1 chr1:7885554-7885848 ENSG00000116288.13
2 chr1:7902783-7902990 ENSG00000116288.13
3 chr1:7913116-7913357 ENSG00000116288.13
```




`config/config.yaml` needs to be filled with the appropriate parameters.
