// Step 2c -- merge the per-replicate null-model chunks into the bundle the simulation consumes.
//
// merge_null_models.R errors on overlapping --rep-offset ranges rather than silently keeping one of
// them, and on a replicate count that does not match --reps. Both matter here: a missing chunk would
// otherwise surface as a simulation that quietly tests some replicates against no null model.
//
// An entry is 11 named coefficients plus a theta scalar, so 100 replicates x 272 genes is a few
// hundred KB -- the bundle is small enough to hand to every simulation task.

process MERGE_NULL_MODELS {
    tag "${meta.id}"

    publishDir "${params.outdir}/${meta.id}/prepared", mode: params.publish_mode

    input:
    tuple val(meta), path(chunks, stageAs: 'chunks/*')
    val  expected_reps

    output:
    tuple val(meta), path('null_precomputations.rds'), emit: null_models

    script:
    """
    n_chunks=\$(ls chunks/*.rds | wc -l)
    echo "merging \${n_chunks} chunk(s) covering ${expected_reps} replicates"

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/merge_null_models.R \\
            --inputs 'chunks/*.rds' \\
            --reps ${expected_reps} \\
            --out null_precomputations.rds
    """

    stub:
    """
    touch null_precomputations.rds
    """
}
