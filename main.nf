#!/usr/bin/env nextflow

// Element-gene power analysis.
//
//   nextflow run . -profile sherlock -params-file config/config.yml
//
// The DAG, and why it has the shape it has:
//
//   samplesheet -> PREPARE_SIM_INPUT ---+-> SPLIT_PAIRS -----------------+
//                                       |                               |
//                                       +-> FIT_NULL_MODELS (x reps)    |
//                                           -> MERGE_NULL_MODELS -------+
//                                                                       v
//                                       POWER_SIMULATION (split x effect size)
//                                            -> collectFile by (sample, effect size)
//                                            -> COMPUTE_POWER -> SUMMARIZE_POWER
//
// FIT_NULL_MODELS hangs off PREPARE_SIM_INPUT rather than off SPLIT_PAIRS because a gene's null
// model is independent of both the target and the effect size: it is fitted on a null simulation
// with no knockdown. One set of fits therefore serves every split and every effect size in a sweep
// -- 100 replicates, not 100 x however many effect sizes. Making it a sibling of SPLIT_PAIRS rather
// than a descendant is what expresses that.

nextflow.enable.dsl = 2

include { PREPARE_SIM_INPUT } from './modules/local/prepare_sim_input'

workflow {
    main:

    // ---- parameter checks -----------------------------------------------------------------
    //
    // Cheap to check here, expensive to discover 40 minutes into a 1,000-task run.
    if (!params.effect_sizes || params.effect_sizes.size() == 0) {
        error "effect_sizes is empty -- nothing to simulate."
    }
    if (params.num_replicates % params.reps_per_chunk != 0) {
        error "num_replicates (${params.num_replicates}) must be a multiple of reps_per_chunk " +
              "(${params.reps_per_chunk}); otherwise the last chunk is short and the power " +
              "denominators differ between pairs."
    }
    if (params.num_replicates % params.reps_per_null_chunk != 0) {
        error "num_replicates (${params.num_replicates}) must be a multiple of " +
              "reps_per_null_chunk (${params.reps_per_null_chunk}), or some replicate will have " +
              "no null model fitted for it."
    }

    // ---- inputs ---------------------------------------------------------------------------
    //
    // Samplesheet paths are resolved against the repository root, not the launch directory. The
    // committed samplesheets use repo-relative paths, so resolving against the launch directory
    // would make a run's validity depend on where it was started from.
    def resolve = { p ->
        def f = file(p)
        f.isAbsolute() ? f : file("${projectDir}/${p}")
    }

    ch_samples = Channel
        .fromPath(resolve(params.samplesheet), checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample?.trim() || !row.sceptre_object?.trim()) {
                error "samplesheet ${params.samplesheet} needs non-empty 'sample' and " +
                      "'sceptre_object' columns; got: ${row}"
            }
            def obj = resolve(row.sceptre_object.trim())
            if (!obj.exists()) {
                error "sample '${row.sample}': sceptre object not found at ${obj}"
            }
            [ [id: row.sample.trim()], obj ]
        }
        .ifEmpty { error "samplesheet ${params.samplesheet} has no rows." }

    // ---- step 1 ---------------------------------------------------------------------------
    PREPARE_SIM_INPUT(ch_samples)
}

workflow.onComplete {
    log.info(
        """
        ${workflow.success ? 'Completed' : 'FAILED'}: ${workflow.runName}
          duration : ${workflow.duration}
          outdir   : ${params.outdir}
          command  : ${workflow.commandLine}
        """.stripIndent()
    )
}
