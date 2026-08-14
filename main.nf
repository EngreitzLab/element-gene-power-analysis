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

    // ---- inputs ---------------------------------------------------------------------------
    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.sceptre_object) {
                error "samplesheet ${params.samplesheet} needs columns 'sample' and 'sceptre_object'; got: ${row}"
            }
            [ [id: row.sample], file(row.sceptre_object, checkIfExists: true) ]
        }

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
