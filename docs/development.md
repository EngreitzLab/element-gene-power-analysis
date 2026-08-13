---
title: Development
nav_order: 6
---

# Development

## Environment

[pixi](https://pixi.sh) manages everything except sceptre.

```sh
pixi install          # runtime environment from pixi.lock
pixi run setup        # install sceptre from the pinned commit
pixi run check-api    # verify the pin against the internals this pipeline uses
pixi run lint         # formatting and convention checks
pixi run install-hooks
```

Three environments, so the one every task activates stays small:

| Environment | Contains | Why separate |
|---|---|---|
| default | `r-base`, `r-optparse`, `nextflow`, and the packages sceptre reaches on our code path (`r-matrix`, `r-bh`, `r-rcpp`, `r-dplyr`, `r-data.table`, `r-purrr`, `r-crayon`) | activated by every task |
| `build` | `r-remotes`, compilers, `r-ggplot2`, `r-cowplot`, `r-scales` | only needed to compile sceptre |
| `dev` | `r-testthat` | only needed to run tests |

### Why `ggplot2` is a build-only dependency

`R CMD INSTALL` enforces a package's DESCRIPTION `Imports` **even for packages absent from its
NAMESPACE**. sceptre's NAMESPACE imports only `BH`, `Matrix` and `Rcpp`, but its DESCRIPTION lists
`ggplot2`, `cowplot`, `scales`, `dplyr`, `data.table`, `purrr` and `crayon`, so all of them must be
present to build it. At *runtime* only the ones actually called matter — and `ggplot2`, `cowplot` and
`scales` appear solely in `plotting_functions.R`, `qq_plot_helpers.R` and
`write_outputs_to_directory()`, none of which this pipeline calls. Hence the split.

`purrr`, by contrast, *is* needed at runtime: `purrr::flatten()` is called in the precomputation
path that `run_discovery_analysis()` goes through.

### No Bioconductor

The pipeline once used `SingleCellExperiment`, but only as a container — a per-gene table, a per-cell
table, two sparse perturbation matrices, and subsetting. Nothing touched `rowRanges` or any genomic
range. `bin/lib/sim_input.R` replaces it with a plain list, which removes sixteen Bioconductor
packages and leaves `nextflow` as the only bioconda dependency.

It also avoids a hard setup failure: `bioconductor-genomeinfodbdata` installs its data through a
conda post-link script, which pixi skips by default, so `GenomeInfoDb` — and therefore
`SingleCellExperiment` — could not be loaded at all without every user first running
`pixi config set --local run-post-link-scripts insecure`.

## sceptre is pinned, and why that matters

sceptre is not on conda-forge or bioconda, so it is installed from a pinned commit
(`SCEPTRE_REF` / `SCEPTRE_SHA` in `pixi.toml`) into a project-local library at `.pixi/rlibs`.

The pin is not housekeeping. This pipeline reaches into sceptre's **unexported S4 slots**: it swaps
the response matrix for simulated counts, narrows the discovery pairs to one target at a time, and
reads the cached negative-binomial precomputations. None of that is a documented interface, so a
sceptre release could change it with no deprecation warning and the pipeline would keep running while
producing wrong numbers.

`bin/check_sceptre_api.R` is the guard. It asserts the presence of every slot the pipeline reads:

`@response_matrix`, `@grna_matrix`, `@covariate_data_frame`, `@covariate_matrix`,
`@grna_target_data_frame`, `@discovery_pairs_with_info`, `@initial_grna_assignment_list`,
`@grna_assignments`, `@cells_in_use`, `@response_precomputations`, `@functs_called`

Run it after any change to the pin. All slot access is confined to `bin/lib/sceptre_io.R`, so a
sceptre upgrade breaks one file rather than five.

Two behaviours of sceptre worth knowing before touching the simulation:

- **`load_row()` dispatches on `odm` and `dgRMatrix` with no `else` branch.** Handing it a dense
  matrix returns `NULL` silently instead of erroring. The simulated counts must be converted to
  `dgRMatrix` before being assigned to `@response_matrix`.
- **With `run_permutations = FALSE`** the conditional-randomisation path reads
  `response_matrix@j/@p/@x` directly, which only a `dgRMatrix` has.

## Repository layout

```
bin/                      standalone executables, one per pipeline step
bin/lib/                  shared code, sourced by script-relative path
  cli.R                   argument parsing, TSV I/O, seeding
  sim_input.R             the simulation-input container and its accessors
  pert_input.R            per-target cell selection
  simulate.R              count simulation and guide-level effect sizes
  sceptre_io.R            the only file that touches sceptre internals
config/config.yml         pipeline parameters
assets/samplesheet.csv    example samplesheet
docs/                     this documentation (published to GitHub Pages)
.githooks/pre-commit      formatting and convention checks
```

`bin/lib/*.R` is sourced via a path resolved from `--file=` in `commandArgs()`, so every script works
identically whether run directly or staged onto `PATH` by a workflow engine.

## Conventions, enforced by the hook

`pixi run install-hooks` points `core.hooksPath` at `.githooks`. On commit, staged files are
**auto-fixed** for:

- CRLF and lone CR → LF
- tabs → spaces, at the width conventional for the language: **2 for R and YAML**, 4 otherwise
  (source files only; `.tsv`/`.csv` are data and `Makefile` needs tabs)
- trailing whitespace
- missing final newline

and **rejected** for things that cannot be fixed without changing meaning:

- any file over **512 KB** — pipeline inputs and outputs belong outside git; `results/` and `*.rds`
  are gitignored, and this is the backstop
- filenames with uppercase in the stem (the extension is exempt, so `compute_power.R` is fine)
- camelCase or PascalCase R identifiers

The identifier check deliberately allows `SCREAMING_SNAKE` constants (`PERT_LEVELS`) and dotted S3
methods (`print.sim_input`), both of which are correct R style. `pixi run lint` runs the same checks
across every tracked file.

R code is indented with **2 spaces**, following tidyverse and Google R style. YAML uses 2 as well
(and forbids tabs outright); shell, Groovy/Nextflow, Python and JSON use 4. The hook expands tabs to
the matching width so a stray tab does not end up fighting the surrounding style, but it does not
reindent existing code.

## Documentation

`docs/` is published to GitHub Pages by `.github/workflows/pages.yml`. Pages are ordinary markdown
with `title` and `nav_order` front matter.

The theme is pinned to `just-the-docs@v0.3.3` because GitHub's `jekyll-build-pages` action uses the
`github-pages` gem (Jekyll 3.9); just-the-docs v0.4+ requires Jekyll 4 and will not build. If the
Pages build ever breaks, replacing `remote_theme` with `theme: jekyll-theme-primer` in
`docs/_config.yml` is a zero-dependency fallback.

## Known gaps

- **Workflow orchestration.** The five steps are complete and run standalone; a Nextflow workflow
  with a SLURM profile is not yet in the repository. `config/config.yml` already carries the
  parameters it will consume.
- **Two-stage replicate allocation** is documented but not orchestrated. The scripts support it
  today via `--rep-offset`.
- **Equivalence with the pre-refactor pipeline** is proven for the upstream statistics — size
  factors, normalised means and raw means are bit-identical, and dropping `@grna_matrix` was verified
  to leave discovery results unchanged — but the simulation itself has not been compared draw-for-draw
  against the old code. Seeding moved inside the replicate loop, so that comparison has to be
  distributional rather than exact.
- **The legacy `R/` directory** is still present so the old Snakemake pipeline remains runnable for
  comparison. `pixi run lint` reports its non-conforming identifiers; both go away when it is removed.
