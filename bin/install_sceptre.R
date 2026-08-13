#!/usr/bin/env Rscript
#
# Install sceptre from a pinned commit into the project-local R library.
#
# sceptre is not available on conda-forge or bioconda, so it cannot be captured in
# pixi.lock. The pipeline depends on *unexported* S4 slots of `sceptre_object`
# (see bin/check_sceptre_api.R), which means an unpinned install can silently change
# pipeline behaviour. The pin lives in pixi.toml ([activation.env] SCEPTRE_SHA) so that
# the environment definition and the pin stay in one place.
#
# Usage:  pixi run setup
#         Rscript bin/install_sceptre.R          # inside an activated environment

lib <- Sys.getenv("R_LIBS_USER")
sha <- Sys.getenv("SCEPTRE_SHA")
ref <- Sys.getenv("SCEPTRE_REF", unset = "(unnamed)")

if (!nzchar(lib)) {
  stop("R_LIBS_USER is not set. Run this via `pixi run setup` so the pixi environment is active.")
}
if (!nzchar(sha)) {
  stop("SCEPTRE_SHA is not set. It is defined in pixi.toml under [activation.env]; ",
       "run this via `pixi run setup`.")
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# Install into the project-local library only, so nothing is written into the shared
# pixi environment and a re-solve of the environment cannot leave a stale sceptre behind.
.libPaths(c(lib, .libPaths()))

installed <- tryCatch(as.character(utils::packageVersion("sceptre")), error = function(e) NA_character_)
installed_sha <- tryCatch(
  utils::packageDescription("sceptre")$RemoteSha,
  error = function(e) NULL
)

if (!is.na(installed) && identical(installed_sha, sha)) {
  message("sceptre ", installed, " (", substr(sha, 1, 8), ") is already installed in ", lib, ".")
  message("Nothing to do. Delete ", file.path(lib, "sceptre"), " to force a reinstall.")
  quit(save = "no", status = 0)
}

message("Installing sceptre ", ref, " (", substr(sha, 1, 8), ") into ", lib, " ...")

remotes::install_github(
  repo = paste0("Katsevich-Lab/sceptre@", sha),
  lib = lib,
  upgrade = "never",      # every dependency is already pinned by pixi.lock
  dependencies = FALSE,   # ditto: do not let remotes pull CRAN builds over the conda ones
  build_vignettes = FALSE,
  quiet = FALSE
)

version <- as.character(utils::packageVersion("sceptre", lib.loc = lib))
message("Installed sceptre ", version, ".")
message("Next: pixi run check-api")
