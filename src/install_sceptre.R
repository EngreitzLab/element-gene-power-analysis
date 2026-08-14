#!/usr/bin/env Rscript
#
# Install sceptre from a pinned commit into the project-local R library, applying the
# local patches in patches/ first.
#
# sceptre is not available on conda-forge or bioconda, so it cannot be captured in
# pixi.lock. The pipeline depends on *unexported* S4 slots of `sceptre_object`
# (see bin/check_sceptre_api.R), which means an unpinned install can silently change
# pipeline behaviour. The pin lives in pixi.toml ([activation.env] SCEPTRE_SHA) so that
# the environment definition and the pin stay in one place.
#
# The patches add optional arguments that let this pipeline hand sceptre work it has already
# done, rather than having it recomputed once per replicate. Each is written against the
# pinned SHA and each has a `NULL` default, so unpatched sceptre behaves identically. See
# patches/*.patch for what and why, and sceptre-upstream-pr.md for the plan to upstream them.
#
# This is a source install: it compiles sceptre's C++, so run it in a job or an `sh_dev`
# shell, never on a login node.
#
# Usage:  pixi run setup
#         Rscript bin/install_sceptre.R          # inside an activated environment
#
# Environment:
#   R_LIBS_USER      target library (set by pixi)
#   SCEPTRE_SHA      commit to install (set by pixi)
#   SCEPTRE_REF      human-readable tag for messages (set by pixi)
#   SCEPTRE_TARBALL  optional path to an already-downloaded source tarball for the pinned
#                    SHA, for installing on a node with no outbound network
#   SCEPTRE_PATCH_DIR  optional directory to take patches from instead of <root>/patches.
#                    Pointing it at an empty directory builds a *stock* sceptre, which is how
#                    step 11 isolates the effect of the sceptre version from the effect of the
#                    patch. Combine it with R_LIBS_USER to install alongside, not over, the
#                    pipeline's own library.

lib <- Sys.getenv("R_LIBS_USER")
sha <- Sys.getenv("SCEPTRE_SHA")
ref <- Sys.getenv("SCEPTRE_REF", unset = "(unnamed)")
tarball_override <- Sys.getenv("SCEPTRE_TARBALL", unset = "")

if (!nzchar(lib)) {
  stop("R_LIBS_USER is not set. Run this via `pixi run setup` so the pixi environment is active.")
}
if (!nzchar(sha)) {
  stop("SCEPTRE_SHA is not set. It is defined in pixi.toml under [activation.env]; ",
       "run this via `pixi run setup`.")
}

# Project root, so patches/ is found whether this is run as `Rscript bin/install_sceptre.R`
# or from an activated environment elsewhere.
project_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
  }
  root <- Sys.getenv("PIXI_PROJECT_ROOT", unset = "")
  if (nzchar(root)) normalizePath(root) else normalizePath(".")
})

patch_dir_override <- Sys.getenv("SCEPTRE_PATCH_DIR", unset = "")
patch_dir <- if (nzchar(patch_dir_override)) {
  patch_dir_override
} else {
  file.path(project_root, "patches")
}
if (!dir.exists(patch_dir)) {
  stop("Patch directory does not exist: ", patch_dir, ".")
}
patch_files <- sort(list.files(patch_dir, pattern = "[.]patch$", full.names = TRUE))

# The patches are part of the pin: an install without them has a different API surface, so
# refusing to proceed is safer than installing a sceptre the pipeline cannot use.
#
# An explicit SCEPTRE_PATCH_DIR is exempt. Pointing it at an empty directory is an unambiguous
# request for a stock build -- there is no way to do that by accident -- and it is what step 11
# needs in order to attribute a change in results to the sceptre version rather than the patch.
if (length(patch_files) == 0) {
  if (!nzchar(patch_dir_override)) {
    stop("No patches found in ", patch_dir, ". The pipeline requires the patched sceptre ",
         "(see patches/ and bin/check_sceptre_api.R). If the patches have been upstreamed, ",
         "remove this check together with the pin bump.")
  }
  message("SCEPTRE_PATCH_DIR=", patch_dir, " contains no patches: building a STOCK sceptre. ",
          "This build will not satisfy `pixi run check-api`.")
}

# Fingerprint the patches as well as the SHA. The SHA alone cannot tell a patched install
# from an unpatched one, nor an install made before a patch was edited from one made after.
# "none" rather than "" when there are no patches: an empty DCF field would make DESCRIPTION
# unparseable, and it also lets check_sceptre_api.R tell a stock build from a missing field.
patch_stamp <- if (length(patch_files) == 0) {
  "none"
} else {
  paste0(basename(patch_files), "@", tools::md5sum(patch_files), collapse = "; ")
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# Install into the project-local library only, so nothing is written into the shared
# pixi environment and a re-solve of the environment cannot leave a stale sceptre behind.
.libPaths(c(lib, .libPaths()))

installed <- tryCatch(as.character(utils::packageVersion("sceptre")), error = function(e) NA_character_)
# packageDescription() returns a bare NA (with a warning), not NULL, when the package is absent, so
# `$` on it would error rather than give NULL. That only shows up when installing into an empty
# library -- which is exactly what SCEPTRE_PATCH_DIR + R_LIBS_USER does for the step 11 comparison.
installed_desc <- suppressWarnings(
  tryCatch(utils::packageDescription("sceptre"), error = function(e) NULL)
)
if (!is.list(installed_desc)) {
  installed_desc <- NULL
}
installed_sha <- installed_desc$RemoteSha
installed_stamp <- installed_desc$LocalPatches

if (!is.na(installed) && identical(installed_sha, sha) && identical(installed_stamp, patch_stamp)) {
  message("sceptre ", installed, " (", substr(sha, 1, 8), " + ", length(patch_files),
          " patch(es)) is already installed in ", lib, ".")
  message("Nothing to do. Delete ", file.path(lib, "sceptre"), " to force a reinstall.")
  quit(save = "no", status = 0)
}
if (!is.na(installed) && identical(installed_sha, sha) && !identical(installed_stamp, patch_stamp)) {
  message("sceptre ", installed, " is installed at the right SHA but the patches have changed. ",
          "Reinstalling.")
}

source(file.path(project_root, "lib", "apply_patch.R"))

message("Installing sceptre ", ref, " (", substr(sha, 1, 8), ") into ", lib, " ...")

work_dir <- tempfile("sceptre-src-")
dir.create(work_dir)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

# --- obtain the source at the pinned SHA ---------------------------------------------------
tarball <- file.path(work_dir, "sceptre.tar.gz")
if (nzchar(tarball_override)) {
  if (!file.exists(tarball_override)) {
    stop("SCEPTRE_TARBALL is set to ", tarball_override, ", which does not exist.")
  }
  message("  using SCEPTRE_TARBALL=", tarball_override)
  file.copy(tarball_override, tarball)
} else {
  url <- paste0("https://codeload.github.com/Katsevich-Lab/sceptre/tar.gz/", sha)
  message("  downloading ", url)
  status <- utils::download.file(url, destfile = tarball, mode = "wb", quiet = TRUE)
  if (!identical(status, 0L) || !file.exists(tarball)) {
    stop("Failed to download the sceptre source from ", url, ". If this node has no outbound ",
         "network, download the tarball elsewhere and set SCEPTRE_TARBALL to its path.")
  }
}

utils::untar(tarball, exdir = work_dir)
src_dir <- file.path(work_dir, paste0("sceptre-", sha))
if (!dir.exists(src_dir)) {
  # A tarball supplied via SCEPTRE_TARBALL may unpack under a different top-level name.
  candidates <- setdiff(list.dirs(work_dir, recursive = FALSE), character(0))
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(candidates) != 1) {
    stop("Expected exactly one unpacked sceptre source directory in ", work_dir,
         " but found ", length(candidates), ".")
  }
  src_dir <- candidates[[1]]
}

# --- apply the local patches ---------------------------------------------------------------
# Each patch is written against SCEPTRE_SHA, so a mismatch here means the pin moved out from
# under it. apply_unified_diff() errors on the first context mismatch rather than sliding the
# hunk to somewhere it happens to fit -- see bin/lib/apply_patch.R. Failing is the point: an
# unapplied patch is not a missing optimisation, it is a missing argument that the pipeline
# passes by name, and the error would not surface until the simulation was already running.
for (patch_file in patch_files) {
  applied <- apply_unified_diff(patch_file, root = src_dir)
  message("  applied ", basename(patch_file), ": ",
          paste0(applied$file, " (", applied$hunks, ")", collapse = ", "))
}

# --- record the provenance in DESCRIPTION --------------------------------------------------
# R CMD INSTALL copies DESCRIPTION verbatim, so these fields survive into the installed
# package and are what the short-circuit above reads on the next run.
desc_path <- file.path(src_dir, "DESCRIPTION")
desc_lines <- readLines(desc_path, warn = FALSE)
desc_lines <- desc_lines[!grepl("^(Remote[A-Za-z]+|LocalPatches):", desc_lines)]
writeLines(
  c(
    desc_lines,
    "RemoteType: github",
    "RemoteHost: api.github.com",
    "RemoteUsername: Katsevich-Lab",
    "RemoteRepo: sceptre",
    paste0("RemoteRef: ", sha),
    paste0("RemoteSha: ", sha),
    paste0("LocalPatches: ", patch_stamp)
  ),
  desc_path
)

# --- install -------------------------------------------------------------------------------
utils::install.packages(
  src_dir,
  lib = lib,
  repos = NULL,
  type = "source",
  dependencies = FALSE,  # every dependency is already pinned by pixi.lock
  INSTALL_opts = "--no-multiarch"
)

version <- as.character(utils::packageVersion("sceptre", lib.loc = lib))
message("Installed sceptre ", version, " with ", length(patch_files), " local patch(es).")
message("Next: pixi run check-api")
