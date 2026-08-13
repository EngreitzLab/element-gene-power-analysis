## Apply a unified diff, in pure R.
##
## Used by bin/install_sceptre.R to patch the pinned sceptre source before compiling it.
##
## WHY NOT `patch`
##
## Everything on this cluster runs inside a stock ubuntu:22.04 container (see
## workflow/slurm_executor/config.sh on why), and that image ships neither `patch` nor `git`.
## Adding a package to the image would mean maintaining a derived image for one command, so the
## ~90 lines below stand in for it.
##
## This applier is deliberately *stricter* than `patch(1)`:
##
##   * no fuzz and no offset search -- a hunk must match at exactly the line it claims. `patch`
##     helpfully slides a hunk to where it fits, which is the wrong behaviour here: our patches are
##     generated against one pinned commit, so a hunk that only matches somewhere else means the
##     pin moved and the patch needs regenerating, not relocating.
##   * every context and removed line must match the original byte for byte.
##   * no file creation, deletion or renaming -- only modification of files that exist.
##
## Anything unexpected is an error. A half-applied patch would produce a sceptre that compiles and
## silently lacks an argument the pipeline passes by name.

#' Apply a unified diff to a source tree.
#'
#' @param patch_path path to a unified diff. Any prose before the first `--- ` header is ignored,
#' so patches can carry an explanatory preamble.
#' @param root directory the paths in the diff are relative to, after stripping the leading
#' `a/` and `b/` components.
#' @return invisibly, a data frame of the files changed and the number of hunks applied to each
apply_unified_diff <- function(patch_path, root) {
  lines <- readLines(patch_path, warn = FALSE)
  n <- length(lines)

  # Skip the human-readable preamble. A `--- ` at column 0 is unambiguous here: inside a hunk body
  # a removed line would be `----` and a context line ` ---`, and the loop below only returns to
  # this state once a hunk has been fully consumed.
  i <- 1L
  while (i <= n && !startsWith(lines[[i]], "--- ")) i <- i + 1L
  if (i > n) {
    stop("No `--- ` file header found in ", patch_path, "; it does not look like a unified diff.",
         call. = FALSE)
  }

  changed_files <- character(0)
  hunk_counts <- integer(0)

  while (i <= n) {
    header_old <- lines[[i]]
    header_new <- if (i + 1L <= n) lines[[i + 1L]] else NA_character_
    if (is.na(header_new) || !startsWith(header_new, "+++ ")) {
      stop("Expected a `+++ ` header after `", header_old, "` in ", basename(patch_path), ".",
           call. = FALSE)
    }
    path_old <- patch_header_path(header_old, "--- ")
    path_new <- patch_header_path(header_new, "+++ ")
    if (!identical(path_old, path_new)) {
      stop("This applier only modifies files in place, but ", basename(patch_path), " maps `",
           path_old, "` to `", path_new, "`.", call. = FALSE)
    }

    target <- file.path(root, path_old)
    if (!file.exists(target)) {
      stop("Cannot apply ", basename(patch_path), ": ", path_old, " does not exist under ", root,
           ".", call. = FALSE)
    }
    original <- readLines(target, warn = FALSE)
    i <- i + 2L

    patched <- character(0)
    cursor <- 1L  # next line of `original` not yet copied or consumed
    n_hunks <- 0L

    while (i <= n && startsWith(lines[[i]], "@@")) {
      hunk <- parse_hunk_header(lines[[i]], path_old, patch_path)
      i <- i + 1L
      n_hunks <- n_hunks + 1L

      if (hunk$old_start < cursor) {
        stop("Hunks for ", path_old, " in ", basename(patch_path),
             " overlap or are out of order (hunk starts at line ", hunk$old_start,
             " but line ", cursor, " has already been consumed).", call. = FALSE)
      }
      if (hunk$old_start > cursor) {
        patched <- c(patched, original[cursor:(hunk$old_start - 1L)])
        cursor <- hunk$old_start
      }

      consumed_old <- 0L
      produced_new <- 0L
      while (i <= n && (consumed_old < hunk$old_len || produced_new < hunk$new_len)) {
        body <- lines[[i]]
        # A blank line in the diff body is a blank *context* line whose leading space some tools
        # strip. Treating it as context rather than erroring keeps patches robust to that.
        tag <- if (nchar(body) == 0L) " " else substr(body, 1L, 1L)
        text <- if (nchar(body) == 0L) "" else substr(body, 2L, nchar(body))

        if (identical(tag, "\\")) {
          # "\ No newline at end of file" -- readLines/writeLines normalise the final newline
          # anyway, so there is nothing to do.
          i <- i + 1L
          next
        }
        if (tag %in% c(" ", "-")) {
          if (cursor > length(original)) {
            stop("Hunk at line ", hunk$old_start, " of ", path_old, " runs past the end of the ",
                 "file (", length(original), " lines). The patch does not match this source.",
                 call. = FALSE)
          }
          if (!identical(original[[cursor]], text)) {
            stop("Context mismatch applying ", basename(patch_path), " to ", path_old,
                 " at line ", cursor, ".\n  expected: ", encodeString(text, quote = "\""),
                 "\n  found:    ", encodeString(original[[cursor]], quote = "\""),
                 "\n  The patch was written against a specific commit; if SCEPTRE_SHA has moved, ",
                 "regenerate it against the new source.", call. = FALSE)
          }
          if (identical(tag, " ")) {
            patched <- c(patched, text)
            produced_new <- produced_new + 1L
          }
          cursor <- cursor + 1L
          consumed_old <- consumed_old + 1L
        } else if (identical(tag, "+")) {
          patched <- c(patched, text)
          produced_new <- produced_new + 1L
        } else {
          stop("Unexpected line in a hunk for ", path_old, " of ", basename(patch_path), ": ",
               encodeString(body, quote = "\""), call. = FALSE)
        }
        i <- i + 1L
      }

      if (consumed_old != hunk$old_len || produced_new != hunk$new_len) {
        stop("Truncated hunk at line ", hunk$old_start, " of ", path_old, " in ",
             basename(patch_path), ": consumed ", consumed_old, "/", hunk$old_len,
             " original and produced ", produced_new, "/", hunk$new_len, " new lines.",
             call. = FALSE)
      }
    }

    if (n_hunks == 0L) {
      stop("No hunks found for ", path_old, " in ", basename(patch_path), ".", call. = FALSE)
    }
    if (cursor <= length(original)) {
      patched <- c(patched, original[cursor:length(original)])
    }
    writeLines(patched, target)
    changed_files <- c(changed_files, path_old)
    hunk_counts <- c(hunk_counts, n_hunks)

    # Advance to the next file section, stepping over the `diff ...` line diff(1) emits between
    # them. See the note above on why a bare `--- ` cannot be confused with hunk content here.
    while (i <= n && !startsWith(lines[[i]], "--- ")) i <- i + 1L
  }

  invisible(data.frame(file = changed_files, hunks = hunk_counts, stringsAsFactors = FALSE))
}


#' Extract the path from a `--- a/path` or `+++ b/path` header.
#'
#' Strips the one-component prefix that `diff -p1`-style patches carry, and any trailing tab-
#' separated timestamp.
patch_header_path <- function(header, prefix) {
  path <- substr(header, nchar(prefix) + 1L, nchar(header))
  path <- sub("\t.*$", "", path)
  path <- trimws(path)
  if (identical(path, "/dev/null")) {
    stop("This applier does not create or delete files, but the patch references /dev/null.",
         call. = FALSE)
  }
  if (!grepl("/", path, fixed = TRUE)) {
    stop("Expected a `a/`- or `b/`-prefixed path in the header `", header, "`.", call. = FALSE)
  }
  sub("^[^/]+/", "", path)
}


#' Parse an `@@ -old_start,old_len +new_start,new_len @@` hunk header.
#'
#' The counts are optional in unified diff format and default to 1 when omitted.
parse_hunk_header <- function(header, path, patch_path) {
  matched <- regmatches(
    header,
    regexec("^@@ -([0-9]+)(,([0-9]+))? \\+([0-9]+)(,([0-9]+))? @@", header)
  )[[1]]
  if (length(matched) == 0L) {
    stop("Malformed hunk header for ", path, " in ", basename(patch_path), ": ",
         encodeString(header, quote = "\""), call. = FALSE)
  }
  old_len <- if (nzchar(matched[[4]])) as.integer(matched[[4]]) else 1L
  new_len <- if (nzchar(matched[[7]])) as.integer(matched[[7]]) else 1L
  list(
    old_start = as.integer(matched[[2]]),
    old_len = old_len,
    new_start = as.integer(matched[[5]]),
    new_len = new_len
  )
}
