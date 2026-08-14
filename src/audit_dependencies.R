#!/usr/bin/env Rscript
#
# Which packages does this pipeline actually need at runtime?
#
# The dependency list in pixi.toml was written before the code was finished, and most of it is not
# ours: sceptre reaches its own dependencies through `::` calls, so what we must ship is
# "everything sceptre touches on the paths we call", which is not the same as its DESCRIPTION
# Imports and not the same as its NAMESPACE either.
#
# Three distinct things get confused here, so this reports them separately:
#
#   NAMESPACE imports   loaded whenever sceptre is loaded, unconditionally. Non-negotiable.
#   `::` references     needed only if that line of code executes. Found by scanning the
#                       package's own function bodies -- present here, but possibly on a path
#                       this pipeline never takes.
#   LinkingTo           C++ headers used at compile time only. Needed to build sceptre, never to
#                       run it, so it belongs in the build environment.
#
# What this cannot tell you is which `::` references lie on *our* path specifically -- that would
# need coverage instrumentation. It narrows the question to a list short enough to reason about.
#
# Usage:  audit_dependencies.R [--pkg sceptre]

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(".")
  }
  for (dir in unique(c(file.path(dirname(here), "lib"), file.path(here, "lib"), here))) {
    if (file.exists(file.path(dir, "cli.R"))) {
      source(file.path(dir, "cli.R"))
      return(invisible(NULL))
    }
  }
  stop("Cannot find lib/cli.R relative to ", here, call. = FALSE)
})

option_list <- list(
  make_option("--pkg", type = "character", default = "sceptre", dest = "pkg",
              help = "Package to audit [default %default].")
)
opts <- parse_args(OptionParser(option_list = option_list))

pkg <- opts$pkg
suppressPackageStartupMessages(library(pkg, character.only = TRUE))

desc <- utils::packageDescription(pkg)
split_field <- function(field) {
  if (is.null(field) || is.na(field)) return(character(0))
  parts <- trimws(strsplit(field, ",")[[1]])
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  setdiff(parts, c("R", ""))
}

imports_field  <- split_field(desc$Imports)
linkingto      <- split_field(desc$LinkingTo)
depends_field  <- split_field(desc$Depends)

## NAMESPACE imports -- loaded unconditionally ======================================================

ns_imports <- sort(unique(names(getNamespaceImports(pkg))))
ns_imports <- setdiff(ns_imports, c("base", pkg))

## `::` references anywhere in the package's own code ==============================================
#
# Deparse every object in the namespace and look for `pkg::`. Deparsing rather than parsing source
# because an installed package keeps no .R files -- only the lazy-load database.

ns <- asNamespace(pkg)
code_text <- unlist(lapply(ls(ns, all.names = TRUE), function(nm) {
  obj <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
  if (is.function(obj)) paste(deparse(obj), collapse = "\n") else NULL
}))

hits <- regmatches(code_text, gregexpr("[A-Za-z][A-Za-z0-9._]*(?=:::?)", code_text, perl = TRUE))
qualified <- sort(unique(unlist(hits)))
qualified <- setdiff(qualified, c("base", pkg))

## REPORT ==========================================================================================

say <- function(...) cat(..., "\n", sep = "")
show_set <- function(label, x) {
  say("  ", label, " (", length(x), ")")
  if (length(x)) say("    ", paste(x, collapse = ", "))
}

say("Dependency audit: ", pkg, " ", desc$Version)
say("")
say("DECLARED")
show_set("Depends  ", depends_field)
show_set("Imports  ", imports_field)
show_set("LinkingTo", linkingto)
say("")
say("ACTUALLY REACHED")
show_set("NAMESPACE imports -- loaded on every library() call", ns_imports)
show_set("referenced by ::  -- loaded only if that code runs", qualified)
say("")

runtime_needed <- sort(unique(c(ns_imports, qualified)))
build_only <- setdiff(union(imports_field, linkingto), runtime_needed)

say("VERDICT")
show_set("runtime", runtime_needed)
show_set("declared but never referenced in code -- build/install only", build_only)
say("")
say("LinkingTo is compile-time by construction: ", paste(linkingto, collapse = ", "))
