## ============================================================
## config.R -- central path configuration
## ------------------------------------------------------------
## The only place in the project where paths are defined.
## No script contains an absolute path; all of them derive their
## inputs and outputs from the variables defined here.
##
## Usage (already at the top of every script):
##     local({ ... })   # finds this file and sources it
##
## The project root ("Shape Analysis") is determined in this order:
##   1. environment variable SHAPE_ANALYSIS_DIR, if set
##      (e.g. in ~/.Renviron:  SHAPE_ANALYSIS_DIR=/path/to/Shape Analysis)
##   2. upwards from the location of this file, until a folder
##      contains the marker folder "Outline Extraction"
##   3. the parent folder of this file
## ============================================================

## ------------------------------------------------------------
## Folder by which the project root is recognised
## ------------------------------------------------------------
SHAPE_MARKER <- "Outline Extraction"

## ------------------------------------------------------------
## Locate this file (config.R)
## ------------------------------------------------------------
## Walks through the call frames from the inside out: the innermost frame
## with $ofile belongs to config.R itself. Falls back to getwd() if the
## file was not loaded with source().
.shape_config_dir <- function() {
  n <- sys.nframe()
  if (n > 0) {
    for (i in rev(seq_len(n))) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
}

## ------------------------------------------------------------
## Determine the project root
## ------------------------------------------------------------
SHAPE_DIR <- local({

  ## -- 1. explicit setting via environment variable --
  env <- Sys.getenv("SHAPE_ANALYSIS_DIR", unset = "")
  if (nzchar(env)) {
    env <- normalizePath(path.expand(env), mustWork = FALSE)
    if (!dir.exists(env)) {
      stop("SHAPE_ANALYSIS_DIR points to a folder that does not exist:\n  ",
           env, call. = FALSE)
    }
    return(env)
  }

  ## -- locate this file --
  ## sys.frames() knows the path when config.R is loaded with source().
  ## Search from the inside out: the innermost frame with $ofile is
  ## config.R itself, outer frames would be the calling script.
  here <- .shape_config_dir()

  ## -- 2. walk upwards until the marker folder appears --
  ## First from the location of this file, then from the working directory.
  for (start in unique(c(here, normalizePath(getwd(), mustWork = FALSE)))) {
    d <- start
    repeat {
      if (dir.exists(file.path(d, SHAPE_MARKER))) return(d)
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }

  ## -- 3. fallback: one level above this file --
  dirname(here)
})

## ------------------------------------------------------------
## Repository folder (this Git repository)
## ------------------------------------------------------------
REPO_DIR <- .shape_config_dir()

## ------------------------------------------------------------
## Derived folders
## ------------------------------------------------------------
DIR_OUTLINE  <- file.path(SHAPE_DIR, "Outline Extraction")
DIR_IMAGES   <- file.path(DIR_OUTLINE, "Images")
DIR_FOURIER  <- file.path(SHAPE_DIR, "Fourier and Classification")
DIR_RESULTS  <- file.path(DIR_FOURIER, "Results")   # all analysis results
DIR_CHECK    <- file.path(SHAPE_DIR, "Harmonic12_Check")
DIR_DATASETS <- local({
  for (nm in c("data sets", "data_sets")) {
    p <- file.path(SHAPE_DIR, nm)
    if (dir.exists(p)) return(p)
  }
  file.path(SHAPE_DIR, "data sets")
})

## ------------------------------------------------------------
## Frequently used single files
## ------------------------------------------------------------
F_OTO        <- file.path(DIR_IMAGES, "oto.RData")
## Metadata with age, written by 2 metadata/age_join.R
F_DATA_AGE   <- file.path(DIR_DATASETS, "Otolith_Shape_DATA_SET_with_Age.csv")
F_SWEEP_CSV  <- file.path(DIR_RESULTS, "Harmonics_sweep",
                          "harmonics_accuracy_kappa.csv")

## ------------------------------------------------------------
## Folder for the figures that go into the thesis
## ------------------------------------------------------------
DIR_FINAL <- file.path(SHAPE_DIR, "results final plots")

## ------------------------------------------------------------
## Colour scheme -- a single source for all figures
## ------------------------------------------------------------
## Every population keeps its colour in every figure, so that the
## reader can carry a colour from one panel to the next.
##
## Checked as a scatter-plot palette, i.e. over ALL pairs, not only
## neighbouring ones (the PCA plot needs this):
##     worst pair, normal vision        20.3   (threshold 15)
##     worst pair, protan/deuteran      10.1   (threshold  8)
##     worst pair, tritan               10.2   (threshold  8)
## Distances are OKLab delta E x100, colour-vision-deficiency simulation
## after Machado, Oliveira & Fernandes (2009), severity 1.0.
##
## CAUTION when changing: a strong mid-tone green merges with the
## orange of WBSS under red-green colour blindness (delta E down to 0.1
## measured). Green must differ from the orange (L = 0.66) in lightness
## or move towards blue-green.
COL_STOCK <- c(
  BAS   = "#07306B",   # dark navy blue
  CBH   = "#5ACE69",   # fresh green
  Downs = "#169597",   # Petrol / Teal
  NSAS  = "#2B5BD1",   # royal blue
  WBSS  = "#D9702D"    # Orange
)

## Order of the populations in ALL figures and tables (legends,
## stacks, axes, panels, columns): WBSS first, as the reference, then
## BAS, NSAS, Downs, CBH. The palette above is checked over ALL pairs,
## so every neighbourhood in stacks and legends stays distinguishable.
## In the figures the groups are called "populations", not "stocks";
## the variable names in the code (Assigned.Stock, COL_STOCK ...) stay.
POP_ORDER   <- c("WBSS", "BAS", "NSAS", "Downs", "CBH")
STOCK_ORDER <- POP_ORDER   # alias used in some figure scripts

## NEATL (n = 2) and NSS (n = 3) deliberately get NO colour of their own.
## Seven categorical hues are more than colour alone can carry -- every
## attempt ended at a CVD distance of 1.5-4.3. For five fish a second
## channel is enough: shape instead of colour. (Both are excluded from
## the analysis in main_analysis.R.)
COL_TINY   <- "#5A6572"
SHAPE_TINY <- c(NEATL = 2, NSS = 0)   # open triangle, open square

## Spawning group: deliberately coded by season.
COL_SPAWN <- c(
  "spring spawner" = "#5ACE69",
  "autumn spawner" = "#D9702D",
  "winter spawner" = "#2B5BD1"
)

## ------------------------------------------------------------
## Font for all figures
## ------------------------------------------------------------
## Times New Roman, so that the figures match the body text of the
## thesis. On a system without this font the block falls back to
## "Times", which is metrically identical and available everywhere.
FONT_PLOT <- local({
  cand <- c("Times New Roman", "Times", "serif")
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    have <- unique(systemfonts::system_fonts()$family)
    hit <- cand[cand %in% have]
    if (length(hit)) return(hit[1])
  }
  "serif"
})

COL_BAR <- "#169597"   # single-colour bars (variable importance)
COL_INK <- "#101821"   # text in figures -- never a series colour
COL_MUT <- "#6B7885"   # secondary text, annotations in figures

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

## Build a path relative to the project root: proj_path("Outline Extraction", "Images")
proj_path <- function(...) file.path(SHAPE_DIR, ...)

## Create a folder if needed and return its path,
## so that out.dir can be used directly in file.path().
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

## Check that the expected data folders really exist.
check_paths <- function(stop_on_missing = FALSE) {
  needed <- c(
    "Project root"            = SHAPE_DIR,
    "Outline Extraction"      = DIR_OUTLINE,
    "Images"                  = DIR_IMAGES,
    "Fourier and Classification" = DIR_FOURIER,
    "data sets"               = DIR_DATASETS
  )
  ok <- vapply(needed, function(p) dir.exists(p) || file.exists(p), logical(1))
  for (i in seq_along(needed)) {
    cat(if (ok[i]) "  OK    " else "  MISSING ", names(needed)[i], ": ", needed[i], "\n", sep = "")
  }
  if (stop_on_missing && !all(ok)) {
    stop("Folders are missing -- please set SHAPE_ANALYSIS_DIR.", call. = FALSE)
  }
  invisible(all(ok))
}

cat("Project root:", SHAPE_DIR, "\n")
