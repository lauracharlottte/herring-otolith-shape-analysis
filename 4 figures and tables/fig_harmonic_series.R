# ============================================================
# Figure 13 -- reconstruction of one otolith outline with an
# increasing number of harmonics (h = 1 to 12)
# ------------------------------------------------------------
# Drawn by 1 outline extraction/harmonic_series_single_otolith.R
# with the thesis settings:
#   - raw outline in the strong blue of the palette (NSAS)
#   - thick lines, so they survive being printed small
#
# Output: results final plots/results chapter/
#     png/fig4_harmonic_series.png          the plate
#     harmonic series tiles/<otolith>/      the 12 single panels
# ============================================================

## ---- load project paths and the shared palette -------------
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- c(file.path(d, "config.R"), file.path(d, "Master-Thesis", "config.R"))
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d))
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    d <- dirname(d)
  }
})

fig.dir <- ensure_dir(file.path(DIR_FINAL, "results chapter"))
rev.png <- ensure_dir(file.path(fig.dir, "png"))

COL_RAW <- COL_STOCK[["NSAS"]]   # strong blue
COL_REC <- COL_STOCK[["WBSS"]]   # orange, as before
LWD_RAW <- 9                     # was 5
LWD_REC <- 14                    # was 9
SERIES_OUT_DIR <- file.path(fig.dir, "harmonic series tiles")

## the otolith used for Figure 4 (SERIES_ID); the
## series script reads its photo name from IMAGE_NAME via the command
## line, so it is passed the same way here
SERIES_ID <- "Geras_2023_SD23__0075_W091_sharpest1_Right"
commandArgs <- function(trailingOnly = FALSE) paste0(SERIES_ID, ".jpg")

source(local = TRUE, file.path(REPO_DIR, "1 outline extraction", "harmonic_series_single_otolith.R"))

rm(commandArgs)
file.copy(montage_path, file.path(rev.png, "fig4_harmonic_series.png"), overwrite = TRUE)
cat("Plate copied to:", file.path(rev.png, "fig4_harmonic_series.png"), "\n")
