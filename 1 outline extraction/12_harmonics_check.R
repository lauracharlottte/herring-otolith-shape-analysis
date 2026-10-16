# ============================================================
# 12-harmonic reconstruction check
# For every sample, reconstructs the outline using only 12
# Elliptic Fourier harmonics and superimposes both the raw
# digitized outline (green) and the 12-harmonic reconstruction
# (red) on top of the original photo, to visually judge
# how much detail is lost/kept at 12 harmonics.
# ============================================================

# -----------------------------
# Libraries
# -----------------------------
library(Momocs)
library(magick)
library(foreach)
library(doParallel)

# -----------------------------
# Settings
# -----------------------------
NB_HARMONICS <- 12
cex <- 1   # point size for the raw outline overlay, as in ex_oto.R

## ---- load project paths ------------------------------------
## Searches upwards from the working directory for config.R (also in a
## subfolder Master-Thesis/) and sources it. Afterwards SHAPE_DIR,
## DIR_IMAGES, DIR_FOURIER, ensure_dir() etc. are available.
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- c(file.path(d, "config.R"), file.path(d, "Master-Thesis", "config.R"))
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) {
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    }
    d <- dirname(d)
  }
})

fix.image.path <- DIR_IMAGES

# Dedicated folder for this check (DIR_CHECK in config.R),
# separate from the "Superimposed/" control images of ex_oto.R.
write.check.path <- ensure_dir(DIR_CHECK)

# -----------------------------
# Load raw outlines
# -----------------------------
load(F_OTO)                                    # loads 'oto'

jpgnames <- paste(names(oto), ".jpg", sep = "")

# ------------------------------------------------------------
# Optional: quick numeric sanity check -- how much harmonic
# power is captured at 12 harmonics, across the whole dataset.
# This needs an Out object built the same way your Fourier
# script does; skip if you don't have 'data'/fac info handy.
# ------------------------------------------------------------
# OTO <- Out(oto)
# harmpow <- calibrate_harmonicpower_efourier(OTO, nb.h = 30)
# cat("Cumulative harmonic power at 12 harmonics (mean across samples):\n")
# print(harmpow$q[NB_HARMONICS])

# -----------------------------
# Reconstruct + superimpose, one image at a time
# -----------------------------
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

results <- foreach(i = seq_along(jpgnames), .errorhandling = "pass",
                   .packages = c("Momocs", "magick")) %dopar% {
                     
                     otoname <- names(oto)[i]
                     jpgname <- jpgnames[i]
                     img_path <- paste(fix.image.path, jpgname, sep = "/")
                     out_path <- file.path(write.check.path, paste0(otoname, ".jpg"))
                     
                     tryCatch({
                       
                       # --- Raw digitized outline (same coordinates used in your
                       #     original outline-extraction superimposition step) ---
                       coo_raw <- oto[[otoname]]
                       
                       # --- 12-harmonic Fourier reconstruction of THIS shape ---
                       # efourier() on a single coo matrix + efourier_i() reconstructs
                       # back into the same (absolute, un-normalized) coordinate space
                       # as the input -- this is what lets us overlay it directly on
                       # the raw outline and the photo, the same way.
                       ef  <- efourier(coo_raw, nb.h = NB_HARMONICS)
                       coo_rec <- efourier_i(ef)
                       
                       # --- Load original photo ---
                       img <- as.raster(magick::image_read(img_path))
                       
                       # --- Scale both outlines into the photo's coordinate space,
                       #     exactly like your original superimposition loop ---
                       x_raw <- coo_raw[, 1] / dim(img)[2] + 1
                       y_raw <- coo_raw[, 2] / dim(img)[1] + 1
                       
                       x_rec <- coo_rec[, 1] / dim(img)[2] + 1
                       y_rec <- coo_rec[, 2] / dim(img)[1] + 1
                       
                       jpeg(filename = out_path, width = dim(img)[2], height = dim(img)[1], units = "px")
                       plot(1:2, type = "n")
                       rasterImage(img, 1, 1, 2, 2)
                       points(x_raw, y_raw, col = "green3", cex = cex, pch = 16)
                       points(x_rec, y_rec, col = "red", cex = cex, pch = 16)
                       legend("topright", legend = c("raw outline", paste0(NB_HARMONICS, "-harmonic reconstruction")),
                              col = c("green3", "red"), pch = 16, bty = "n")
                       dev.off()
                       
                       "ok"
                       
                     }, error = function(e){
                       cat("Error on", otoname, ":", conditionMessage(e), "\n")
                       paste("ERROR:", conditionMessage(e))
                     })
                   }

stopCluster(cl)

n_ok    <- sum(unlist(results) == "ok", na.rm = TRUE)
n_error <- length(results) - n_ok

cat("\n12-harmonic reconstruction check finished.\n")
cat("Successful:", n_ok, "of", length(jpgnames), "\n")
if(n_error > 0){
  cat("Errors:", n_error, "-- check console output above for names.\n")
}
cat("Superimposed comparison images saved to:\n")
cat(normalizePath(write.check.path, mustWork = FALSE), "\n")
