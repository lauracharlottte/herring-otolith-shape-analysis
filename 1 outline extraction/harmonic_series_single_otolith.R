## ============================================================
## Harmonic series for ONE otolith
## ------------------------------------------------------------
## Reconstructs the same otolith outline with 1, 2, ... 12
## harmonics and superimposes each reconstruction on the original
## photo -- one image per harmonic, all in a dedicated folder.
## Shows how much shape detail every added harmonic buys.
##
## Used for Fig. 13 (see 4 figures and tables/fig_harmonic_series.R).
##
## NOTE on the outline: oto.RData only holds the *_Right.jpg
## extractions from ex_oto.R. For any other photo (a _Left one,
## for instance) this script does the binarization and outline
## extraction itself, with the exact parameters ex_oto.R uses.
## ============================================================

library(Momocs)
library(imager)
library(magick)

## ------------------------------------------------------------
## 0. Settings
## ------------------------------------------------------------

## The photo, as it is named in Outline Extraction/Images/.
## Can also be passed on the command line, which is handy for a
## second specimen without editing the file:
##   Rscript "1 outline extraction/harmonic_series_single_otolith.R" <photo.jpg>
IMAGE_NAME <- "Geras_2023_SD23__0065_W081_sharpest1_Left.jpg"

local({
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1 && nzchar(a[1])) IMAGE_NAME <<- a[1]
})

## Harmonics to render -- one output image each
HARMONICS <- 1:12

## Colours. Colours, line widths and the output folder can be set
## before source()-ing this file (4 figures and tables/fig_harmonic_series.R
## does that); otherwise these defaults apply.
if (!exists("COL_RAW")) COL_RAW <- "#2271B5"   # raw digitized outline
if (!exists("COL_REC")) COL_REC <- "#D9702D"   # harmonic reconstruction

## How the reconstruction is drawn:
##   "line"   -> continuous curve; thick enough to survive being
##               shrunk to a few cm in the thesis  <- use this one
##   "points" -> dotted, like the 12-harmonic check images
REC_STYLE  <- "line"
REC_NB_PTS <- 720      # points on the reconstructed outline; high
                       # so the "line" style comes out smooth
CEX_RAW    <- 1.0      # point size, raw outline
CEX_REC    <- 1.4      # point size, reconstruction ("points" style)

## Line widths. These are multiplied by the panel resolution, so a
## panel rendered at twice the size gets twice the line width and
## the figure looks the same at any output size.
## The reconstruction is the thick line, the raw outline a thin one
## underneath it. That works here because blue and orange are far
## enough apart to stay separable: where the fit is good the orange
## simply covers the blue, where it is off you see a thin blue
## contour next to a fat orange curve -- which is the whole point
## of the plate.
if (!exists("LWD_REC")) LWD_REC <- 9.0   # reconstruction
if (!exists("LWD_RAW")) LWD_RAW <- 5.0   # raw outline ("line" style only)

## ------------------------------------------------------------
## Figure layout -- this is what makes the series usable as a
## thesis plate rather than a folder of screenshots
## ------------------------------------------------------------

## Crop every panel to the outline's bounding box plus this much
## padding (fraction of the otolith's size). The photo is mostly
## black background; cropping it away makes the otolith itself
## 2-3x larger on the page at the same figure width.
## Set to NA to keep the full photo.
CROP_PAD <- 0.05

## Panel resolution as a multiple of the photo's own pixels.
## 1 = native (~1900 x 3000 px here), which is far above 300 dpi
## for anything you would print. Lower it to 0.5 for smaller files.
PANEL_SCALE <- 1

## "png" is lossless -- better for thin lines than jpeg, and what
## you want to hand to LaTeX/Word. "jpg" for smaller files.
PANEL_FORMAT <- "png"

## Big harmonic number burned into each panel, so the plate is
## readable without a legend under every tile.
LABEL       <- TRUE
LABEL_COL   <- "white"
LABEL_CEX   <- 6.5     # scaled by panel resolution, like the lines;
                       # 6.5 = about 9 pt when the 4-column plate is
                       # printed 11-12 cm wide (2.6 came out at ~3 pt)

## Add the cumulative harmonic power to the label, e.g.
## "h = 5   97.4 %". This is the number the choice of 12 harmonics
## rests on (99 % of shape information), so the plate then shows
## the criterion and the visual result side by side.
LABEL_POWER <- FALSE

## Harmonics the power is normalised over. 30 is what
## calibrate_harmonicpower_efourier() uses in the main pipeline --
## keep it there, or the percentages stop being comparable with the
## ones in Fourier_analysis_ex1.R. (Normalising over 1:12 instead
## would make h = 12 come out as 100 % by construction.)
POWER_NB_H <- 30

## Colour key: "first" = only on panel 1 (put the rest in the
## figure caption), "all", or "none".
LEGEND     <- "none"   # the key is in the thesis caption; a legend
                       # small enough to fit panel 1 is unreadable in print
LEGEND_CEX <- 1.9      # scaled by panel resolution, like the lines

## The contact sheet
MAKE_MONTAGE   <- TRUE
MONTAGE_COLS   <- 4      # 4 -> 4x3 plate; 6 -> 6x2, a wider, shorter figure
## px per tile in the montage. NA = the panels' own resolution.
## This photo is 703 x 1297 px, so anything above that is
## interpolation, not detail -- 4 native tiles across a 16 cm
## figure are already ~380 dpi, well past what print resolves.
MONTAGE_TILE_W <- NA
MONTAGE_PAD    <- 10     # px of white between tiles

## Binarization parameters -- identical to ex_oto.R, do not change
## unless you change them there too.
xrange    <- 50
yrange    <- 50
adjust    <- 0
threshold <- 0.5

## ---- load project paths ------------------------------------
## Searches upwards from the working directory for config.R
## and sources it. Afterwards SHAPE_DIR,
## DIR_IMAGES, DIR_FOURIER, ensure_dir() etc. are available.
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- file.path(d, "config.R")
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) {
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    }
    d <- dirname(d)
  }
})

img.dir  <- DIR_IMAGES
oto_name <- sub("\\.jpg$", "", IMAGE_NAME)

## Folder at the project root, one subfolder
## per specimen so a second otolith does not overwrite the first.
## SERIES_OUT_DIR, if set before source(), replaces Harmonic_Series/.
out.dir <- ensure_dir(file.path(
  if (exists("SERIES_OUT_DIR")) SERIES_OUT_DIR else file.path(SHAPE_DIR, "Harmonic_Series"),
  oto_name))

img_path <- file.path(img.dir, IMAGE_NAME)
if (!file.exists(img_path)) {
  stop("Photo not found:\n  ", img_path, call. = FALSE)
}

cat("Photo      :", IMAGE_NAME, "\n")
cat("Harmonics  :", paste(range(HARMONICS), collapse = " to "), "\n")
cat("Output     :", out.dir, "\n\n")

## ------------------------------------------------------------
## 1. Get the raw outline
## ------------------------------------------------------------
## Fast path: the specimen is already in oto.RData (all _Right
## photos). Otherwise binarize + extract it here, exactly the way
## ex_oto.R does, and keep the binary image next to the others.
coo_raw <- NULL

if (file.exists(F_OTO)) {
  load(F_OTO)                                   # loads 'oto'
  if (oto_name %in% names(oto)) {
    coo_raw <- oto[[oto_name]]
    cat("Outline taken from oto.RData\n")
  }
}

if (is.null(coo_raw)) {

  cat("Not in oto.RData -- extracting the outline from the photo\n")

  binary.dir  <- ensure_dir(file.path(img.dir, "Binary"))
  binary_path <- file.path(binary.dir, IMAGE_NAME)

  if (!file.exists(binary_path)) {
    a <- load.image(img_path)
    im_width  <- width(a)
    im_height <- height(a)

    ## Grayscale, white rectangle over the centre so dark spots on
    ## the otolith surface cannot break the contour, invert, threshold.
    a <- a %>% grayscale() %>%
      draw_rect(im_width / 2 - xrange, im_height / 2 - yrange,
                im_width / 2 + xrange, im_height / 2 + yrange,
                color = 1, opacity = 1, filled = TRUE) %>%
      imchange(where = ~ xs > 0, fo = ~ 1 - .) %>%
      threshold("auto", approx = TRUE, adjust = 1.4 + adjust)

    a %>% as.cimg() %>% imager::save.image(file = binary_path)
    cat("  binary image written to Images/Binary/\n")
  } else {
    cat("  binary image already there, reusing it\n")
  }

  imported <- import_jpg(binary_path, threshold = threshold)
  coo_raw  <- imported[[1]]
}

coo_raw <- as.matrix(coo_raw)
cat("Outline points:", nrow(coo_raw), "\n")
if (nrow(coo_raw) < 100) {
  stop("Only ", nrow(coo_raw), " outline points -- the extraction went wrong. ",
       "Check Images/Binary/", IMAGE_NAME, call. = FALSE)
}
if (max(HARMONICS) > nrow(coo_raw) / 2 - 1) {
  stop("Too few outline points for ", max(HARMONICS), " harmonics.", call. = FALSE)
}

## ------------------------------------------------------------
## 2. One EFD fit, reconstructed at each harmonic count
## ------------------------------------------------------------
## The coefficients of harmonic k do not depend on how many
## harmonics were fitted, so a single fit at max(HARMONICS) can be
## truncated to 1, 2, ... n. Same approach Momocs itself uses in
## calibrate_reconstructions_efourier().
##
## No normalisation here: efourier_i() then returns the outline in
## the photo's own coordinate space, which is what lets it be laid
## over the picture directly.
ef <- efourier(coo_raw, nb.h = max(HARMONICS))

## Cumulative harmonic power for this otolith, the same way
## calibrate_harmonicpower_efourier() computes it: fit far more
## harmonics than are drawn, drop the first one (the fitting
## ellipse carries almost all the power on its own and would
## flatten the whole curve), then accumulate.
## pow_cum[h] is the share reached at h harmonics; h = 1 is the
## dropped baseline and has no value.
pow_cum <- rep(NA_real_, max(HARMONICS))
if (LABEL_POWER) {
  nb_pow <- min(POWER_NB_H, floor(nrow(coo_raw) / 2) - 1)
  if (nb_pow > max(HARMONICS)) {
    hp  <- harm_pow(efourier(coo_raw, nb.h = nb_pow))[-1]   # drop harmonic 1
    cum <- cumsum(hp) / sum(hp) * 100
    pow_cum[2:max(HARMONICS)] <- cum[seq_len(max(HARMONICS) - 1)]
  } else {
    message("Not enough outline points for a ", POWER_NB_H,
            "-harmonic power reference -- labels stay plain.")
  }
}

img   <- as.raster(magick::image_read(img_path))
im_h  <- dim(img)[1]
im_w  <- dim(img)[2]

## Both outlines are mapped into the [1, 2] x [1, 2] box that the
## photo is drawn into -- same mapping as ex_oto.R.
to_plot_x <- function(x) x / im_w + 1
to_plot_y <- function(y) y / im_h + 1

x_raw <- to_plot_x(coo_raw[, 1])
y_raw <- to_plot_y(coo_raw[, 2])

## ------------------------------------------------------------
## 2b. Crop window
## ------------------------------------------------------------
## Cropping happens through the plot window, not by cutting the
## photo: xlim/ylim select a sub-box of the [1,2] square while the
## photo is still drawn across the full square and clipped to the
## plot region. One coordinate system, no offset arithmetic.
##
## The box comes from the RAW outline, so it is identical for every
## harmonic -- all panels line up and the otolith does not jump
## around in the montage.
if (is.na(CROP_PAD)) {
  xlim <- c(1, 2)
  ylim <- c(1, 2)
} else {
  xr <- range(x_raw); yr <- range(y_raw)
  xlim <- xr + c(-1, 1) * CROP_PAD * diff(xr)
  ylim <- yr + c(-1, 1) * CROP_PAD * diff(yr)
  ## never ask for more than the photo actually has
  xlim <- c(max(1, xlim[1]), min(2, xlim[2]))
  ylim <- c(max(1, ylim[1]), min(2, ylim[2]))
}

## Panel size in pixels. diff(xlim) * im_w is exactly the width of
## the crop in original photo pixels, so the aspect ratio is kept
## and nothing is stretched (mar = 0 + xaxs/yaxs = "i" make the
## plot region equal to the device).
panel_w <- max(1, round(diff(xlim) * im_w * PANEL_SCALE))
panel_h <- max(1, round(diff(ylim) * im_h * PANEL_SCALE))

## Everything drawn on top scales with the panel, so the figure
## looks identical whatever PANEL_SCALE is set to.
k       <- panel_w / 700
lwd_rec <- LWD_REC * k
lwd_raw <- LWD_RAW * k
cex_txt <- max(1, k)

cat("Panel size :", panel_w, "x", panel_h, "px",
    if (!is.na(CROP_PAD)) sprintf("(cropped to %.0f%% of the photo area)",
                                  100 * diff(xlim) * diff(ylim)) else "", "\n\n")

open_panel <- function(path) {
  if (PANEL_FORMAT == "png") {
    png(filename = path, width = panel_w, height = panel_h, units = "px")
  } else {
    jpeg(filename = path, width = panel_w, height = panel_h, units = "px",
         quality = 95)
  }
}

files <- character(0)

for (h in HARMONICS) {

  coo_rec <- efourier_i(ef, nb.h = h, nb.pts = REC_NB_PTS)
  x_rec   <- to_plot_x(coo_rec[, 1])
  y_rec   <- to_plot_y(coo_rec[, 2])

  out_path <- file.path(out.dir,
                        sprintf("%s_h%02d.%s", oto_name, h, PANEL_FORMAT))

  open_panel(out_path)
  op <- par(mar = c(0, 0, 0, 0),          # photo fills the frame, no axes
             family = FONT_PLOT)          # Times New Roman, like every thesis figure
  plot(xlim, ylim, type = "n", axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")
  rasterImage(img, 1, 1, 2, 2)            # clipped to the plot region

  ## Raw outline first, reconstruction on top: where the two agree
  ## you see the reconstruction, where they part you see both.
  if (REC_STYLE == "line") {
    lines(c(x_raw, x_raw[1]), c(y_raw, y_raw[1]), col = COL_RAW, lwd = lwd_raw)
  } else {
    points(x_raw, y_raw, col = COL_RAW, cex = CEX_RAW, pch = 16)
  }

  if (REC_STYLE == "line") {
    ## close the loop so no gap shows at the starting point
    lines(c(x_rec, x_rec[1]), c(y_rec, y_rec[1]), col = COL_REC, lwd = lwd_rec)
  } else {
    points(x_rec, y_rec, col = COL_REC, cex = CEX_REC * k, pch = 16)
  }

  if (LABEL) {
    ## italic h, as in the caption ("h = 1 to 12")
    lab <- bquote(bolditalic("h") * bold(" = " * .(h)))
    if (LABEL_POWER && !is.na(pow_cum[h])) {
      lab <- bquote(bolditalic("h") * bold(" = " * .(h) * "     " *
                                            .(sprintf("%.1f%%", pow_cum[h]))))
    }
    text(xlim[1] + 0.03 * diff(xlim), ylim[2] - 0.02 * diff(ylim),
         labels = lab, adj = c(0, 1),
         col = LABEL_COL, cex = LABEL_CEX * k, font = 2)
  }

  show_legend <- LEGEND == "all" || (LEGEND == "first" && h == HARMONICS[1])
  if (show_legend) {
    legend("bottomright",
           legend = c("raw outline", "reconstruction"),
           col    = c(COL_RAW, COL_REC),
           lwd    = if (REC_STYLE == "line") c(lwd_raw, lwd_rec) else c(NA, NA),
           pch    = if (REC_STYLE == "line") c(NA, NA) else c(16, 16),
           bty = "n", text.col = "white",
           cex = LEGEND_CEX * k, pt.cex = LEGEND_CEX * k)
  }

  par(op)
  dev.off()

  files <- c(files, out_path)
  cat("  h =", formatC(h, width = 2), "->", basename(out_path), "\n")
}

## ------------------------------------------------------------
## 3. Contact sheet of the whole series
## ------------------------------------------------------------
if (MAKE_MONTAGE && length(files) > 1) {

  montage_path <- file.path(out.dir, paste0(oto_name, "_series.png"))
  n_col <- min(MONTAGE_COLS, length(files))
  n_row <- ceiling(length(files) / n_col)

  ## geometry here is padding only -- the tiles are already scaled,
  ## and magick rejects a combined size+offset string.
  tiles <- magick::image_read(files)
  if (!is.na(MONTAGE_TILE_W)) {
    tiles <- magick::image_scale(tiles, paste0(MONTAGE_TILE_W, "x"))
  }

  tiles %>%
    magick::image_montage(tile = paste0(n_col, "x"),
                          geometry = sprintf("+%d+%d", MONTAGE_PAD, MONTAGE_PAD),
                          bg = "white") %>%
    magick::image_write(montage_path, format = "png")

  info <- magick::image_info(magick::image_read(montage_path))
  cat("\nContact sheet:", basename(montage_path), "--",
      info$width, "x", info$height, "px\n")

  ## What that plate actually measures on paper. A thesis text
  ## column is about 16 cm; the question is only whether the height
  ## still fits the page once the caption is there.
  for (w_cm in c(12, 16)) {
    h_cm <- w_cm * info$height / info$width
    cat(sprintf("  at %2.0f cm wide -> %.1f cm tall, %.0f dpi\n",
                w_cm, h_cm, info$width / (w_cm / 2.54)))
  }
  cat("  (", n_col, " x ", n_row, " tiles; set MONTAGE_COLS <- 6 for a ",
      "wider, shorter plate)\n", sep = "")

  ## Caption skeleton -- the colour key lives here rather than in
  ## twelve little legends inside the plate.
  cat("\nCaption to start from:\n")
  cat("  Elliptic Fourier reconstruction of one herring otolith",
      "(", oto_name, ") with 1 to", max(HARMONICS), "harmonics.\n")
  cat("  Blue: outline digitized from the photograph. Orange:",
      "reconstruction from the given\n  number of harmonics.")
  if (LABEL_POWER && any(!is.na(pow_cum))) {
    cat(" Percentages give the cumulative harmonic power\n",
        " (normalised over harmonics 2 to ", nb_pow, ").", sep = "")
  }
  cat("\n")
}

cat("\nDone.", length(files), "images in:\n ", out.dir, "\n")
