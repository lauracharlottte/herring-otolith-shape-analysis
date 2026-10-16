# ============================================================
# Appendix -- examples of otoliths removed as outliers
# ------------------------------------------------------------
# Six control images (extracted outline drawn on the photograph,
# from ex_oto.R) of otoliths removed by the D^2 rule, one or two per
# cause found in the visual check of the control images of all 129
# removed otoliths. The R plot frame and axes of the
# control images are cropped away; every panel is scaled to the same
# height.
# Output: Shape Analysis/Appendix/{png,pdf}/figA_outlier_examples
# ============================================================
suppressPackageStartupMessages(library(magick))
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- c(file.path(d, "config.R"), file.path(d, "Master-Thesis", "config.R"))
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) stop("config.R not found -- please setwd() into the repository folder.", call. = FALSE)
    d <- dirname(d)
  }
})
out <- proj_path("Appendix")
invisible(ensure_dir(file.path(out, "png"))); invisible(ensure_dir(file.path(out, "pdf")))

## panel order a-f: (a, b) outline continues onto the neighbouring otolith,
## (c, d) tissue or fibre attached, (e) only part of the outline traced,
## (f) crystalline otolith with debris
ids <- c("Geras_2025_SD23__0086_W016_sharpest2_Right", "Geras_2023_SD21__0135_W044_sharpest1_Right",
         "Geras_2025_SD22__0033_W040_sharpest2_Right", "Geras_2023_SD21__0146_W057_sharpest1_Right",
         "Geras_2023_SD23__0070_W086_sharpest1_Right", "Geras_2023_SD22__0189_W052_sharpest2_Right")

crop_plot <- function(im) {
  ## the photograph is the dark rectangle inside the white plot frame:
  ## keep the rows and columns that are mostly dark
  g <- as.integer(image_data(image_convert(im, colorspace = "gray"), "gray"))[, , 1]
  dark <- g < 90
  rows <- which(rowMeans(dark) > 0.25); cols <- which(colMeans(dark) > 0.25)
  y0 <- min(rows); y1 <- max(rows); x0 <- min(cols); x1 <- max(cols)
  image_crop(im, sprintf("%dx%d+%d+%d", x1 - x0 - 5, y1 - y0 - 5, x0 + 3, y0 + 3))
}

H <- 900
tiles <- lapply(seq_along(ids), function(i) {
  im <- image_read(file.path(DIR_IMAGES, "Superimposed", paste0(ids[i], ".jpg")))
  im <- image_scale(crop_plot(im), paste0("x", H))
  im <- image_annotate(im, letters[i], size = 70, color = "white", weight = 700,
                       gravity = "northwest", location = "+22+12", font = "Times New Roman")
  image_border(im, "white", "14x14")
})
row1 <- image_append(do.call(c, tiles[1:3]))
row2 <- image_append(do.call(c, tiles[4:6]))
w <- max(image_info(row1)$width, image_info(row2)$width)
row1 <- image_extent(row1, paste0(w, "x", image_info(row1)$height), color = "white")
row2 <- image_extent(row2, paste0(w, "x", image_info(row2)$height), color = "white")
fig <- image_append(c(row1, row2), stack = TRUE)
image_write(fig, file.path(out, "png", "figA_outlier_examples.png"), format = "png", density = 300)
image_write(fig, file.path(out, "pdf", "figA_outlier_examples.pdf"), format = "pdf")
cat("figA_outlier_examples:", image_info(fig)$width, "x", image_info(fig)$height, "px\n")
