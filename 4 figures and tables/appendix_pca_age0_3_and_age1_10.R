# ============================================================
# Fig. A3 -- PCAs of the age subsets 0-3 and 1-10
# ------------------------------------------------------------
# Same design as the two-panel PCA in the Results chapter
# (fig_pca.R): a age 0-3, b age 1-10. Both are
# separate PCAs; their axes are flipped where needed to match the
# PCA of the full data set, as for the age 1-3 panel.
# Output: Shape Analysis/Appendix/{pdf,png}/figA_pca_age0_3_and_age1_10
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(Momocs)
})

## ---- load project paths and the shared palette -------------
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- file.path(d, "config.R")
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d))
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    d <- dirname(d)
  }
})

RESULTS   <- DIR_RESULTS
fig.dir <- ensure_dir(file.path(SHAPE_DIR, "Appendix"))
rev.pdf <- ensure_dir(file.path(fig.dir, "pdf"))
rev.png <- ensure_dir(file.path(fig.dir, "png"))
FAM <- FONT_PLOT

## ---- data: the same objects as Figs. 15 and 16 ----------------
load(file.path(RESULTS, "dataBIO.RData"))   # data
e6 <- local({ e <- new.env(); load(file.path(RESULTS, "eOTO.RData"), envir = e); e$eOTO })
stopifnot(nrow(e6$coe) == nrow(data))
age6 <- suppressWarnings(as.numeric(as.character(e6$fac$Age)))
eA <- Momocs::slice(e6, which(!is.na(age6) & age6 >= 0 & age6 <= 3))
eB <- Momocs::slice(e6, which(!is.na(age6) & age6 >= 1 & age6 <= 10))

prep <- function(eO) {
  eO$fac$Assigned.Stock <- factor(as.character(eO$fac$Assigned.Stock),
    levels = intersect(POP_ORDER, unique(as.character(eO$fac$Assigned.Stock))))
  eO
}
e6 <- prep(e6); eA <- prep(eA); eB <- prep(eB)

pc6 <- PCA(e6)
pcA <- PCA(eA)
pcB <- PCA(eB)
## align the signs of both panels with the full data set (as for age 1-3)
align <- function(pc, nm) {
  for (k in 1:2) if (cor(pc$rotation[, k], pc6$rotation[, k]) < 0) {
    pc$x[, k] <- -pc$x[, k]; pc$rotation[, k] <- -pc$rotation[, k]
    cat("PC", k, " of the ", nm, " PCA flipped to match the full data set\n", sep = "")
  }
  pc
}
pcA <- align(pcA, "age 0-3"); pcB <- align(pcB, "age 1-10")

## ---- axis titles "PC1 (34.5%)" ---------------------------------
layer_axis_titles <- function(p, cex = 1.3, col = "grey25") {
  old <- par(xpd = NA); on.exit(par(old))
  ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  u  <- par("usr"); dx <- diff(u[1:2]) / 100; dy <- diff(u[3:4]) / 100
  xpi <- diff(u[1:2]) / par("pin")[1]; ypi <- diff(u[3:4]) / par("pin")[2]
  lab <- sprintf("PC%d (%.1f%%)", p$axes, ve[p$axes])
  w <- strwidth(lab, units = "inches", cex = cex) * 1.06
  h <- strheight(lab, units = "inches", cex = cex) * 1.6
  ## no white box behind the titles: it cut through the ellipses; the
  ## points are light enough for the titles to stay readable
  text(u[2] - dx, -dy, lab[1], adj = c(1, 1), cex = cex, col = col)
  text(dx, u[4] - dy, lab[2], adj = c(1, 1), srt = 90, cex = cex, col = col)
  invisible(p)
}

## group means as large points with a white rim, drawn on top
layer_group_means <- function(p, cex = 2.6) {
  f <- p$f
  for (i in seq_along(levels(f))) {
    xy <- p$xy[f == levels(f)[i], , drop = FALSE]
    if (!nrow(xy)) next
    m <- colMeans(xy)
    points(m[1], m[2], pch = 21, bg = p$colors_groups[i], col = "white",
           cex = cex, lwd = 1.6)
  }
  invisible(p)
}

## panel letter and subset name, top left on a white background
layer_panel_label <- function(lab, cex = 1.35) {
  old <- par(xpd = NA); on.exit(par(old))
  u <- par("usr"); dx <- diff(u[1:2]) / 100; dy <- diff(u[3:4]) / 100
  w <- strwidth(lab, cex = cex); h <- strheight(lab, cex = cex)
  rect(u[1] + dx * 0.5, u[4] - dy * 1.5 - h * 1.5, u[1] + dx * 2.5 + w, u[4] - dy * 0.5,
       col = grDevices::adjustcolor("white", alpha.f = 0.85), border = NA)
  text(u[1] + dx * 1.5, u[4] - dy * 1.5, lab, adj = c(0, 1), cex = cex, col = COL_INK)
}

## panel letter and subset name, top left on a white background
layer_panel_label <- function(lab, cex = 1.35) {
  old <- par(xpd = NA); on.exit(par(old))
  u <- par("usr"); dx <- diff(u[1:2]) / 100; dy <- diff(u[3:4]) / 100
  w <- strwidth(lab, cex = cex); h <- strheight(lab, cex = cex)
  rect(u[1] + dx * 0.5, u[4] - dy * 1.5 - h * 1.5, u[1] + dx * 2.5 + w, u[4] - dy * 0.5,
       col = grDevices::adjustcolor("white", alpha.f = 0.85), border = NA)
  text(u[1] + dx * 1.5, u[4] - dy * 1.5, lab, adj = c(0, 1), cex = cex, col = COL_INK)
}

CONF_90_2D <- 2 * pnorm(sqrt(qchisq(0.90, df = 2))) - 1   # = 0.968, see panel()

cols <- unname(COL_STOCK[levels(e6$fac$Assigned.Stock)])

## one panel; the panel title is drawn
## above it in a strip of its own (see below)
panel <- function(pc, legend) {
  p <- plot_PCA(pc, ~Assigned.Stock, zoom = 2.2,
                palette = function(n) cols[seq_len(n)],
                chull = FALSE, eigen = FALSE, morphospace = FALSE,
                axesnames = FALSE, axesvar = FALSE, points = FALSE,
                legend = FALSE)
  ## Background outlines on a fixed grid relative to the plot window
  ## (7 columns x 7 rows, kept inside the edges), not on the range of the
  ## scores ("range"): with "range" every PCA got its own grid and the
  ## outlines at the edges were cut off.
  u  <- par("usr")
  fx <- seq(0.07, 0.93, length.out = 7)      # grid as fractions of the window
  fy <- seq(0.07, 0.93, length.out = 7)
  g  <- expand.grid(fx = fx, fy = fy)
  ## cells left free because text sits there: eigenvalue inset (bottom
  ## right), "PC1" (right, on the axis), "PC2" (top, on the axis) and,
  ## only in the panel that has it, the legend (top right)
  near <- function(a, b) abs(a - b) < 1e-9
  skip <- (near(g$fx, fx[7]) & (near(g$fy, fy[1]) | near(g$fy, fy[4]))) |
          (near(g$fx, fx[4]) & (near(g$fy, fy[7]) | near(g$fy, fy[6])))
  if (legend) skip <- skip | (near(g$fx, fx[7]) & (near(g$fy, fy[7]) | near(g$fy, fy[6])))
  g   <- g[!skip, ]
  pos <- cbind(x = u[1] + diff(u[1:2]) * g$fx, y = u[3] + diff(u[3:4]) * g$fy)
  layer_morphospace_PCA(p, position = pos, rotate = 0, size = 0.8, col = "#999999",
                        flipx = FALSE, flipy = FALSE, draw = TRUE)
  layer_points(p, cex = 0.55, transp = 0.45)
  ## Ellipses that really enclose 90 % of the otoliths of a population.
  ## Momocs' conf_ell() uses the one-dimensional normal quantile as radius
  ## (conf = 0.9 -> 1.645 SD, which in two dimensions encloses only ~74 %).
  ## A 2-D 90 % ellipse needs radius sqrt(qchisq(0.9, 2)) = 2.146 SD; this
  ## conf value makes conf_ell() draw exactly that radius.
  layer_ellipses(p, conf = CONF_90_2D, lwd = 2, alpha = 0)
  layer_group_means(p)
  layer_eigen(p, nb_max = 5, cex = 1.3)
  layer_axis_titles(p)
  if (legend) layer_legend(p, cex = 1.15)
  invisible(p)
}

## panel title strip: the letter in bold, n in italics
title_strip <- function(lab, cex = 1.35) {
  par(family = FAM, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(0.015, 0.45, lab, adj = c(0, 0.5), cex = cex, col = COL_INK)
}
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
labs <- list(
  bquote(bold("a)") ~ "  Age 0\u20133 (" * italic(n) * " = " * .(fmt_n(nrow(eA$coe))) * ")"),
  bquote(bold("b)") ~ "  Age 1\u201310 (" * italic(n) * " = " * .(fmt_n(nrow(eB$coe))) * ")"))
pcs <- list(pcA, pcB)

name <- "figA_pca_age0_3_and_age1_10"
W <- 9; H <- 6.5    # inches per panel
HT  <- 0.4          # inches, title strip above each panel
GAP <- 0.2          # inches, white gap between the two panels
tmp <- tempfile("pca_panels_"); dir.create(tmp)
out_both <- function(file, h, fun) {
  cairo_pdf(file.path(tmp, paste0(file, ".pdf")), width = W, height = h, family = FAM)
  fun(); invisible(dev.off())
  ragg::agg_png(file.path(tmp, paste0(file, ".png")), width = W, height = h,
                units = "in", res = 300)
  fun(); invisible(dev.off())
}
for (i in 1:2) {
  out_both(sprintf("t%d", i), HT, function() title_strip(labs[[i]]))
  out_both(sprintf("p%d", i), H,  function() { par(family = FAM); panel(pcs[[i]], legend = i == 1) })
}
out_both("gap", GAP, function() { par(mar = c(0, 0, 0, 0)); plot.new() })

## Pieces, top to bottom. Each panel is drawn on a page of its own and the
## pages are stacked: drawing several PCAs into one device made
## layer_morphospace_PCA put out a second, offset grid of outlines.
pieces <- c("t1", "p1", "gap", "t2", "p2")
magick::image_read(file.path(tmp, paste0(pieces, ".png"))) %>%
  magick::image_append(stack = TRUE) %>%
  magick::image_write(file.path(rev.png, paste0(name, ".png")), format = "png",
                      density = "300x300")

py <- '
import sys
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
out, files = sys.argv[1], sys.argv[2:]
pages = [PdfReader(f).pages[0] for f in files]
w = float(pages[0].mediabox.width)
H = sum(float(p.mediabox.height) for p in pages)
page = PageObject.create_blank_page(width=w, height=H)
y = H
for p in pages:                      # top to bottom
    y -= float(p.mediabox.height)
    page.merge_transformed_page(p, Transformation().translate(0, y))
wr = PdfWriter(); wr.add_page(page); wr.write(out)
'
res <- system2("python3", c("-c", shQuote(py), shQuote(file.path(rev.pdf, paste0(name, ".pdf"))),
                            shQuote(file.path(tmp, paste0(pieces, ".pdf")))))
if (res != 0) stop("stacking the PDF panels failed (needs Python with pypdf)", call. = FALSE)
unlink(tmp, recursive = TRUE)
cat("Done:", fig.dir, "\n")
