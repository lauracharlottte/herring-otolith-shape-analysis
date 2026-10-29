# ============================================================
# Fig. A2 -- sample flow from the collected herring to the analysed otoliths
# ------------------------------------------------------------
# The numbers are written into the figure by hand; they are the counts
# reported in Sections 2.1-2.6 and 3.1-3.3 (sample composition from
# stock_composition.R, outliers from main_analysis.R).
# Usage: Rscript appendix_sample_flow.R [output.png]
# ============================================================
library(grid)

out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- "appendix_sample_flow.png"

navy  <- "#07306B"
grey  <- "#5A5A5A"
fill_main <- "#E8EDF5"
fill_excl <- "#F3F3F3"
fam <- "Times New Roman"

main <- list(
  c("3,411 herring",  "collected during GERAS 2023–2025"),
  c("3,383 herring",  "with a genetic population assignment"),
  c("3,163 outlines", "of usable right otoliths"),
  c("3,158 outlines", "of the five analysed populations"),
  c("3,029 otoliths", "analysed (WBSS, BAS, NSAS, Downs, CBH)")
)
excl <- list(
  c("− 28", "genetic assignment failed or ambiguous"),
  c("− 220", "no usable right otolith: only the",
    "left otolith photographed, otolith damaged",
    "or outline not extracted (Section 2.4)"),
  c("− 5", "assigned to NEATL (2) or NSS (3);",
    "too rare for the analyses (Section 2.2)"),
  c("− 129", "outliers (D² > 20.5; Section 2.6)")
)

W <- 16 / 2.54; H <- 13 / 2.54
ragg::agg_png(out, width = W, height = H, units = "in", res = 600)
grid.newpage()
pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 13)))

mx <- 4.3; bw <- 7.0; bh <- 1.55          # main column
ex <- 12.2; ew <- 7.0                      # exclusion column
ys <- seq(12.0, 1.0, length.out = 5)

box <- function(x, y, w, h, fill, col) {
  grid.roundrect(unit(x, "native"), unit(y, "native"),
                 unit(w, "native"), unit(h, "native"),
                 r = unit(0.08, "in"),
                 gp = gpar(fill = fill, col = col, lwd = 1.4))
}
txt <- function(lab, x, y, size, col, face = "plain", just = "centre") {
  grid.text(lab, unit(x, "native"), unit(y, "native"), just = just,
            gp = gpar(fontfamily = fam, fontsize = size, col = col,
                      fontface = face))
}
arrow_ <- function(x0, y0, x1, y1, col) {
  grid.segments(unit(x0, "native"), unit(y0, "native"),
                unit(x1, "native"), unit(y1, "native"),
                arrow = arrow(length = unit(0.07, "in"), type = "closed"),
                gp = gpar(col = col, fill = col, lwd = 1.3))
}

for (i in seq_along(main)) {
  box(mx, ys[i], bw, bh, fill_main, navy)
  txt(main[[i]][1], mx, ys[i] + 0.28, 12, navy, "bold")
  txt(main[[i]][2], mx, ys[i] - 0.30, 10, navy)
  if (i < length(main)) {
    arrow_(mx, ys[i] - bh / 2, mx, ys[i + 1] + bh / 2, navy)
    # exclusion branch from the midpoint of the arrow
    ym <- (ys[i] + ys[i + 1]) / 2
    e  <- excl[[i]]
    eh <- 0.55 + 0.38 * (length(e) - 1)
    arrow_(mx, ym, ex - ew / 2, ym, grey)
    box(ex, ym, ew, eh, fill_excl, grey)
    lines <- c(paste0(e[1], "  ", e[2]), e[-(1:2)])
    yy <- ym + 0.19 * (length(lines) - 1) - 0.38 * (seq_along(lines) - 1)
    for (k in seq_along(lines))
      txt(lines[k], ex - ew / 2 + 0.25, yy[k], 9.5, grey,
          just = "left")
  }
}
invisible(dev.off())
