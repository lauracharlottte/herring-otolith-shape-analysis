# ============================================================
# Figure 12 -- choice of the number of harmonics
# ------------------------------------------------------------
# Cumulative Fourier power and classification accuracy and kappa
# against the number of harmonics, as one figure with three
# vertically stacked panels on a shared x axis:
#
#     a  cumulative Fourier power (boxplots per harmonic rank)
#     b  accuracy
#     c  Cohen's kappa
#
# Output: results final plots/results chapter/{pdf,png,data}
#
# Built as one ggplot with three facets; the facets share
# the x axis, so harmonic 12 sits at the same position in every
# panel.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(Momocs)
})

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

RESULTS   <- DIR_RESULTS
fig.dir <- ensure_dir(file.path(DIR_FINAL, "results chapter"))
rev.pdf <- ensure_dir(file.path(fig.dir, "pdf"))
rev.png <- ensure_dir(file.path(fig.dir, "png"))
rev.dat <- ensure_dir(file.path(fig.dir, "data"))

FAM <- FONT_PLOT
theme_set(theme_bw(base_size = 11, base_family = FAM))

NB_H_CAL <- 30    # denominator: harmonics 2..30, as in the pipeline
H_USED   <- 12    # the value the thesis adopts
THRESH   <- 99

## ---- panel a: cumulative Fourier power ------------------------
load(file.path(RESULTS, "dataBIO.RData"))   # data
load(file.path(RESULTS, "oto_all.RData"))   # oto
oto <- oto[names(oto) %in% data$ID]
cat("otoliths:", length(oto), "\n")

pw  <- t(sapply(oto, function(o)
  harm_pow(Momocs:::efourier.default(o, nb.h = NB_H_CAL))))
cum <- t(apply(pw[, 2:NB_H_CAL], 1, function(x) cumsum(x) / sum(x))) * 100
colnames(cum) <- 2:NB_H_CAL
med <- apply(cum, 2, median)
cat("median at h12:", round(med["12"], 2), "%\n")

## ---- panels b and c: the harmonic sweep -----------------------
## Guard: the sweep must have been run on the analysed sample.
sw_sum <- file.path(dirname(F_SWEEP_CSV), "summary.txt")
sw_n <- if (file.exists(sw_sum)) suppressWarnings(as.integer(sub(".*:\\s*", "",
          grep("^Specimens", readLines(sw_sum), value = TRUE)[1]))) else NA
if (!isTRUE(sw_n == nrow(data)))
  stop("the harmonic sweep was run on n = ", sw_n, ", not on the ", nrow(data),
       " analysed otoliths -- run 3 analysis/harmonics_sweep.R first.", call. = FALSE)
sw <- read.csv(F_SWEEP_CSV) %>% mutate(Accuracy = 100 * Accuracy)

## ---- one long table, one facet per panel ---------------------
## The facet labels double as y-axis titles (strips on the left).
PANELS <- c("Cumulative Fourier\npower (%)", "Accuracy (%)", "Cohen's κ")
pf <- function(i) factor(PANELS[i], levels = PANELS)

box_df <- as.data.frame(cum) %>%
  pivot_longer(everything(), names_to = "h", values_to = "y") %>%
  mutate(h = as.integer(h), Panel = pf(1))
line_df <- bind_rows(
  data.frame(h = sw$Harmonics, y = sw$Accuracy, Panel = pf(2)),
  data.frame(h = sw$Harmonics, y = sw$Kappa,    Panel = pf(3)))
## accuracy green (CBH), kappa strong blue (NSAS)
col_line <- c(COL_STOCK[["CBH"]], COL_STOCK[["NSAS"]])
names(col_line) <- PANELS[2:3]

acc12 <- sw$Accuracy[sw$Harmonics == H_USED]
cat("accuracy at h12:", round(acc12, 2), "%\n")

## y breaks per panel: power in steps of 25, accuracy in steps of 1 %
## so the value at 12 harmonics can be read off, kappa left to ggplot
y_breaks <- function(l) {
  if (diff(l) > 50) return(seq(0, 100, 25))
  if (diff(l) > 2)  return(seq(ceiling(l[1]), floor(l[2]), 1))
  scales::breaks_extended()(l)
}

## panel letters, bottom right of every panel (top left is taken by
## the boxplots and the 99 % line in panel a)
lab_df <- data.frame(Panel = pf(1:3), lab = c("a", "b", "c"),
                     y = c(-Inf, -Inf, -Inf))

x_max <- max(sw$Harmonics)

fig <- ggplot() +
  geom_vline(xintercept = H_USED, colour = COL_INK, linewidth = 0.4) +
  ## a -- 99 % criterion and boxplots
  geom_hline(data = data.frame(Panel = pf(1)), aes(yintercept = THRESH),
             linetype = "22", colour = COL_STOCK[["WBSS"]], linewidth = 0.6) +
  geom_boxplot(data = box_df, aes(h, y, group = h),
               fill = "#DCE6EA", outlier.size = 0.35, outlier.colour = "grey55",
               linewidth = 0.3, width = 0.68) +
  geom_text(data = data.frame(Panel = pf(1)), aes(x = x_max, y = THRESH - 6),
            label = paste0(THRESH, "%"), family = FAM, size = 3.1,
            colour = COL_STOCK[["WBSS"]], hjust = 1) +
  geom_label(data = data.frame(Panel = pf(1)), aes(x = H_USED + 0.4, y = 45),
             label = paste(H_USED, "harmonics"), hjust = 0, family = FAM, size = 3.1,
             colour = COL_INK, fill = "white",
             label.padding = unit(1.2, "mm")) +
  ## b, c -- sweep lines
  geom_line(data = line_df, aes(h, y, colour = Panel), linewidth = 0.6) +
  geom_point(data = line_df, aes(h, y, colour = Panel, shape = Panel), size = 1.6) +
  geom_text(data = lab_df, aes(x = Inf, y = y, label = lab),
            hjust = 1.7, vjust = -0.8, family = FAM, fontface = "bold", size = 3.8,
            colour = COL_INK) +
  scale_colour_manual(values = col_line, guide = "none") +
  scale_shape_manual(values = setNames(c(16, 17), PANELS[2:3]), guide = "none") +
  scale_x_continuous(breaks = seq(2, x_max, 2), limits = c(1.4, x_max + 0.6),
                     expand = expansion(0)) +
  scale_y_continuous(breaks = y_breaks) +
  facet_wrap(~Panel, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = "Number of harmonics", y = NULL) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
        panel.border       = element_rect(colour = "grey80", linewidth = 0.4),
        panel.spacing      = unit(3, "mm"),
        strip.placement    = "outside",
        strip.background   = element_blank(),
        strip.text.y.left  = element_text(angle = 90, colour = COL_INK, size = 11),
        axis.text          = element_text(colour = COL_INK),
        axis.title         = element_text(colour = COL_INK))

name <- "fig_harmonics_power_kappa_accuracy"
ggsave(file.path(rev.pdf, paste0(name, ".pdf")), fig,
       width = 16, height = 18, units = "cm", device = cairo_pdf)
ggsave(file.path(rev.png, paste0(name, ".png")), fig,
       width = 16, height = 18, units = "cm", dpi = 600, device = ragg::agg_png)

## the numbers behind the figure
write.csv(data.frame(harmonic = as.integer(names(med)),
                     median_cum_power = round(med, 3),
                     q25 = round(apply(cum, 2, quantile, .25), 3),
                     q75 = round(apply(cum, 2, quantile, .75), 3)),
          file.path(rev.dat, paste0(name, "_power.csv")), row.names = FALSE)
write.csv(sw, file.path(rev.dat, paste0(name, "_sweep.csv")), row.names = FALSE)

cat("Done:", file.path(fig.dir), "\n")
