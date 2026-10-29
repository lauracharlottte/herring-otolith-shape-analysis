# ============================================================
# Figures 10, 11, 16 and 17 of the thesis
# ------------------------------------------------------------
#   Figure 10  number and proportion of otoliths per population
#   Figure 11  population composition by survey year and ICES Subdivision
#   Figure 16  mean otolith outlines, full data set
#   Figure 17  mean otolith outlines, otoliths aged 1-3 years
#
# Output: results final plots/{pdf,png,data}
# Vector output (PDF) for LaTeX, a 600 dpi PNG for Word, and the
# numbers behind every figure as CSV.
#
# Colours and population order come from config.R (COL_STOCK,
# POP_ORDER). Requires 3 analysis/main_analysis.R to have run.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(Momocs)
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
out.dir  <- ensure_dir(DIR_FINAL)
## results final plots/ is split by format: pdf/ for LaTeX, png/ for
## Word and quick looks, data/ for the numbers behind each figure.
pdf.dir  <- ensure_dir(file.path(out.dir, "pdf"))
png.dir  <- ensure_dir(file.path(out.dir, "png"))
data.dir <- ensure_dir(file.path(out.dir, "data"))

## In-figure titles, subtitles and notes. Off for the thesis: the LaTeX
## caption carries that text, and repeating it inside the image doubles it.
## Switch on for talks or a quick look:
##   FIG_TITLES <- TRUE; source("fig_composition_and_mean_shapes.R")
if (!exists("FIG_TITLES")) FIG_TITLES <- FALSE
tt <- function(x) if (FIG_TITLES) x else NULL   # wraps every title/subtitle/caption

## ---- typography --------------------------------------------
## Every figure is set in FONT_PLOT (Times New Roman, see
## config.R) so the plates match the body text of the thesis.
## theme_set covers axes, strips and legends; geom_text carries
## its own family argument and has to be told separately, hence
## FAM being passed explicitly at every label layer below.
FAM <- FONT_PLOT
theme_set(theme_bw(base_size = 11, base_family = FAM))
cat("figure font:", FAM, "\n")

## Save one figure in two formats under the same basename.
##   .pdf  the deliverable: vector, so it stays sharp at any size,
##         with Times New Roman embedded. Goes straight into LaTeX
##         via \includegraphics and opens in Preview.
##   .png  600 dpi raster, for anything that handles PDF badly --
##         Word, PowerPoint or a quick look.
##         ragg is used because the stock png() device on macOS
##         does not resolve system font families.
save_fig <- function(plot, name, width, height) {
  pdir <- pdf.dir
  gdir <- png.dir
  ggsave(file.path(pdir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "cm", device = cairo_pdf)
  ggsave(file.path(gdir, paste0(name, ".png")), plot,
         width = width, height = height, units = "cm", dpi = 600,
         device = ragg::agg_png)
  cat("  ", name, ".pdf / .png\n", sep = "")
}

# ============================================================
# Figure 11 -- population composition by year and ICES Subdivision
# ============================================================
# 100 % stacked bars, one panel per ICES subdivision, three years
# per panel, the cell n printed above every bar.
#
# Percentages are within a year x subdivision cell, never within
# a year alone: sampling effort moved between subdivisions across
# the three surveys (SD23 went 67 -> 476 -> 289 otoliths), so a
# year-pooled percentage would mostly report where the survey
# sampled rather than how the stocks changed.
# ------------------------------------------------------------
cat("Figure 11: population composition by year and SD\n")

## Population order from config.R (POP_ORDER): WBSS sits on the
## baseline, then BAS, NSAS, Downs and CBH on top.
STACK <- STOCK_ORDER

## The analysed sample. NEATL (n = 2) and NSS (n = 3) are excluded
## from the analysis itself in main_analysis.R, so dataBIO
## holds exactly the five populations -- checked here.
load(file.path(RESULTS, "dataBIO.RData"))   # data
stopifnot(all(data$Assigned.Stock %in% STACK))

comp <- data %>%
  transmute(Year  = as.integer(Year),
            SD    = as.integer(Subdivision),
            Stock = factor(Assigned.Stock, levels = STACK)) %>%
  count(Year, SD, Stock, .drop = FALSE) %>%
  group_by(Year, SD) %>%
  mutate(cell_n = sum(n), percent = 100 * n / cell_n) %>%
  ungroup()

cat("  n =", sum(comp$n), "\n")

bars <- comp %>%
  filter(Stock %in% STACK) %>%
  mutate(Stock = factor(Stock, levels = STACK),
         SD    = factor(SD, labels = paste("SD", sort(unique(SD)))))

n_lab <- distinct(bars, Year, SD, cell_n)

fig1 <- ggplot(bars, aes(factor(Year), percent, fill = Stock)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.5,
           position = position_stack(reverse = TRUE)) +
  ## italic n, as statistical notation requires; parse = TRUE turns
  ## the string into a plotmath expression rather than literal text
  geom_text(data = n_lab,
            aes(factor(Year), 103, label = paste0("italic(n)==", cell_n)),
            parse = TRUE, inherit.aes = FALSE, size = 2.9,
            colour = COL_MUT, vjust = 0, family = FAM) +
  facet_wrap(~SD, nrow = 1) +
  scale_fill_manual(values = COL_STOCK[STACK], name = NULL) +
  scale_y_continuous(limits = c(0, 111), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "Proportion of otoliths within each\nyear and ICES Subdivision") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
    panel.border       = element_rect(colour = "grey80", linewidth = 0.4),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", colour = COL_INK, size = 10.5),
    axis.text          = element_text(colour = COL_INK),
    axis.title.y       = element_text(colour = COL_INK, margin = margin(r = 8)),
    legend.position    = "bottom",
    legend.key.size    = unit(4, "mm"),
    legend.text        = element_text(colour = COL_INK)
  )

save_fig(fig1, "fig11_composition_by_year_and_SD", width = 18, height = 8.5)

## the numbers behind the figure, so the caption can be checked
write.csv(
  comp %>% mutate(percent = round(percent, 1)) %>%
    select(Year, SD, Stock, n, percent, cell_n),
  file.path(data.dir, "fig11_composition_by_year_and_SD_data.csv"), row.names = FALSE)


# ============================================================
# Figure 10 -- number and proportion of otoliths per population
# ============================================================
# The same analysed otoliths without the year and subdivision split:
# one bar per stock, so the headline composition can be read off
# directly instead of being reconstructed from twelve cells.
#
# Figure 11 and Table 2 carry the year x subdivision detail;
# this figure shows at a glance how unbalanced the material is
# overall, which the class-imbalance analysis rests on.
#
# Horizontal bars, not a single stacked bar: with one bar per
# stock every value sits on a common baseline and can be
# compared directly.
# ------------------------------------------------------------
cat("Figure 10: overall composition\n")

## Percentages over the analysed sample, the same base as
## Figure 11 and Table 2.
overall_all <- data %>%
  count(Stock = Assigned.Stock, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

overall <- overall_all %>%
  mutate(Stock = factor(Stock, levels = rev(STOCK_ORDER))) %>%
  arrange(Stock)

cat("  n =", sum(overall$n), "in", nrow(overall), "populations\n")

fill_overall <- COL_STOCK

fig2 <- ggplot(overall, aes(percent, Stock, fill = Stock)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%s (%.1f%%)", formatC(n, format = "d", big.mark = ","), percent)),
            hjust = -0.1, size = 3.1, colour = COL_INK, family = FAM) +
  scale_fill_manual(values = fill_overall, guide = "none") +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0, 0)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = NULL,
       title = tt("Proportion of otoliths per population")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
    panel.border       = element_rect(colour = "grey80", linewidth = 0.4),
    axis.text.y        = element_text(colour = COL_INK, face = "bold", size = 10),
    axis.text.x        = element_text(colour = COL_INK),
    axis.ticks.y       = element_blank(),
    plot.title         = element_text(colour = COL_INK, size = 11,
                                      hjust = 0, margin = margin(b = 9)),
    plot.title.position = "plot"
  )

save_fig(fig2, "fig10_overall_composition", width = 13, height = 7)

## the numbers behind the figure
write.csv(overall_all %>% mutate(percent = round(percent, 1)) %>% arrange(desc(n)),
          file.path(data.dir, "fig10_overall_composition_data.csv"), row.names = FALSE)



# ============================================================
# Figures 16 and 17 -- mean otolith outlines
# ============================================================
# The mean outline of every population drawn against the mean
# outline of WBSS (filled), for the full data set and for the
# otoliths aged 1-3 years.

## Mean outline of every population over the reference (WBSS) as a
## filled template. Outlines closed -- Momocs leaves out the last segment.
mean_shape_figure <- function(eO, file, ref = "WBSS", title = NULL) {
  ms  <- MSHAPES(eO, fac = "Assigned.Stock")$shp
  lv  <- intersect(POP_ORDER, names(ms))
  shp <- bind_rows(lapply(lv, function(g) {
    s <- as.data.frame(ms[[g]]); s <- rbind(s, s[1, ]); s$Group <- g; s }))
  shp$Group <- factor(shp$Group, levels = lv)
  oth <- setdiff(lv, ref)
  pd  <- bind_rows(lapply(oth, function(g) bind_rows(
           dplyr::filter(shp, Group == ref) %>% dplyr::mutate(Panel = g),
           dplyr::filter(shp, Group == g)   %>% dplyr::mutate(Panel = g)))) %>%
    dplyr::mutate(Panel = factor(Panel, levels = oth))
  f <- ggplot(pd, aes(x, y)) +
    geom_polygon(data = ~ dplyr::filter(.x, Group == ref),
                 fill = COL_STOCK[[ref]], alpha = 0.18, colour = NA) +
    geom_path(aes(colour = Group), linewidth = 0.6) +
    facet_wrap(~Panel, ncol = 2,
               labeller = as_labeller(setNames(paste(ref, "and", oth), oth))) +
    scale_colour_manual(values = COL_STOCK[lv], limits = lv, name = "Population") +
    coord_equal() +
    labs(title = tt(title)) +
    theme_void(base_size = 11, base_family = FAM) +
    theme(legend.position = "bottom",
          legend.key.width = unit(8, "mm"),
          legend.title  = element_text(colour = COL_INK, size = 10),
          legend.text   = element_text(colour = COL_INK, size = 10),
          strip.text    = element_text(colour = COL_INK, size = 10.5,
                                       margin = margin(b = 3)),
          panel.spacing = unit(4, "mm"),
          plot.title    = element_text(colour = COL_INK, size = 11, hjust = 0.5,
                                       margin = margin(b = 6)),
          plot.margin   = margin(6, 6, 6, 6))
  save_fig(f, file, width = 17, height = 11)
  write.csv(shp %>% dplyr::rename(Population = Group) %>%
              dplyr::mutate(x = round(x, 5), y = round(y, 5)),
            file.path(data.dir, paste0(file, "_data.csv")), row.names = FALSE)
  invisible(f)
}

## -- Figure 16: full data set -----------------------------------------
cat("Figure 16: mean shapes against WBSS, full data set\n")
e6 <- local({ e <- new.env(); load(file.path(RESULTS, "eOTO.RData"), envir = e); e$eOTO })
stopifnot(nrow(e6$coe) == nrow(data))
mean_shape_figure(e6, "fig16_mean_shapes", title = "WBSS (filled) against each population")

## -- Figure 17: age 1-3 -----------------------------------------------
## Fourier coefficients are computed per otolith, so the subset is the
## same rows of eOTO.
cat("Figure 17: mean shapes against WBSS, age 1-3\n")
age6 <- suppressWarnings(as.numeric(as.character(e6$fac$Age)))
e8   <- Momocs::slice(e6, which(!is.na(age6) & age6 >= 1 & age6 <= 3))
stopifnot(nrow(e8$coe) == sum(!is.na(data$Age) & data$Age >= 1 & data$Age <= 3))
cat("  n =", nrow(e8$coe), "\n")
mean_shape_figure(e8, "fig17_mean_shapes_age1_3",
                  title = "WBSS (filled) against each population, age 1-3")

cat("\nDone. Figures in:", out.dir, "\n")
