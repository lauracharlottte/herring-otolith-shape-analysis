# ============================================================
# Figures A2-A4 -- confusion matrices of all LDA models in the thesis
# ------------------------------------------------------------
# Every model that appears in the classification tables of the
# Results chapter, drawn as a confusion matrix:
#
#   figA_confusion_five_populations  a full data set, b age 0-3,
#                                    c age 1-3, d age 1-10
#   figA_confusion_imbalance         a observed priors (= full data set),
#                                    b equal priors, c downsampling,
#                                    d upsampling
#   figA_confusion_two_groups        a WBSS vs rest, b BAS vs rest,
#                                    c spawning groups
#
# Source: caret::confusionMatrix(fit)$table of the saved fits in
# Fourier and Classification/Results -- the cross-validated
# matrix averaged over all 400 runs (4-fold CV, 100 repeats), in
# per cent of all otoliths. Columns = genetic (true) population,
# rows = population assigned by the LDA.
#
# Cells show the percentage of the otoliths of a population (column)
# that were assigned to each population (row), so every column adds
# up to 100 and the diagonal is "correctly classified" in the
# Results tables. The small number below is the mean number of
# otoliths per run. Colour = column percentage.
#
# Output: Shape Analysis/Appendix/{pdf,png}/ and
#         Appendix/tables/confusion_matrices.csv (all cells, long format)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(caret); library(grid)
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
app.dir <- ensure_dir(file.path(SHAPE_DIR, "Appendix"))
app.pdf <- ensure_dir(file.path(app.dir, "pdf"))
app.png <- ensure_dir(file.path(app.dir, "png"))
app.tab <- ensure_dir(file.path(app.dir, "tables"))
FAM <- FONT_PLOT

## ---- the models ------------------------------------------------
## n_eval = otoliths in the held-out folds of one run; for down- and
## upsampling the test folds were not resampled, so it is 3,029 too.
load_fit <- function(path) { e <- new.env(); load(path, envir = e); e[[ls(e)[1]]] }
models <- list(
  full     = list(file = "00_full/LDA_results.RData",            title = "Full data set"),
  age0_3   = list(file = "age0_3/LDA_results.RData",             title = "Age 0–3"),
  age1_3   = list(file = "age1_3/LDA_results.RData",             title = "Age 1–3"),
  age1_10  = list(file = "age1_10/LDA_results.RData",            title = "Age 1–10"),
  observed = list(file = "00_full/LDA_results.RData",            title = "Observed priors"),
  equal    = list(file = "Imbalance/fit_Equal_priors.RData",     title = "Equal priors"),
  down     = list(file = "Imbalance/fit_Downsampling.RData",     title = "Downsampling"),
  up       = list(file = "Imbalance/fit_Upsampling.RData",       title = "Upsampling"),
  wbss     = list(file = "WBSS_vs_Rest/LDA_results.RData",       title = "WBSS vs. rest"),
  bas      = list(file = "BAS_vs_Rest/LDA_results.RData",        title = "BAS vs. rest"),
  spawn    = list(file = "SpawningComponent/LDA_results.RData",  title = "Spawning groups")
)
N_FULL <- 3029

SPAWN_LAB <- c("spring spawner" = "Spring", "autumn spawner" = "Autumn",
               "winter spawner" = "Winter")

long <- list()
for (k in names(models)) {
  fit <- load_fit(file.path(RESULTS, models[[k]]$file))
  tab <- unclass(caret::confusionMatrix(fit)$table)       # % of all, rows = prediction
  lv  <- rownames(tab)
  ord <- if (all(lv %in% POP_ORDER)) intersect(POP_ORDER, lv)
         else if (all(lv %in% names(SPAWN_LAB))) names(SPAWN_LAB)
         else c(setdiff(lv, "Rest"), "Rest")
  tab <- tab[ord, ord]
  n_eval <- if (k %in% c("down", "up")) N_FULL else nrow(fit$trainingData)
  colpct <- sweep(tab, 2, colSums(tab), "/") * 100
  lab <- function(x) if (all(x %in% names(SPAWN_LAB))) unname(SPAWN_LAB[x]) else x
  df <- expand.grid(assigned = ord, true = ord, stringsAsFactors = FALSE)
  df$pct_of_all  <- tab[cbind(df$assigned, df$true)]
  df$pct_of_true <- colpct[cbind(df$assigned, df$true)]
  df$n_mean      <- df$pct_of_all * n_eval / 100
  df$assigned <- factor(lab(df$assigned), levels = rev(lab(ord)))
  df$true     <- factor(lab(df$true),     levels = lab(ord))
  df$model <- models[[k]]$title; df$n <- n_eval
  df$accuracy <- sum(diag(tab))
  long[[k]] <- df
}

## ---- one panel -----------------------------------------------
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
## blue ramp from the thesis palette: white -> light blue -> NSAS blue -> BAS navy
## (explicit stops, so the blend stays blue instead of turning grey-violet)
RAMP <- c("#FFFFFF", "#DCE6F8", "#93AEEA", "#2B5BD1", "#07306B")
panel <- function(df, letter, base = 11) {
  df$txt_col <- ifelse(df$pct_of_true > 50, "white", COL_INK)
  diag_df <- df[as.character(df$assigned) == as.character(df$true), ]
  big <- nlevels(df$true) <= 3
  ggplot(df, aes(true, assigned)) +
    geom_tile(aes(fill = pct_of_true), colour = "grey80", linewidth = 0.3) +
    geom_tile(data = diag_df, fill = NA, colour = COL_INK, linewidth = 1.2) +
    annotate("rect", xmin = 0.5, xmax = nlevels(df$true) + 0.5,
             ymin = 0.5, ymax = nlevels(df$true) + 0.5,
             fill = NA, colour = COL_INK, linewidth = 1.2) +
    geom_text(aes(label = sprintf("%.1f", pct_of_true), colour = txt_col),
              family = FAM, size = if (big) 4.4 else 3.7, vjust = -0.15) +
    geom_text(aes(label = sprintf("(%s)", fmt_n(round(n_mean))), colour = txt_col),
              family = FAM, size = if (big) 3.1 else 2.6, vjust = 1.35) +
    scale_fill_gradientn(colours = RAMP, limits = c(0, 100), guide = "none") +
    scale_colour_identity() +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(clip = "off") +   # bold edges at the border are not cut in half
    labs(x = "Genetic population", y = "Assigned by the LDA",
         title = bquote(bold(.(letter)) ~ " " * .(df$model[1]) ~
                          "(" * italic(n) * " = " * .(fmt_n(df$n[1])) * ")")) +
    theme_minimal(base_size = base, base_family = FAM) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = COL_INK),
          axis.title = element_text(colour = COL_MUT, size = base - 1),
          plot.title = element_text(colour = COL_INK, size = base, hjust = 0,
                                    margin = margin(b = 4)),
          plot.title.position = "plot",
          plot.margin = margin(6, 10, 6, 6))
}

## ---- figures: panels placed with grid viewports ------------------
draw_grid <- function(plots, nrow, ncol) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  for (i in seq_along(plots)) {
    r <- (i - 1) %/% ncol + 1; c <- (i - 1) %% ncol + 1
    print(plots[[i]], vp = viewport(layout.pos.row = r, layout.pos.col = c))
  }
}
save_fig <- function(name, plots, nrow, ncol, w, h) {
  cairo_pdf(file.path(app.pdf, paste0(name, ".pdf")), width = w, height = h, family = FAM)
  draw_grid(plots, nrow, ncol); invisible(dev.off())
  ragg::agg_png(file.path(app.png, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 300)
  draw_grid(plots, nrow, ncol); invisible(dev.off())
}
L <- c("a)", "b)", "c)", "d)")
save_fig("figA_confusion_five_populations",
         Map(panel, long[c("full", "age0_3", "age1_3", "age1_10")], L), 2, 2, 9, 9)
save_fig("figA_confusion_imbalance",
         Map(panel, long[c("observed", "equal", "down", "up")], L), 2, 2, 9, 9)
save_fig("figA_confusion_two_groups",
         Map(panel, long[c("wbss", "bas", "spawn")], L[1:3]), 1, 3, 11, 4.2)

## ---- all cells as one table --------------------------------------
out <- do.call(rbind, lapply(long[setdiff(names(long), "observed")], function(d)
  data.frame(model = d$model, n = d$n, true_population = d$true,
             assigned_population = d$assigned,
             percent_of_true_population = round(d$pct_of_true, 2),
             percent_of_all_otoliths = round(d$pct_of_all, 3),
             mean_number_of_otoliths = round(d$n_mean, 1),
             accuracy_percent = round(d$accuracy, 2))))
write.csv(out, file.path(app.tab, "confusion_matrices.csv"), row.names = FALSE)
cat("Done:", app.dir, "\n")
