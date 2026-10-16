# ============================================================
# Additional analyses
# ------------------------------------------------------------
# main_analysis.R covers the full model, three age windows, the
# two one-vs-rest models and the class-imbalance experiment.
# This script runs the remaining analyses on exactly the same
# otoliths, so every number in the results shares one data set:
#
#   Calibration/     harmonic power + reconstruction plots
#   SpawningComponent/
#   no_CBH/
#   age1_3_maturity12/
#   Harmonics_18/    12 vs 18, paired on identical folds
#
# The harmonics sweep (2-35) is a separate script because it
# takes minutes rather than seconds: harmonics_sweep.R
#
# Everything is read from Results/, so the sample is by
# construction the same one the main model used.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(Momocs); library(caret); library(MASS); library(doParallel)
})
select <- dplyr::select
theme_set(theme_bw())

NB_H       <- 12
CV_FOLDS   <- 4
CV_REPEATS <- 100
SEED       <- 2021
PLOT_DPI   <- 600

## ---- load project paths ------------------------------------
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

RESULTS <- DIR_RESULTS

# ------------------------------------------------------------
# Colours and font come from config.R (COL_STOCK, COL_SPAWN,
# COL_BAR, FONT_PLOT) -- the same scheme as in all thesis figures.
# ------------------------------------------------------------
theme_set(theme_bw(base_family = FONT_PLOT))
COL_UNKNOWN <- "#07306B"

cols_for <- function(lv, map, fallback = COL_UNKNOWN) {
  out <- unname(map[lv]); out[is.na(out)] <- fallback; out
}
pal_fixed <- function(cols) function(n) cols[seq_len(n)]

# ------------------------------------------------------------
# The main-model sample (five stocks, outliers removed), exactly
# as main_analysis.R left it in dataBIO.RData
# ------------------------------------------------------------
load(file.path(RESULTS, "dataBIO.RData"))     # data
load(file.path(RESULTS, "dataSHAPE.RData"))   # 45 descriptors
load(file.path(RESULTS, "oto_all.RData"))     # oto, still with outliers

oto  <- oto[names(oto) %in% data$ID]
data <- data %>% arrange(match(ID, names(oto)))
stopifnot(nrow(data) == length(oto), nrow(data) == nrow(dataSHAPE))
data$Age <- suppressWarnings(as.numeric(as.character(data$Age)))
cat("base sample:", nrow(data), "otoliths\n\n")

slide_to_bottom <- function(coo) coo_slide(coo, id = which.min(coo[, 2]))
build_out <- function(o, d) {
  O <- Out(o, fac = d); O$coo <- lapply(O$coo, slide_to_bottom)
  O <- Out(O$coo, fac = d); O$fac <- droplevels(O$fac); O
}

cl <- makeCluster(detectCores()); registerDoParallel(cl)
trcntr <- trainControl(method = "repeatedcv", number = CV_FOLDS,
                       repeats = CV_REPEATS, allowParallel = TRUE)

metrics <- function(fit) {
  tb <- confusionMatrix(fit)$table
  dg <- diag(tb); cs <- colSums(tb); rs <- rowSums(tb)
  list(acc = mean(fit$resample$Accuracy), accSD = sd(fit$resample$Accuracy),
       kap = mean(fit$resample$Kappa),    kapSD = sd(fit$resample$Kappa),
       sens = 100 * dg / cs, fdr = 100 * (1 - dg / rs), tb = tb)
}

# ============================================================
# 1. Calibration plots on the analysed sample
# ============================================================
cat("=== 1. calibration ===\n")
cal.dir <- ensure_dir(file.path(RESULTS, "Calibration"))
OTO_all <- build_out(oto, data)

set.seed(SEED)
harmpow <- calibrate_harmonicpower_efourier(x = OTO_all, nb.h = 30)
cat("harmonics for 90/95/99/99.9 % shape power:",
    paste(harmpow$minh, collapse = " / "), "\n")
ggsave(file.path(cal.dir, "calibrate_harmonicpower_efourier.tiff"),
       plot = harmpow$gg, width = 10, height = 5,
       dpi = PLOT_DPI, compression = "lzw")

set.seed(SEED)
rec <- calibrate_reconstructions_efourier(OTO_all, range = c(1, 4, 8, 10, NB_H, 20))
## Momocs labels these panels 1..6 instead of the harmonic counts
rec$data$id <- as.factor(rec$data$id)
levels(rec$data$id) <- c(1, 4, 8, 10, NB_H, 20)
ggsave(file.path(cal.dir, "calibrate_reconstructions_efourier.tiff"),
       plot = rec, width = 10, height = 5, dpi = PLOT_DPI, compression = "lzw")

## cumulative shape power at 12 harmonics, for the caption
hp <- harmpow$q
writeLines(c(
  paste("Sample                :", nrow(data)),
  paste("Harmonics for 90 %    :", harmpow$minh[1]),
  paste("Harmonics for 95 %    :", harmpow$minh[2]),
  paste("Harmonics for 99 %    :", harmpow$minh[3]),
  paste("Harmonics for 99.9 %  :", harmpow$minh[4]),
  paste("Used in the thesis    :", NB_H)),
  file.path(cal.dir, "harmonic_power.txt"))
cat("  written\n\n")

# ============================================================
# 2. Generic sub-analysis runner (LDA + CV + all four figures)
# ============================================================
run_one <- function(folder, idx, grp, title, pca_name = "PCA",
                    focal = NULL, spawn_panel = TRUE) {

  dir <- ensure_dir(file.path(RESULTS, folder))
  d   <- data[idx, , drop = FALSE]
  grp <- droplevels(factor(grp[idx]))
  lv  <- levels(grp)
  cat("---", folder, "| n =", nrow(d), "| classes:",
      paste(lv, collapse = ", "), "\n")

  ## Colours by what the groups are: spawning components from COL_SPAWN,
  ## populations from COL_STOCK. Before 2026-09-11 the spawning
  ## components were looked up in COL_STOCK, found nothing, and all
  ## three fell back to the same navy in the mean-shape and PCA plots.
  cols_grp <- if (!is.null(focal)) ifelse(lv == focal, COL_STOCK[[focal]], "#2271B5") else
              if (all(lv %in% names(COL_SPAWN))) cols_for(lv, COL_SPAWN) else
              cols_for(lv, COL_STOCK)
  names(cols_grp) <- lv

  ## EFD recomputed on the subset, as every subset script does.
  ## The figures show populations in POP_ORDER (WBSS first), but only
  ## through this plotting copy: `grp` itself goes into the model
  ## unchanged, because caret draws its folds class by class in level
  ## order and a reordered factor would change the cross-validation.
  d$Group <- if (all(lv %in% POP_ORDER)) factor(grp, levels = intersect(POP_ORDER, lv)) else grp
  OTO  <- build_out(oto[idx], d)
  eOTO <- efourier(OTO, nb.h = NB_H, smooth.it = 0, norm = TRUE)
  sh   <- as.data.frame(eOTO$coe) %>% select(-A1, -B1, -C1)

  ## -- mean shapes --
  ## Populations in POP_ORDER (WBSS first); spawning components keep
  ## their order. Outlines closed -- Momocs leaves a one-step gap.
  ms <- MSHAPES(eOTO, fac = "Group")$shp
  is_pop <- all(names(ms) %in% POP_ORDER)
  lv_ms  <- if (is_pop) intersect(POP_ORDER, names(ms)) else intersect(lv, names(ms))
  leg_ms <- if (is_pop) "Population" else
            if (all(names(ms) %in% names(COL_SPAWN))) "Spawning component" else "Group"
  dms <- bind_rows(lapply(lv_ms, function(g) {
    sh <- as.data.frame(ms[[g]]); sh <- rbind(sh, sh[1, ]); sh$Group <- g; sh }))
  dms$Group <- factor(dms$Group, levels = lv_ms)
  ggplot(dms) +
    geom_path(aes(x, y, color = Group, linetype = Group), lwd = 1) +
    scale_color_manual(values = cols_grp) +
    labs(color = leg_ms, linetype = leg_ms) +
    theme_void(base_family = FONT_PLOT) + coord_equal() +
    theme(legend.position = "inside", legend.position.inside = c(0.6, 0.6),
          legend.key.size = unit(5, "mm"))
  ggsave(file.path(dir, "mshapes.tiff"), width = 10, height = 5,
         dpi = PLOT_DPI, compression = "lzw")

  ## -- PCA, Momocs morphospace style --
  pca_ <- PCA(eOTO)
  panels <- c(pca_name)
  if (spawn_panel) panels <- c(panels, "PCA_SpringvsAutumn")
  for (fn in panels) {
    grpcol <- if (fn == pca_name) "Group" else "Spawning.Component"
    ttl <- if (fn == pca_name) paste("PCA --", title) else
      paste("PCA --", title, "-- spawning component")
    lv_i <- levels(droplevels(factor(OTO$fac[[grpcol]])))
    cols_i <- if (fn == pca_name) unname(cols_grp[lv_i]) else cols_for(lv_i, COL_SPAWN)
    tiff(file.path(dir, paste0(fn, ".tiff")), width = 12, height = 10,
         units = "in", res = PLOT_DPI, compression = "lzw")
    par(family = FONT_PLOT)
    p <- plot_PCA(pca_, as.formula(paste0("~", grpcol)), zoom = 2.2,
                  palette = pal_fixed(cols_i), chull = FALSE, eigen = FALSE,  # inset via layer_eigen()
                  morphospace = FALSE, axesnames = FALSE, axesvar = FALSE,
                  points = FALSE)
    layer_morphospace_PCA(p, position = "range", nb = 12, nr = 12, nc = 10,
                          rotate = 0, size = 0.8, col = "#999999",
                          flipx = FALSE, flipy = FALSE, draw = TRUE)
    layer_axesvar(p, cex = 1.5); layer_stars(p, alpha = 0.8)
    layer_points(p, cex = 0.8);  layer_eigen(p, nb_max = 5, cex = 1)
    layer_title(p, title = ttl, cex = 1)
    layer_ellipses(p, conf = 0.9, lwd = 1, alpha = 0)
    dev.off()
  }

  ## -- LDA projection + variable importance --
  dcl <- cbind(Group = grp, sh)
  lin <- lda(Group ~ ., dcl)
  tiff(file.path(dir, "LDA.tiff"), width = 10, height = 8, units = "in",
       res = PLOT_DPI, compression = "lzw")
  par(family = FONT_PLOT)
  plot(lin, col = unname(cols_grp)[as.integer(grp)])
  dev.off()

  vi <- data.frame(Variable = rownames(lin$scaling),
                   Importance = rowSums(lin$scaling^2)) %>% arrange(desc(Importance))
  write.csv(vi, file.path(dir, "variable_importance.csv"), row.names = FALSE)
  ggplot(head(vi, 20), aes(reorder(Variable, Importance), Importance)) +
    geom_bar(stat = "identity", fill = COL_BAR) + coord_flip() +
    labs(x = "Variable", y = "Importance", title = "LDA variable importance")
  ggsave(file.path(dir, "Importance.tiff"), width = 6, height = 8,
         dpi = PLOT_DPI, compression = "lzw")

  ## -- cross-validated LDA --
  set.seed(SEED)
  fit <- train(Group ~ ., data = dcl, method = "lda",
               preProcess = c("scale", "center"), trControl = trcntr)
  m <- metrics(fit)
  save(fit, file = file.path(dir, "LDA_results.RData"))

  sink(file.path(dir, "LDA_results.txt"))
  cat("===", title, "=== (n =", nrow(d), ")\n\n"); print(fit)
  cat("\n=== Confusion Matrix ===\n"); print(confusionMatrix(fit))
  cat("\nClass sizes:\n"); print(table(grp))
  cat("\nCorrectly classified / sensitivity (%):\n"); print(round(m$sens, 1))
  cat("\nFalsely assigned / false discovery rate (%):\n"); print(round(m$fdr, 1))
  sink()

  cat("    accuracy", round(m$acc, 4), "| kappa", round(m$kap, 4), "\n")
  c(m, list(label = title, n = nrow(d), k = nlevels(grp)))
}

res <- list()

# ============================================================
# 3. Spawning component
# ============================================================
cat("=== 2. spawning component ===\n")
sc_ok <- !is.na(data$Spawning.Component) & nzchar(trimws(data$Spawning.Component))
res$spawn <- run_one("SpawningComponent", which(sc_ok),
                     data$Spawning.Component, "Spawning component",
                     pca_name = "PCA_SpawningComponent", spawn_panel = FALSE)
capture.output({
  cat("=== Spawning component, n =", sum(sc_ok), "===\n\n")
  print(table(data$Spawning.Component[sc_ok]))
  cat("\nSpawning component by assigned stock:\n")
  print(table(data$Spawning.Component[sc_ok], data$Assigned.Stock[sc_ok]))
}, file = file.path(RESULTS, "SpawningComponent", "summary.txt"))
cat("\n")

# ============================================================
# 4. CBH excluded
# ============================================================
cat("=== 3. CBH excluded ===\n")
res$noCBH <- run_one("no_CBH", which(data$Assigned.Stock != "CBH"),
                     data$Assigned.Stock, "CBH excluded")
cat("\n")

# ============================================================
# 5. Age 1-3 restricted to maturity stages 1-2
# ============================================================
cat("=== 4. age 1-3 x maturity 1-2 ===\n")
mat <- read.csv(file.path(DIR_DATASETS, "Otolith_Shape_DATA_SET_with_Age_Maturity.csv"),
                fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE) %>%
  transmute(Year = as.integer(Year), Subdivision = as.integer(Subdivision),
            IdentNr.Picture = as.integer(IdentNr.Picture),
            Maturity = suppressWarnings(as.numeric(as.character(Maturity)))) %>%
  filter(!is.na(IdentNr.Picture)) %>%
  distinct(Year, Subdivision, IdentNr.Picture, .keep_all = TRUE)

data$Maturity <- data %>%
  select(Year, Subdivision, IdentNr.Picture) %>%
  left_join(mat, by = c("Year", "Subdivision", "IdentNr.Picture")) %>%
  pull(Maturity)
cat("maturity matched for", sum(!is.na(data$Maturity)), "of", nrow(data), "\n")

idx_am <- which(!is.na(data$Age) & data$Age >= 1 & data$Age <= 3 &
                !is.na(data$Maturity) & data$Maturity %in% c(1, 2))
res$agemat <- run_one("age1_3_maturity12", idx_am, data$Assigned.Stock,
                      "Age 1-3, maturity 1-2")
cat("\n")

stopCluster(cl); registerDoSEQ()

# ============================================================
# 6. Summary table for the three new analyses
# ============================================================
pick <- function(v, nm) if (nm %in% names(v)) unname(v[nm]) else NA_real_
t3 <- bind_rows(lapply(res, function(r) data.frame(
  Analysis = r$label, n = r$n, k = r$k,
  Accuracy = 100 * r$acc, AccuracySD = 100 * r$accSD,
  Kappa = r$kap, KappaSD = r$kapSD)))
write.csv(t3, file.path(RESULTS, "table3_remaining.csv"), row.names = FALSE)
cat("\n--- Table 3: analyses on the main-model sample (n =", nrow(data), ") ---\n")
print(t3, digits = 4)
cat("\nDone. Everything in:", RESULTS, "\n")
