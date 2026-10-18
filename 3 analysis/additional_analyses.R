# ============================================================
# Additional analyses
# ------------------------------------------------------------
# main_analysis.R covers the full model, three age windows, the
# two one-vs-rest models and the class-imbalance experiment.
# This script runs the spawning-group model on exactly the same
# otoliths, so every number in the results shares one data set:
#
#   SpawningComponent/   LDA with repeated cross-validation
#   table3_remaining.csv accuracy and kappa of this model
#
# The harmonics sweep (2-35) is a separate script because it
# takes minutes rather than seconds: harmonics_sweep.R
#
# Everything is read from Results/, so the sample is by
# construction the same one the main model used.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(Momocs); library(caret); library(MASS); library(doParallel)
})
select <- dplyr::select

NB_H       <- 12
CV_FOLDS   <- 4
CV_REPEATS <- 100
SEED       <- 2021

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
# The main-model sample (five stocks, outliers removed), exactly
# as main_analysis.R left it in dataBIO.RData
# ------------------------------------------------------------
load(file.path(RESULTS, "dataBIO.RData"))     # data
load(file.path(RESULTS, "dataSHAPE.RData"))   # 45 descriptors
load(file.path(RESULTS, "oto_all.RData"))     # oto, still with outliers

oto  <- oto[names(oto) %in% data$ID]
data <- data %>% arrange(match(ID, names(oto)))
stopifnot(nrow(data) == length(oto), nrow(data) == nrow(dataSHAPE))
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
# 1. Sub-analysis runner (cross-validated LDA)
# ============================================================
run_one <- function(folder, idx, grp, title) {

  dir <- ensure_dir(file.path(RESULTS, folder))
  d   <- data[idx, , drop = FALSE]
  grp <- droplevels(factor(grp[idx]))
  cat("---", folder, "| n =", nrow(d), "| classes:",
      paste(levels(grp), collapse = ", "), "\n")

  ## EFD recomputed on the subset, as every subset script does.
  ## `grp` goes into the model in its original level order, because
  ## caret draws its folds class by class in level order and a
  ## reordered factor would change the cross-validation.
  eOTO <- efourier(build_out(oto[idx], d), nb.h = NB_H, smooth.it = 0, norm = TRUE)
  sh   <- as.data.frame(eOTO$coe) %>% select(-A1, -B1, -C1)
  dcl  <- cbind(Group = grp, sh)

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
# 2. Spawning component
# ============================================================
cat("=== spawning component ===\n")
sc_ok <- !is.na(data$Spawning.Component) & nzchar(trimws(data$Spawning.Component))
res$spawn <- run_one("SpawningComponent", which(sc_ok),
                     data$Spawning.Component, "Spawning component")
capture.output({
  cat("=== Spawning component, n =", sum(sc_ok), "===\n\n")
  print(table(data$Spawning.Component[sc_ok]))
  cat("\nSpawning component by assigned stock:\n")
  print(table(data$Spawning.Component[sc_ok], data$Assigned.Stock[sc_ok]))
}, file = file.path(RESULTS, "SpawningComponent", "summary.txt"))
cat("\n")

stopCluster(cl); registerDoSEQ()

# ============================================================
# 3. Summary table
# ============================================================
t3 <- bind_rows(lapply(res, function(r) data.frame(
  Analysis = r$label, n = r$n, k = r$k,
  Accuracy = 100 * r$acc, AccuracySD = 100 * r$accSD,
  Kappa = r$kap, KappaSD = r$kapSD)))
write.csv(t3, file.path(RESULTS, "table3_remaining.csv"), row.names = FALSE)
cat("\n--- Spawning-group model (n =", nrow(data), ") ---\n")
print(t3, digits = 4)
cat("\nDone. Everything in:", RESULTS, "\n")
