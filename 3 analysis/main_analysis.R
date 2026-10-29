# ============================================================
# Main analysis: outlines -> shape descriptors -> classification
# ------------------------------------------------------------
# Runs the whole evaluation in one pass and writes everything
# into one folder:
#
#     Fourier and Classification/Results/
#
# Requires oto_all.RData (the extracted outlines, see
# 1 outline extraction/ex_oto.R) to be present in that folder.
#
# NEATL (n = 2) and NSS (n = 3) are excluded from EVERY analysis,
# see EXCLUDE_STOCKS below.
#
# Produces:
#   00_full/          the five stocks, main model
#   age0_3/ age1_3/ age1_10/
#   WBSS_vs_Rest/ BAS_vs_Rest/
#   Imbalance/        priors, down-/upsampling
#   table1_performance.csv   <- numbers for LaTeX table 1
#   table2_imbalance.csv     <- numbers for LaTeX table 2
#   summary.txt
#
# The thesis figures are produced by the scripts in
# 4 figures and tables/.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
  library(Momocs); library(caret); library(MASS); library(doParallel)
})
select <- dplyr::select

# -----------------------------
# Settings
# -----------------------------
NB_H        <- 12          # harmonics, as in the thesis
CV_FOLDS    <- 4
CV_REPEATS  <- 100
SEED        <- 2021
## NEATL (n = 2) and NSS (n = 3) are left out of every analysis.
## With two and three fish these classes cannot be estimated in
## 4-fold cross-validation -- most folds hold none of them -- they
## are never classified correctly, and in a 7-class LDA their two-
## and three-fish means can claim whole discriminant axes. They are
## removed before outlier detection, so the five
## stocks define the study population from the first step on.
EXCLUDE_STOCKS <- c("NEATL", "NSS")

## ---- load project paths ------------------------------------
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

RESULTS <- ensure_dir(DIR_RESULTS)
stopifnot(file.exists(file.path(RESULTS, "oto_all.RData")))

# ============================================================
# 1. Outlines + metadata
# ============================================================
load(file.path(RESULTS, "oto_all.RData"))       # oto
cat("Outlines loaded:", length(oto), "\n")

id_photo <- tibble(
  ID              = names(oto),
  IdentNr.Picture = as.integer(sub(".*__(\\d+)_.*", "\\1", names(oto))),
  Year            = as.integer(sub(".*_(\\d{4})_SD.*", "\\1", names(oto))),
  SD              = as.integer(sub(".*_SD(\\d+).*", "\\1", names(oto))),
  Side            = sub(".*_", "", names(oto)))

data <- read.csv(F_DATA_AGE, row.names = NULL, fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)
valid.stocks <- setdiff(c("BAS","CBH","Downs","NEATL","NSAS","NSS","WBSS"),
                        EXCLUDE_STOCKS)
data <- data[data$Assigned.Stock %in% valid.stocks, ]
data <- data[!is.na(data$IdentNr.Picture), ]
data <- data[!is.na(data$Spawning.Component), ]
data$Subdivision <- as.integer(data$Subdivision)
data$Year        <- as.integer(data$Year)

data <- data %>%
  left_join(id_photo, by = c("IdentNr.Picture" = "IdentNr.Picture",
                             "Subdivision" = "SD", "Year" = "Year")) %>%
  filter(ID %in% names(oto)) %>%
  arrange(match(ID, names(oto)))
oto <- oto[names(oto) %in% data$ID]
stopifnot(nrow(data) == length(oto))
cat("after metadata join:", nrow(data), "\n")

slide_to_bottom <- function(coo) coo_slide(coo, id = which.min(coo[, 2]))
build_out <- function(oto, data) {
  O <- Out(oto, fac = data); O$coo <- lapply(O$coo, slide_to_bottom)
  O <- Out(O$coo, fac = data); O$fac <- droplevels(O$fac); O
}
OTO <- build_out(oto, data)

# ============================================================
# 2. Outliers (Mahalanobis, first 5 PCs, p < 0.001)
# ============================================================
eOTO_prelim <- efourier(OTO, nb.h = NB_H, smooth.it = 0, norm = TRUE)
pca_all     <- PCA(eOTO_prelim)
sc          <- as.data.frame(pca_all$x[, 1:5])
md          <- mahalanobis(sc, colMeans(sc), cov(sc))
outliers    <- rownames(pca_all$x)[md > qchisq(0.999, df = 5)]
cat("outliers:", length(outliers), "\n")
write.csv(data.frame(ID = outliers), file.path(RESULTS, "outliers.csv"), row.names = FALSE)

oto  <- oto[!(names(oto) %in% outliers)]
data <- data %>% filter(ID %in% names(oto)) %>% arrange(match(ID, names(oto)))
OTO  <- build_out(oto, data)
cat("analysed:", nrow(data), "\n\n")

# ============================================================
# 3. EFD + shape descriptors
# ============================================================
set.seed(SEED)
harmpow <- calibrate_harmonicpower_efourier(x = OTO, nb.h = 30)
cat("harmonics for 99 % / 99.9 % shape power:",
    harmpow$minh[3], "/", harmpow$minh[4], " -- using:", NB_H, "\n")

eOTO      <- efourier(OTO, nb.h = NB_H, smooth.it = 0, norm = TRUE)
dataSHAPE <- as.data.frame(eOTO$coe) %>% select(-A1, -B1, -C1)
stopifnot(ncol(dataSHAPE) == NB_H * 4 - 3)
data$Age <- suppressWarnings(as.numeric(as.character(data$Age)))

save(data,      file = file.path(RESULTS, "dataBIO.RData"))
save(dataSHAPE, file = file.path(RESULTS, "dataSHAPE.RData"))
save(eOTO,      file = file.path(RESULTS, "eOTO.RData"))
save(OTO,       file = file.path(RESULTS, "OTO.RData"))

# ============================================================
# 4. Helpers
# ============================================================
cl <- makeCluster(detectCores()); registerDoParallel(cl)
trcntr <- trainControl(method = "repeatedcv", number = CV_FOLDS,
                       repeats = CV_REPEATS, allowParallel = TRUE)

## pull the numbers out of a caret confusion matrix
metrics <- function(fit) {
  tb <- confusionMatrix(fit)$table
  dg <- diag(tb); cs <- colSums(tb); rs <- rowSums(tb)
  list(acc = mean(fit$resample$Accuracy), accSD = sd(fit$resample$Accuracy),
       kap = mean(fit$resample$Kappa),    kapSD = sd(fit$resample$Kappa),
       sens = 100 * dg / cs, fdr = 100 * (1 - dg / rs), tb = tb)
}

## one complete sub-analysis: cross-validated LDA + output
run_analysis <- function(label, folder, idx, group_col = "Assigned.Stock",
                         group_vec = NULL) {

  dir <- ensure_dir(file.path(RESULTS, folder))
  d   <- data[idx, , drop = FALSE]
  sh  <- dataSHAPE[idx, , drop = FALSE]
  grp <- if (is.null(group_vec)) droplevels(factor(d[[group_col]]))
         else droplevels(factor(group_vec[idx]))

  cat("\n=====", label, "-- n =", nrow(d), "| classes:",
      paste(levels(grp), collapse = ", "), "\n")

  dcl <- cbind(Group = grp, sh)

  ## cross-validated LDA
  set.seed(SEED)
  fit <- train(Group ~ ., data = dcl, method = "lda",
               preProcess = c("scale", "center"), trControl = trcntr)
  m <- metrics(fit)

  save(fit, file = file.path(dir, "LDA_results.RData"))
  sink(file.path(dir, "LDA_results.txt"))
  cat("===", label, "===\n\n"); print(fit)
  cat("\n=== Confusion Matrix ===\n"); print(confusionMatrix(fit))
  cat("\nClass sizes:\n"); print(table(grp))
  cat("\nCorrectly classified / sensitivity (%):\n"); print(round(m$sens, 1))
  cat("\nFalsely assigned / false discovery rate (%):\n"); print(round(m$fdr, 1))
  sink()

  cat("  Accuracy", round(m$acc, 4), "| Kappa", round(m$kap, 4), "\n")
  c(m, list(label = label, n = nrow(d), k = nlevels(grp)))
}

# ============================================================
# 5. The six analyses for table 1
# ============================================================
res <- list()
res$full   <- run_analysis("Full data set", "00_full", seq_len(nrow(data)))
res$age0_3 <- run_analysis("Age 0-3", "age0_3",
                           which(!is.na(data$Age) & data$Age >= 0 & data$Age <= 3))
res$age1_3 <- run_analysis("Age 1-3", "age1_3",
                           which(!is.na(data$Age) & data$Age >= 1 & data$Age <= 3))
res$age1_10<- run_analysis("Age 1-10", "age1_10",
                           which(!is.na(data$Age) & data$Age >= 1 & data$Age <= 10))

g_wbss <- factor(ifelse(data$Assigned.Stock == "WBSS", "WBSS", "Rest"),
                 levels = c("WBSS", "Rest"))
res$wbss <- run_analysis("WBSS vs. rest", "WBSS_vs_Rest", seq_len(nrow(data)),
                         group_vec = g_wbss)

g_bas <- factor(ifelse(data$Assigned.Stock == "BAS", "BAS", "Rest"),
                levels = c("BAS", "Rest"))
res$bas <- run_analysis("BAS vs. rest", "BAS_vs_Rest", seq_len(nrow(data)),
                        group_vec = g_bas)

## assemble table 1
pick <- function(v, nm) if (nm %in% names(v)) unname(v[nm]) else NA_real_
t1 <- bind_rows(lapply(res, function(r) data.frame(
  Analysis  = r$label, n = r$n, k = r$k,
  Accuracy  = 100 * r$acc, AccuracySD = 100 * r$accSD,
  Kappa     = r$kap,       KappaSD    = r$kapSD,
  Sens_WBSS = pick(r$sens, "WBSS"), Sens_BAS = pick(r$sens, "BAS"),
  FDR_WBSS  = pick(r$fdr,  "WBSS"), FDR_BAS  = pick(r$fdr,  "BAS"))))
write.csv(t1, file.path(RESULTS, "table1_performance.csv"), row.names = FALSE)
cat("\n--- Classification performance (thesis Tab. 4) ---\n"); print(t1, digits = 4)

# ============================================================
# 6. Class-imbalance analysis (thesis Tab. 3)
# ============================================================
imb.dir <- ensure_dir(file.path(RESULTS, "Imbalance"))
## The same five stocks as every other analysis, now that NEATL and
## NSS go up front. "Observed priors" below is therefore the full
## model of table 1 and has to reproduce its numbers exactly.
dimb    <- cbind(Assigned.Stock = droplevels(factor(data$Assigned.Stock)),
                 dataSHAPE)
STOCKS  <- levels(dimb$Assigned.Stock); K <- length(STOCKS)
cat("\nImbalance data set: n =", nrow(dimb), "| classes:",
    paste(STOCKS, collapse = ", "), "\n")
print(table(dimb$Assigned.Stock))

PRIOR_EQUAL <- rep(1/K, K)
clusterExport(cl, "PRIOR_EQUAL", envir = environment())

fit_imb <- function(sampling = NULL, prior = NULL) {
  ctl <- trainControl(method = "repeatedcv", number = CV_FOLDS,
                      repeats = CV_REPEATS, sampling = sampling,
                      allowParallel = TRUE)
  a <- list(Assigned.Stock ~ ., data = dimb, method = "lda",
            preProcess = c("scale","center"), trControl = ctl)
  if (!is.null(prior)) a$prior <- prior
  set.seed(SEED); do.call(train, a)
}

variants <- list(
  "Observed priors"   = list(NULL,   NULL),
  "Equal priors"      = list(NULL,   PRIOR_EQUAL),
  "Downsampling"      = list("down", NULL),
  "Upsampling"        = list("up",   NULL))

t2 <- list()
for (nm in names(variants)) {
  cat("\n-- Imbalance:", nm, "--\n")
  f <- fit_imb(variants[[nm]][[1]], variants[[nm]][[2]])
  m <- metrics(f)
  save(f, file = file.path(imb.dir, paste0("fit_", gsub("[^A-Za-z0-9]", "_", nm), ".RData")))
  sink(file.path(imb.dir, paste0("confusion_", gsub("[^A-Za-z0-9]", "_", nm), ".txt")))
  cat("===", nm, "===\n\n"); print(confusionMatrix(f))
  cat("\nCorrectly classified / sensitivity (%):\n"); print(round(m$sens, 1))
  cat("\nFalsely assigned / false discovery rate (%):\n"); print(round(m$fdr, 1))
  cat("\nN in final model:", f$finalModel$N, "\n")
  cat("Priors used:\n"); print(round(f$finalModel$prior, 3))
  sink()
  t2[[nm]] <- data.frame(
    Model = nm, N = f$finalModel$N,
    Accuracy = 100*m$acc, AccuracySD = 100*m$accSD,
    Kappa = m$kap, KappaSD = m$kapSD,
    t(setNames(m$sens, paste0("Sens_", names(m$sens)))),
    t(setNames(m$fdr,  paste0("FDR_",  names(m$fdr)))),
    check.names = FALSE, row.names = NULL)
  cat("   Accuracy", round(m$acc,4), "| Kappa", round(m$kap,4), "| N", f$finalModel$N, "\n")
}
t2 <- bind_rows(t2)
write.csv(t2, file.path(RESULTS, "table2_imbalance.csv"), row.names = FALSE)
cat("\n--- Table 2 ---\n"); print(t2, digits = 4)

stopCluster(cl); registerDoSEQ()

# ============================================================
# 7. summary.txt
# ============================================================
sink(file.path(RESULTS, "summary.txt"))
cat("Otolith shape analysis -- main analysis\n")
cat("Excluded from all analyses:", paste(EXCLUDE_STOCKS, collapse = ", "),
    "(too few fish to estimate)\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Outlines total    :", nrow(data) + length(outliers), "\n")
cat("Outliers removed  :", length(outliers), "\n")
cat("Analysed          :", nrow(data), "\n")
cat("Harmonics         :", NB_H, "\n")
cat("Cross-validation  : ", CV_FOLDS, "-fold x ", CV_REPEATS, " repeats\n", sep = "")
cat("\nStocks:\n"); print(table(data$Assigned.Stock))
cat("\nYear x subdivision:\n"); print(table(data$Year, data$Subdivision))
cat("\n\n=== Table 1: LDA performance ===\n"); print(t1, digits = 4)
cat("\n\n=== Table 2: class-imbalance treatment ===\n"); print(t2, digits = 4)
cat("\n\nColumn meanings\n")
cat("  n / N        number of otoliths\n")
cat("  k            number of classes\n")
cat("  Accuracy     mean over ", CV_FOLDS * CV_REPEATS,
    " cross-validation runs, in per cent\n", sep = "")
cat("  AccuracySD   standard deviation across those runs\n")
cat("  Kappa        accuracy corrected for chance agreement\n")
cat("  Sens_X       correctly classified: of the otoliths that ARE stock X,\n")
cat("               the percentage the model found\n")
cat("  FDR_X        falsely assigned: of the otoliths LABELLED stock X,\n")
cat("               the percentage that belonged to another stock\n")
sink()

cat("\nDone. Everything in:\n ", RESULTS, "\n")
