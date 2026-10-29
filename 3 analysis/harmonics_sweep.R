# ============================================================
# Choice of the number of harmonics
# ------------------------------------------------------------
# Sweep over 2-35 harmonics -> Results/Harmonics_sweep/
# Reads the main-model sample from Results/, so it matches
# the main model exactly.
#
# A note on features: the
# sweep uses PC scores up to 99 % cumulative variance, not the raw
# EFD descriptors of the main model. At high h the raw
# coefficients make the within-class covariance singular for the
# small stocks inside a fold and every resample fails silently.
# The sweep therefore shows the SHAPE of the curve, not numbers
# comparable to the main model.
#
# The figure (Fig. 12b, c) is drawn by
# 4 figures and tables/fig_number_of_harmonics.R from the saved CSV.
#
# Runtime: about ten minutes.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Momocs); library(caret); library(MASS); library(doParallel)
})
select <- dplyr::select

SWEEP     <- 2:35          # B
CV_FOLDS   <- 4
CV_REPEATS <- 100
SEED       <- 2021

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
RESULTS <- DIR_RESULTS

load(file.path(RESULTS, "dataBIO.RData"))
load(file.path(RESULTS, "oto_all.RData"))
oto  <- oto[names(oto) %in% data$ID]
data <- data %>% arrange(match(ID, names(oto)))
stopifnot(nrow(data) == length(oto))
cat("sample:", nrow(data), "otoliths\n")

slide_to_bottom <- function(coo) coo_slide(coo, id = which.min(coo[, 2]))
OTO <- Out(oto, fac = data)
OTO$coo <- lapply(OTO$coo, slide_to_bottom)
OTO <- Out(OTO$coo, fac = data)
OTO$fac <- droplevels(OTO$fac)
stock <- droplevels(factor(OTO$fac$Assigned.Stock))

cl <- makeCluster(detectCores()); registerDoParallel(cl)

# ============================================================
# B. Sweep over 2-35 harmonics
# ============================================================
cat("\n=== B. sweep", min(SWEEP), "-", max(SWEEP), "harmonics ===\n")
out.B <- ensure_dir(file.path(RESULTS, "Harmonics_sweep"))
ctrl_s <- trainControl(method = "repeatedcv", number = CV_FOLDS,
                       repeats = CV_REPEATS, allowParallel = TRUE)

sweep <- list()
for (h in SWEEP) {
  r <- tryCatch({
    eO <- efourier(OTO, nb.h = h, smooth.it = 0, norm = TRUE)
    pc <- PCA(eO)
    ve <- cumsum(pc$eig) / sum(pc$eig)
    nc <- which(ve >= 0.99)[1]; if (is.na(nc)) nc <- ncol(pc$x)
    nc <- max(nc, 2)
    sc <- as.data.frame(pc$x[, 1:nc, drop = FALSE]); sc$stock <- stock
    set.seed(SEED)
    f <- train(stock ~ ., data = sc, method = "lda", trControl = ctrl_s)
    if (all(is.na(f$results$Accuracy))) stop("every resample failed")
    data.frame(Harmonics = h, nComponents = nc,
               Accuracy = f$results$Accuracy, Kappa = f$results$Kappa)
  }, error = function(e) {
    message("  SKIPPED h = ", h, " -- ", conditionMessage(e)); NULL })
  if (!is.null(r)) {
    sweep[[as.character(h)]] <- r
    cat("  h =", sprintf("%2d", h), "| PCs", sprintf("%2d", r$nComponents),
        "| acc", round(r$Accuracy, 4), "| kappa", round(r$Kappa, 4), "\n")
  }
}
sweep <- bind_rows(sweep)
write.csv(sweep, file.path(out.B, "harmonics_accuracy_kappa.csv"), row.names = FALSE)
save(sweep, file = file.path(out.B, "harmonics_accuracy_kappa.RData"))

stopCluster(cl); registerDoSEQ()

sink(file.path(out.B, "summary.txt"))
cat("Harmonics sweep on the analysed sample\n")
cat("Run:", format(Sys.time()), "\n\n")
cat("Specimens :", nrow(data), "\n")
cat("Features  : PC scores up to 99 % cumulative variance\n")
cat("CV        : ", CV_FOLDS, "-fold x ", CV_REPEATS, " repeats\n\n", sep = "")
print(sweep, digits = 5, row.names = FALSE)
cat("\nBest accuracy at h =", sweep$Harmonics[which.max(sweep$Accuracy)],
    "| best kappa at h =", sweep$Harmonics[which.max(sweep$Kappa)], "\n")
cat("Accuracy at h = 12:", round(sweep$Accuracy[sweep$Harmonics == 12], 4), "\n")
sink()

cat("\nDone.\n  ", out.B, "\n")
