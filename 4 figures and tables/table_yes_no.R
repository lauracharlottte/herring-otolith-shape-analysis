## ============================================================
## Table 5 -- LDA output reduced to a yes/no decision
## ------------------------------------------------------------
## Does the otolith belong to the focal population (WBSS or BAS)
## or not? For each approach, from the cross-validated confusion
## matrix that main_analysis.R saved (caret::confusionMatrix(fit),
## averaged over 400 runs, cells in per cent of all otoliths;
## rows = assigned, columns = true population):
##
##   true positives (%)  diagonal cell of the focal population divided
##                       by its column total (= correctly classified,
##                       Table 4)
##   true negatives (%)  all cells outside the focal row and column,
##                       divided by the column totals of the other groups
##   accuracy (%)        diagonal cell + cells outside the focal row and
##                       column, i.e. the share of correct yes/no decisions
##
## Writes  results final plots/tables/table5_yes_no.csv
## ============================================================

suppressPackageStartupMessages({ library(dplyr); library(caret) })

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

yes_no <- function(folder, label, focal) {
  e <- new.env(); load(file.path(RESULTS, folder, "LDA_results.RData"), envir = e)
  fit <- e$fit
  tb  <- confusionMatrix(fit)$table            # per cent of all otoliths
  i   <- match(focal, colnames(tb))
  cs  <- colSums(tb)
  y   <- fit$trainingData$.outcome
  data.frame(
    Focal          = focal,
    Analysis       = label,
    n_total        = length(y),
    n_focal        = sum(y == focal),
    TruePositives  = 100 * tb[i, i] / cs[i],
    TrueNegatives  = 100 * sum(tb[-i, -i]) / sum(cs[-i]),
    Accuracy       = tb[i, i] + sum(tb[-i, -i]),
    row.names = NULL)
}

approaches <- c("00_full" = "Full data set", "age0_3" = "Age 0-3",
                "age1_3" = "Age 1-3", "age1_10" = "Age 1-10")

tab <- bind_rows(
  bind_rows(lapply(names(approaches), function(f) yes_no(f, approaches[[f]], "WBSS"))),
  yes_no("WBSS_vs_Rest", "WBSS vs. rest", "WBSS"),
  bind_rows(lapply(names(approaches), function(f) yes_no(f, approaches[[f]], "BAS"))),
  yes_no("BAS_vs_Rest", "BAS vs. rest", "BAS"))

out <- tab %>% mutate(across(c(TruePositives, TrueNegatives, Accuracy), ~ round(.x, 1)))
f <- file.path(ensure_dir(file.path(DIR_FINAL, "tables")), "table5_yes_no.csv")
write.csv(out, f, row.names = FALSE)
print(out, row.names = FALSE)
cat("\nwritten:", f, "\n")
