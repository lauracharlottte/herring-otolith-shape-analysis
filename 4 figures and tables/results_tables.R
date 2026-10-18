## ============================================================
## Results tables 1-3 (LaTeX), generated from the result files
## ------------------------------------------------------------
## Writes ONE plain-text file with the LaTeX code of all three tables of
## the results chapter, in chapter order:
##   results final plots/tables/Tables results chapter (LaTeX code).txt
##   Table 1  population composition of all genetically assigned herring
##            (written by stock_composition.R,
##            copied in here)
##   Table 2  LDA performance per analysis
##   Table 3  class-imbalance treatments
## Plain .txt so the blocks can be copied into the thesis directly.
##
## Every number comes straight from Results/table{1,2,3}_*.csv.
## Values are rounded once, from full precision: percentages to one decimal,
## kappa to two.
##
## Units: accuracy and the per-population columns are percentages of
## otoliths; Cohen's kappa is a dimensionless coefficient (-1 to 1) and
## therefore has no unit -- it must not be written as a percentage.
##
## Needs \usepackage{booktabs} and \usepackage{adjustbox} in the thesis
## preamble. adjustbox shrinks a table only if it is wider than the text
## (Table 3 at a 12 pt base font); otherwise it changes nothing.
## ============================================================

suppressPackageStartupMessages(library(dplyr))

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

RESULTS   <- DIR_RESULTS
out.dir <- ensure_dir(file.path(DIR_FINAL, "tables"))

t1 <- read.csv(file.path(RESULTS, "table1_performance.csv"), stringsAsFactors = FALSE)
t2 <- read.csv(file.path(RESULTS, "table2_imbalance.csv"),   stringsAsFactors = FALSE)
t3 <- read.csv(file.path(RESULTS, "table3_remaining.csv"),   stringsAsFactors = FALSE)

## The spawning-group model is reported as the last row of the
## classification table.
t1 <- bind_rows(t1, t3 %>% filter(Analysis == "Spawning component"))

## group sizes for the table notes, from the analysed sample
load(file.path(RESULTS, "dataBIO.RData"))                  # data
n_pop   <- table(data$Assigned.Stock)
n_small <- min(n_pop); n_large <- max(n_pop)             # CBH, WBSS
sp      <- table(data$Spawning.Component, data$Assigned.Stock)
sp_txt  <- function(comp) {
  pops <- intersect(POP_ORDER, colnames(sp)[sp[comp, ] > 0])
  paste0(sub(" spawner", " spawners", comp), " (",
         paste(pops, collapse = " and "), ", $n$ = ", formatC(sum(sp[comp, ]), format = "d", big.mark = "{,}"), ")")
}
SPAWN_NOTE <- paste0("Spawning component: three groups -- ",
  paste(vapply(c("autumn spawner", "spring spawner", "winter spawner"), sp_txt, ""),
        collapse = ", "), ".")

## ---- formatting ----------------------------------------------
pct <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")   # 3029 -> 3,029, as in the text
pm  <- function(m, s, d) sprintf(paste0("%.", d, "f $\\pm$ %.", d, "f"), m, s)
lab <- function(x) gsub("vs\\. ", "vs.\\\\ ",                      # vs. rest -> vs.\ rest
                         gsub("(\\d)-(\\d)", "\\1--\\2", x))     # 0-3 -> 0--3 (en dash)
row <- function(...) paste0("    ", paste(c(...), collapse = " & "), " \\\\")

note <- function(txt) c(
  "  \\vspace{0.4em}",
  "  \\begin{minipage}{\\linewidth}",
  "    \\footnotesize",
  paste0("    \\textit{Note.} ", txt),
  "  \\end{minipage}")

CC_FA <- paste0(
  "\\emph{Correctly classified}: percentage of the otoliths of a population ",
  "that were assigned to that population. \\emph{Falsely assigned}: percentage ",
  "of the otoliths assigned to a population that actually belonged to another ",
  "population, i.e.\\ the false positives among all assignments to that population.")

CV_NOTE <- paste0(
  "Accuracy and $\\kappa$ are means $\\pm$ standard deviation over 400 ",
  "cross-validation runs (4-fold cross-validation, 100 repeats). Accuracy is given ",
  "in per cent; Cohen's $\\kappa$ is a dimensionless coefficient that corrects the ",
  "accuracy for agreement expected by chance (0 = chance level, 1 = perfect ",
  "agreement).")
CM_NOTE <- paste0("Correctly classified and falsely assigned were calculated from the ",
  "cross-validated confusion matrix, averaged over all runs. The two measures ",
  "refer to different totals (the otoliths of a population and the otoliths ",
  "assigned to it) and therefore do not add up to 100\\%.")

## ============================================================
## Table 1 -- LDA performance per analysis (thesis: Table 2)
## ============================================================
## Correctly classified / falsely assigned for ALL five populations,
## computed from the saved models exactly as main_analysis.R does
## (metrics(): averaged cross-validated confusion matrix; correctly
## classified = diagonal / column total, falsely assigned = 1 - diagonal /
## row total). table1_performance.csv only holds WBSS and BAS -- those two
## are checked against it, so both routes give the same numbers.
## Layout: population names once as column heads, correctly classified in
## the upper block, falsely assigned in the lower one -- 10 columns instead
## of 14, so the table fits the page at normal table size.
pops_hdr <- paste(POP_ORDER, collapse = " & ")
FOLDER <- c("Full data set" = "00_full", "Age 0-3" = "age0_3", "Age 1-3" = "age1_3",
            "Age 1-10" = "age1_10", "WBSS vs. rest" = "WBSS_vs_Rest",
            "BAS vs. rest" = "BAS_vs_Rest", "Spawning component" = "SpawningComponent")
stopifnot(all(t1$Analysis %in% names(FOLDER)))
pop_metrics <- function(folder) {
  e <- new.env(); load(file.path(RESULTS, folder, "LDA_results.RData"), envir = e)
  tb <- caret::confusionMatrix(e$fit)$table
  dg <- diag(tb)
  list(sens = 100 * dg / colSums(tb), fdr = 100 * (1 - dg / rowSums(tb)))
}
mets <- lapply(FOLDER[t1$Analysis], pop_metrics)
## NA where a population is not a group of that model ("Rest", spawning)
cc <- t(vapply(mets, function(m) unname(m$sens[POP_ORDER]), numeric(5)))
fa <- t(vapply(mets, function(m) unname(m$fdr[POP_ORDER]),  numeric(5)))
colnames(cc) <- colnames(fa) <- POP_ORDER
same <- function(a, b) all(abs(a - b) < 1e-6 | (is.na(a) & is.na(b)))
stopifnot(same(cc[, "WBSS"], t1$Sens_WBSS), same(cc[, "BAS"], t1$Sens_BAS),
          same(fa[, "WBSS"], t1$FDR_WBSS),  same(fa[, "BAS"], t1$FDR_BAS))

block_label <- function(n_before, txt)
  row(rep("", n_before), paste0("\\multicolumn{5}{c}{\\textit{", txt, " (\\%)}}"))

body1 <- c(
  block_label(4, "Correctly classified"),
  vapply(seq_len(nrow(t1)), function(i) with(t1[i, ], row(
    lab(Analysis), fmt_n(n), k, pm(Accuracy, AccuracySD, 1), pm(Kappa, KappaSD, 2),
    pct(cc[i, ]))), ""),
  "    \\midrule",
  block_label(4, "Falsely assigned"),
  vapply(which(rowSums(!is.na(fa)) > 0), function(i) row(
    lab(t1$Analysis[i]), "", "", "", "", pct(fa[i, ])), ""))

tab1 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0("  \\caption[Classification performance of the LDA per analysis]{",
         "Classification performance of the linear discriminant analysis (LDA) for ",
         "each analysis. ", CC_FA, "}"),
  "  \\label{tab:classification}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{4pt}",
  "  \\begin{adjustbox}{max width=\\linewidth}",
  "  \\begin{tabular}{lrrccccccc}",
  "    \\toprule",
  paste0("    Analysis & $n$ & $k$ & Accuracy (\\%) & $\\kappa$ & ", pops_hdr, " \\\\"),
  "    \\midrule",
  body1,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\end{adjustbox}",
  note(paste0(CV_NOTE, " ", CM_NOTE, " $n$ = number of otoliths, $k$ = number of groups. In the ",
              "two-group models only the focal population is given; dashes mark ",
              "values that do not apply. ", SPAWN_NOTE)),
  "\\end{table}")

## ============================================================
## Table 2 -- class-imbalance treatments (full data set; thesis: Table 3)
## ============================================================
sens2 <- paste0("Sens_", POP_ORDER); fdr2 <- paste0("FDR_", POP_ORDER)
stopifnot(all(c(sens2, fdr2) %in% names(t2)))
body2 <- c(
  block_label(3, "Correctly classified"),
  vapply(seq_len(nrow(t2)), function(i) row(
    t2$Model[i], fmt_n(t2$N[i]),
    pm(t2$Accuracy[i], t2$AccuracySD[i], 1), pm(t2$Kappa[i], t2$KappaSD[i], 2),
    pct(unlist(t2[i, sens2]))), ""),
  "    \\midrule",
  block_label(3, "Falsely assigned"),
  vapply(seq_len(nrow(t2)), function(i) row(
    t2$Model[i], "", "", "", pct(unlist(t2[i, fdr2]))), ""))

tab2 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0("  \\caption[Effect of the class-imbalance treatments on the LDA]{",
         "Effect of different treatments of class imbalance on the LDA classification ",
         "of the five populations (full data set, $n$ = ", formatC(nrow(data), format = "d", big.mark = "{,}"), "). ", CC_FA, "}"),
  "  \\label{tab:imbalance}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{4pt}",
  "  \\begin{adjustbox}{max width=\\linewidth}",
  "  \\begin{tabular}{lrccccccc}",
  "    \\toprule",
  paste0("    Treatment & $N$ & Accuracy (\\%) & $\\kappa$ & ", pops_hdr, " \\\\"),
  "    \\midrule",
  body2,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\end{adjustbox}",
  note(paste0(
    "Observed priors: LDA prior probabilities equal to the proportions of the ",
    "populations in the sample -- this model is identical to the full data set in ",
    "Table~\\ref{tab:classification}. Equal priors: 0.2 for each population. ",
    "Down- and upsampling: within each training fold, every population was reduced ",
    "to the size of the smallest or increased to the size of the largest population; ",
    "$N$ is the number of otoliths in the resampled final model (5 $\\times$ ", n_small,
    " and 5 $\\times$ ", formatC(n_large, format = "d", big.mark = "{,}"), "). The held-out folds used for the evaluation were ",
    "not resampled. ", CV_NOTE, " ", CM_NOTE)),
  "\\end{table}")

## ---- Table 1: copied in from stock_composition.R ---------------
gen_tex <- file.path(DIR_FOURIER, "Composition", "table_composition_genetic_total.tex")
if (!file.exists(gen_tex))
  stop('run stock_composition.R first', call. = FALSE)
if (file.mtime(gen_tex) < file.mtime(file.path(RESULTS, "dataBIO.RData")))
  warning("table_composition_genetic_total.tex is older than dataBIO.RData -- ",
          'run stock_composition.R first', call. = FALSE)
tab0 <- readLines(gen_tex)

## ---- write one plain-text file ---------------------------------
banner <- function(title) c("", strrep("=", 72), title, strrep("=", 72), "")
txt <- c(
  "TABLES FOR THE RESULTS CHAPTER -- LaTeX code",
  "",
  "Copy each block into the thesis. Needs \\usepackage{booktabs} and \\usepackage{adjustbox} in the preamble.",
  "Generated from the result files by Master-Thesis/4 figures and tables/results_tables.R",
  "(Table 1 by stock_composition.R). Do not change numbers by hand -- run the scripts again.",
  "Table numbers = order in the Results chapter; LaTeX numbers them automatically.",
  banner("TABLE 1 -- Population composition of all genetically assigned herring (Section 3.1)"), tab0,
  banner("TABLE 2 -- Classification performance per analysis (Section 3.7)"), tab1,
  banner("TABLE 3 -- Class-imbalance treatments (Section 3.7)"), tab2)
f <- file.path(out.dir, "Tables results chapter (LaTeX code).txt")
writeLines(txt, f)
cat("written:", f, "\n")
