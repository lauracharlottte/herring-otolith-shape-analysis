# ============================================================
# Population composition by year and ICES subdivision
# ------------------------------------------------------------
# Table 2 of the thesis (all genetically assigned herring) and the
# numbers in Section 3.1 on how many fish of each population entered
# the shape analysis.
#
# Percentages are computed WITHIN a year x subdivision cell,
# never within a year alone. Sampling effort moved between
# subdivisions across the three surveys, and SD 23 is WBSS-rich,
# so a percentage pooled over subdivisions mostly reports where
# the survey sampled, not how the populations changed.
#
# The analysed shape sample is compared with the full genetically
# assigned set, to show that dropping fish without a usable
# photograph did not bias the composition.
#
# Writes to: Fourier and Classification/Composition/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

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
out.dir <- ensure_dir(file.path(DIR_FOURIER, "Composition"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
VALID <- c("BAS", "CBH", "Downs", "NEATL", "NSAS", "NSS", "WBSS")

## Column order = POP_ORDER from config.R: WBSS, BAS, NSAS, Downs, CBH.
STACK <- POP_ORDER

## NEATL (n = 2) and NSS (n = 3) are excluded from the shape
## analysis itself (main_analysis.R), so the analysed sample
## holds none of them. They remain only in the genetic-total table.
TINY <- c("NEATL", "NSS")

# ------------------------------------------------------------
# 1. The two data sets
# ------------------------------------------------------------
## A -- the otoliths that actually entered the shape analysis
load(file.path(RESULTS, "dataBIO.RData"))          # data
A <- data %>%
  transmute(Year  = as.integer(Year),
            SD    = as.integer(Subdivision),
            Stock = Assigned.Stock)

## B -- every fish with a clean genetic stock assignment,
##      whether or not a usable right otolith photo exists
B <- read.csv(F_DATA_AGE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE) %>%
  filter(Assigned.Stock %in% VALID, !is.na(Subdivision), !is.na(Year)) %>%
  transmute(Year  = as.integer(Year),
            SD    = as.integer(Subdivision),
            Stock = Assigned.Stock)

cat("analysed shape sample :", nrow(A), "\n")
cat("genetically assigned  :", nrow(B), "\n")
cat("coverage              :", sprintf("%.1f %%", 100 * nrow(A) / nrow(B)), "\n\n")

# ------------------------------------------------------------
# 2. Long table: counts and within-cell percentages
# ------------------------------------------------------------
compose <- function(df, set_label) {
  df %>%
    mutate(Stock = factor(Stock, levels = c(STACK, TINY))) %>%
    count(Year, SD, Stock, .drop = FALSE) %>%
    group_by(Year, SD) %>%
    mutate(cell_n  = sum(n),
           percent = round(100 * n / cell_n, 1)) %>%
    ungroup() %>%
    mutate(Set = set_label, .before = 1)
}

long <- bind_rows(compose(A, "shape sample"),
                  compose(B, "genetic total"))

write.csv(long, file.path(out.dir, "stock_composition_long.csv"), row.names = FALSE)

## Does dropping un-photographed fish shift the composition?
## (a) within each year x subdivision cell
bias <- long %>%
  filter(Stock %in% STACK) %>%
  select(Set, Year, SD, Stock, percent) %>%
  pivot_wider(names_from = Set, values_from = percent) %>%
  mutate(diff_pp = round(`shape sample` - `genetic total`, 1))
cat("largest shift within a year x subdivision cell:",
    sprintf("%.1f pp\n", max(abs(bias$diff_pp), na.rm = TRUE)))
print(bias %>% slice_max(abs(diff_pp), n = 1, with_ties = FALSE) %>% as.data.frame(), row.names = FALSE)
write.csv(bias, file.path(out.dir, "coverage_bias_check.csv"), row.names = FALSE)

## (b) per population, pooled over all cells: share of the fish that
## entered the analyses, and shift of the population proportion
retention <- full_join(count(A, Stock, name = "analysed"),
                       count(B %>% filter(Stock %in% STACK), Stock, name = "assigned"),
                       by = "Stock") %>%
  mutate(retained_pct = round(100 * analysed / assigned, 1),
         share_analysed_pct = 100 * analysed / sum(analysed),
         share_assigned_pct = 100 * assigned / sum(assigned),
         shift_pp = round(share_analysed_pct - share_assigned_pct, 1)) %>%
  mutate(across(starts_with("share"), ~ round(.x, 1))) %>%
  arrange(match(Stock, STACK))
cat("\nper population (Section 3.1):\n")
print(as.data.frame(retention), row.names = FALSE)
cat("largest shift of a population proportion:",
    sprintf("%.1f pp\n", max(abs(retention$shift_pp))))
write.csv(retention, file.path(out.dir, "retention_by_population.csv"), row.names = FALSE)

cat("\nWritten to:", out.dir, "\n")

# ------------------------------------------------------------
# 3. LaTeX table (Table 2: all genetically assigned herring)
# ------------------------------------------------------------
## booktabs only -- no multirow, no threeparttable, so it drops
## into any preamble. The year is printed once per block and left
## blank on the following rows. NEATL and NSS get their own columns.
SET <- "genetic total"
tex_stocks <- c(STACK, TINY)
tiny_note <- c("    absent from that year--subdivision sample. NEATL and NSS are listed for",
               "    completeness; they were excluded from the shape analyses because of their",
               "    very small numbers.")
note_scope <- "fish with a genetic population assignment, regardless of whether a usable photograph of the right sagitta was available"

tex <- long %>%
  filter(Set == SET) %>%
  mutate(cell = ifelse(n == 0, "--", sprintf("%s (%.1f)", formatC(n, format = "d", big.mark = ","), percent))) %>%
  select(Year, SD, Stock, cell) %>%
  pivot_wider(names_from = Stock, values_from = cell) %>%
  arrange(Year, SD)

tot_cell <- long %>% filter(Set == SET) %>% distinct(Year, SD, cell_n)
tex <- left_join(tex, tot_cell, by = c("Year", "SD"))

## grand total row: counts over the whole set, percentages of the whole set
gt <- long %>%
  filter(Set == SET) %>%
  group_by(Stock) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct = 100 * n / sum(n),
         cell = ifelse(n == 0, "--", sprintf("%s (%.1f)", formatC(n, format = "d", big.mark = ","), pct)))
gt_cells <- setNames(gt$cell, as.character(gt$Stock))[tex_stocks]
gt_total <- sum(gt$n)

body <- character(0)
prev_year <- ""
for (i in seq_len(nrow(tex))) {
  y <- as.character(tex$Year[i])
  if (y != prev_year && i > 1) body <- c(body, "\\addlinespace")
  ycol <- if (y == prev_year) "" else y
  prev_year <- y
  body <- c(body, paste0(
    paste(c(ycol, as.character(tex$SD[i]),
            unlist(tex[i, tex_stocks]),
            format(tex$cell_n[i], big.mark = ",")),
          collapse = " & "), " \\\\"))
}

cap_scope <- "all genetically assigned herring"
cap_short <- "Population composition of all genetically assigned herring"
tab_label <- "tab:genetic-total"

lines <- c(
"\\begin{table}[htbp]",
"  \\centering",
paste0("  \\caption[", cap_short, "]{%"),
paste0("    Population composition of ", cap_scope, " by survey year and"),
"    ICES Subdivision (SD). Population identity was determined genetically",
"    by mixed-stock assignment (posterior probability of population membership);",
"    the otolith shape data were not used for this assignment. Cells give the",
"    number of individuals and, in parentheses, the percentage of the",
"    year--subdivision sample; the last row gives the totals and their",
"    percentage of the whole data set.}",
paste0("  \\label{", tab_label, "}"),
"  \\footnotesize",
"  \\setlength{\\tabcolsep}{4pt}",
"  \\begin{adjustbox}{max width=\\linewidth}",
paste0("  \\begin{tabular}{ll", strrep("r", length(tex_stocks)), "r}"),
"    \\toprule",
paste0("    & & \\multicolumn{", length(tex_stocks),
       "}{c}{Genetically assigned population, $n$ (\\%)} & \\\\"),
paste0("    \\cmidrule(lr){3-", 2 + length(tex_stocks), "}"),
paste0("    Year & SD & ", paste(tex_stocks, collapse = " & "), " & Total \\\\"),
"    \\midrule",
paste0("    ", body),
"    \\midrule",
paste0("    \\multicolumn{2}{l}{All} & ",
       paste(gt_cells, collapse = " & "), " & ", formatC(gt_total, format = "d", big.mark = ","), " \\\\"),
"    \\bottomrule",
"  \\end{tabular}",
"  \\end{adjustbox}",
"",
"  \\vspace{0.5em}",
"  \\begin{minipage}{\\textwidth}",
"    \\footnotesize",
paste0("    \\textit{Note.} $n = ", formatC(gt_total, format = "d", big.mark = "{,}"), "$ ", note_scope, ". Dashes mark populations"),
tiny_note,
"  \\end{minipage}",
"\\end{table}")

tex_file <- "table_composition_genetic_total.tex"
writeLines(lines, file.path(out.dir, tex_file))
cat("\nLaTeX table written to:", file.path(out.dir, tex_file), "\n")
