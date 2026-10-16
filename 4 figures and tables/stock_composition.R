# ============================================================
# Stock composition by year and ICES subdivision
# ------------------------------------------------------------
# Descriptive figure + table for the opening of the results
# chapter: which stocks were caught where, in which year.
#
# Two things this script is careful about:
#
#  1. Percentages are computed WITHIN a year x subdivision cell,
#     never within a year alone. Sampling effort moved between
#     subdivisions across the three surveys (SD23: 67 -> 476 ->
#     289 otoliths), and SD23 is WBSS-rich, so a percentage
#     pooled over subdivisions mostly reports where the survey
#     sampled, not how the stocks changed.
#
#  2. It reports the ANALYSED shape sample (five stocks; NEATL and
#     NSS are excluded from the analysis) and the
#     FULL genetically assigned set (n = 3386) side by side, so
#     the text can state that dropping otoliths without a usable
#     photograph did not bias the composition.
#
# Colours are the fixed per-population scheme from config.R, so
# the reader can carry a colour from this figure into every
# later plate.
#
# Writes to: Fourier and Classification/Composition/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
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

## Stack order = POP_ORDER from config.R: WBSS first, so it sits on
## the baseline -- it is the focal population and the largest class,
## and a segment anchored to the axis can be read off directly in all
## twelve bars without mental subtraction. Then BAS, NSAS, Downs, CBH.
## The config.R palette is checked over all pairs, so every
## neighbouring pair in the stack stays distinguishable.
STACK <- POP_ORDER

## NEATL (n = 2) and NSS (n = 3) are excluded from the shape
## analysis itself (main_analysis.R), so the analysed sample
## holds none of them. They remain only in the genetic-total set,
## and only in its table -- never in the figure.
TINY <- c("NEATL", "NSS")

## Colours: COL_STOCK from config.R.

INK  <- "#101821"   # all text -- never the series colour
MUTED <- "#6B7885"

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

## Wide count table, ready to drop into the thesis
wide <- long %>%
  filter(Set == "shape sample", Stock %in% STACK) %>%
  mutate(entry = ifelse(n == 0, "–", sprintf("%d (%.1f)", n, percent))) %>%
  select(Year, SD, Stock, entry) %>%
  pivot_wider(names_from = Stock, values_from = entry) %>%
  left_join(count(A, Year, SD, name = "Total"), by = c("Year", "SD")) %>%
  arrange(Year, SD)

write.csv(wide, file.path(out.dir, "table_composition_counts.csv"), row.names = FALSE)
cat("--- counts (percent within year x subdivision) ---\n")
print(as.data.frame(wide), row.names = FALSE)

## Does dropping un-photographed fish shift the composition?
bias <- long %>%
  filter(Stock %in% STACK) %>%
  select(Set, Year, SD, Stock, percent) %>%
  pivot_wider(names_from = Set, values_from = percent) %>%
  mutate(diff_pp = round(`shape sample` - `genetic total`, 1))
cat("\nlargest composition shift caused by photo availability:",
    sprintf("%.1f pp\n", max(abs(bias$diff_pp), na.rm = TRUE)))
write.csv(bias, file.path(out.dir, "coverage_bias_check.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 3. The figure
# ------------------------------------------------------------
## 100 % stacked bars, one panel per subdivision, three bars per
## panel. Composition is the message, so the bars are normalised;
## the sample size sits above each bar so the reader can see how
## thin a percentage is before trusting it.
plot_df <- long %>%
  filter(Set == "shape sample", Stock %in% STACK) %>%
  mutate(Stock = factor(Stock, levels = STACK),
         SD    = factor(SD, labels = paste("SD", sort(unique(SD)))))

n_lab <- plot_df %>%
  distinct(Year, SD, cell_n)

p <- ggplot(plot_df, aes(x = factor(Year), y = percent, fill = Stock)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.5,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = n_lab,
            aes(x = factor(Year), y = 103, label = paste0("n=", cell_n)),
            inherit.aes = FALSE, size = 2.9, colour = MUTED, vjust = 0) +
  facet_wrap(~SD, nrow = 1) +
  scale_fill_manual(values = COL_STOCK, name = NULL) +
  scale_y_continuous(limits = c(0, 111), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Share of otoliths in the subdivision") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
    panel.border       = element_rect(colour = "grey80", linewidth = 0.4),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", colour = INK, size = 10.5),
    axis.text          = element_text(colour = INK),
    axis.title.y       = element_text(colour = INK, margin = margin(r = 8)),
    legend.position    = "bottom",
    legend.key.size    = unit(4, "mm"),
    legend.text        = element_text(colour = INK)
  )

ggsave(file.path(out.dir, "stock_composition.tiff"), p,
       width = 18, height = 8.5, units = "cm", dpi = 600, compression = "lzw")
ggsave(file.path(out.dir, "stock_composition.pdf"), p,
       width = 18, height = 8.5, units = "cm")

cat("\nWritten to:", out.dir, "\n")

# ------------------------------------------------------------
# 4. LaTeX table
# ------------------------------------------------------------
## booktabs only -- no multirow, no threeparttable, so it drops
## into any preamble. The year is printed once per block and left
## blank on the following rows.
##
## SET <- "genetic total" writes the same table for all
## genetically assigned fish instead of the analysed ones.
if (!exists("SET")) SET <- "shape sample"

## The analysed sample has only the five stocks; the genetic total
## also carries NEATL and NSS, which get their columns there.
tex_stocks <- if (SET == "shape sample") STACK else c(STACK, TINY)
tiny_note <- if (SET == "shape sample") {
  c("    absent from that year--subdivision sample. NEATL (2 fish) and NSS (3 fish)",
    "    were excluded from the shape analyses because of their very small numbers.")
} else {
  c("    absent from that year--subdivision sample. NEATL and NSS are listed for",
    "    completeness; they were excluded from the shape analyses because of their",
    "    very small numbers.")
}

## The note has to describe whichever set SET selects, or the
## genetic-total variant claims an outline the fish may not have.
## ("genetic population assignment" = Assigned.Stock is one of VALID)
note_scope <- if (SET == "shape sample") {
  "otoliths with both a genetic population assignment and a usable outline of the right sagitta"
} else {
  "fish with a genetic population assignment, regardless of whether a usable photograph of the right sagitta was available"
}

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

## caption and label follow the set, so the two versions of this table
## never describe the wrong sample or share a \\label
cap_scope <- if (SET == "shape sample") "the analysed otolith sample" else
             "all genetically assigned herring"
cap_short <- if (SET == "shape sample") "Population composition of the analysed sample" else
             "Population composition of all genetically assigned herring"
tab_label <- if (SET == "shape sample") "tab:genetic-composition" else "tab:genetic-total"

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

## one file per set, so the genetic-total version no longer has to be
## renamed by hand
tex_file <- if (SET == "shape sample") "table_composition.tex" else
            "table_composition_genetic_total.tex"
writeLines(lines, file.path(out.dir, tex_file))
cat("\nLaTeX table written to:", file.path(out.dir, tex_file), "\n")
