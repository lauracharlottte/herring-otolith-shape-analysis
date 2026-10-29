# ============================================================
# Appendix tables -- age composition of the analysed otoliths
# ------------------------------------------------------------
#   Table A1  number of otoliths per population and age group
#             (in parentheses: % of the otoliths of each population)
#   Table A2  number of otoliths per age group and ICES Subdivision,
#             with the share of each subdivision for all ages and for
#             the otoliths aged 1-3 years
# Age code -9 (not read) counts as unknown.
# Output: Shape Analysis/Appendix/tables/table_A1_population_age.csv,
#         table_A2_age_subdivision.csv
# ============================================================
suppressPackageStartupMessages(library(dplyr))
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- file.path(d, "config.R")
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) stop("config.R not found -- please setwd() into the repository folder.", call. = FALSE)
    d <- dirname(d)
  }
})
RESULTS <- DIR_RESULTS
out <- ensure_dir(proj_path("Appendix", "tables"))
load(file.path(RESULTS, "dataBIO.RData"))   # data, the analysed sample

age <- suppressWarnings(as.numeric(as.character(data$Age))); age[age < 0] <- NA
grp <- cut(age, c(-Inf, 0, 1, 2, 3, 5, Inf), labels = c("0", "1", "2", "3", "4–5", "6–11"))
grp <- factor(ifelse(is.na(grp), "Unknown", as.character(grp)),
              levels = c("0", "1", "2", "3", "4–5", "6–11", "Unknown"))
pop <- factor(data$Assigned.Stock, levels = POP_ORDER)
sd  <- factor(as.integer(data$Subdivision))

## Table A1: population x age group, n (% of the population)
t1 <- table(pop, grp)
a1 <- sapply(levels(grp)[1:6], function(g) sprintf("%s (%.1f)", format(t1[, g], big.mark = ",", trim = TRUE),
                                                   100 * t1[, g] / rowSums(t1)))
a1 <- data.frame(Population = levels(pop), a1, Total = rowSums(t1), check.names = FALSE)
write.csv(a1, file.path(out, "table_A1_population_age.csv"), row.names = FALSE)
print(a1, row.names = FALSE)

## Table A2: age group x subdivision, plus the subdivision shares of all ages and of age 1-3
t2 <- table(grp, sd)
a2 <- data.frame(Age = levels(grp), unclass(t2), Total = rowSums(t2), check.names = FALSE)
names(a2)[2:5] <- paste("SD", levels(sd))
all13 <- colSums(t2[c("1", "2", "3"), ])
share <- function(x) sprintf("%s (%.1f)", format(x, big.mark = ",", trim = TRUE), 100 * x / sum(x))
a2[] <- lapply(a2, as.character)
extra <- as.data.frame(rbind(c("All ages", share(colSums(t2)), sum(t2)),
                             c("Age 1–3", share(all13), sum(all13))), stringsAsFactors = FALSE)
names(extra) <- names(a2)
a2 <- rbind(a2, extra)
write.csv(a2, file.path(out, "table_A2_age_subdivision.csv"), row.names = FALSE)
print(a2, row.names = FALSE)
