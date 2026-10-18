# ============================================================
# Figure 14 -- outlier removal
# ------------------------------------------------------------
# Shows the quantity the removal rule is based on: the
# Mahalanobis distance D^2 of every otolith over
# PC1-PC5, one column of jittered points per population, with the
# threshold (chi-square, 5 df, p = 0.001 -> D^2 = 20.5) as a line.
#
# Every removed otolith lies above the line, and the removal
# rates of the populations can be compared directly.
#
# Output: results final plots/results chapter/{pdf,png,data}
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(Momocs)
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
fig.dir <- ensure_dir(file.path(DIR_FINAL, "results chapter"))
rev.pdf <- ensure_dir(file.path(fig.dir, "pdf"))
rev.png <- ensure_dir(file.path(fig.dir, "png"))
rev.dat <- ensure_dir(file.path(fig.dir, "data"))

FAM <- FONT_PLOT
theme_set(theme_bw(base_size = 11, base_family = FAM))
H_USED <- 12

## ---- the pre-removal sample, rebuilt as in main_analysis.R ---
load(file.path(RESULTS, "dataBIO.RData"))   # data, the analysed sample
out_ids <- read.csv(file.path(RESULTS, "outliers.csv"), stringsAsFactors = FALSE)$ID
oto_all <- local({
  e <- new.env()
  load(file.path(RESULTS, "oto_all.RData"), envir = e)
  e$oto
})

id_photo <- tibble(
  ID              = names(oto_all),
  IdentNr.Picture = as.integer(sub(".*__(\\d+)_.*", "\\1", names(oto_all))),
  Year            = as.integer(sub(".*_(\\d{4})_SD.*", "\\1", names(oto_all))),
  SD              = as.integer(sub(".*_SD(\\d+).*", "\\1", names(oto_all))))

bio_all <- read.csv(F_DATA_AGE, row.names = NULL,
                    fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE) %>%
  filter(Assigned.Stock %in% STOCK_ORDER,
         !is.na(IdentNr.Picture), !is.na(Spawning.Component)) %>%
  mutate(Subdivision = as.integer(Subdivision), Year = as.integer(Year)) %>%
  left_join(id_photo, by = c("IdentNr.Picture" = "IdentNr.Picture",
                             "Subdivision"     = "SD",
                             "Year"            = "Year")) %>%
  filter(ID %in% names(oto_all)) %>%
  mutate(removed = ID %in% out_ids)
stopifnot(sum(bio_all$removed) == length(out_ids),
          sum(!bio_all$removed) == nrow(data))

## ---- pre-removal PCA and the D^2 rule ------------------------
slide_to_bottom <- function(coo) coo_slide(coo, id = which.min(coo[, 2]))
d5 <- bio_all %>% arrange(match(ID, names(oto_all)))
O5 <- Out(oto_all[d5$ID], fac = d5)
O5$coo <- lapply(O5$coo, slide_to_bottom)
O5 <- Out(O5$coo, fac = d5); O5$fac <- droplevels(O5$fac)
pca5 <- PCA(efourier(O5, nb.h = H_USED, smooth.it = 0, norm = TRUE))

sc5  <- as.data.frame(pca5$x[, 1:5])
thr5 <- qchisq(0.999, df = 5)
md5  <- mahalanobis(sc5, colMeans(sc5), cov(sc5))
## must be the very PCA the pipeline removed outliers in
stopifnot(setequal(rownames(pca5$x)[md5 > thr5], out_ids))
ve5 <- 100 * pca5$eig / sum(pca5$eig)

pts5 <- data.frame(ID = rownames(pca5$x), Stock = as.character(O5$fac$Assigned.Stock),
                   PC1 = sc5$PC1, PC2 = sc5$PC2, D2 = md5,
                   removed = md5 > thr5)
write.csv(pts5 %>% filter(removed) %>%
            mutate(D2 = round(D2, 2), PC1 = round(PC1, 4), PC2 = round(PC2, 4)) %>%
            select(-removed) %>% arrange(Stock, desc(D2)),
          file.path(rev.dat, "fig5_outliers_d2_data.csv"), row.names = FALSE)

cnt5 <- pts5 %>% group_by(Stock) %>%
  summarise(n = n(), r = sum(removed), .groups = "drop")

## Did the removal affect some populations more than others?
## Chi-square test of independence on removed vs retained per population
## (Section 2.5; reported in Section 3.3)
chi5 <- chisq.test(cbind(cnt5$r, cnt5$n - cnt5$r))
cat("removed per population:",
    paste(sprintf("%s %.1f%%", cnt5$Stock, 100 * cnt5$r / cnt5$n), collapse = ", "), "\n")
cat("overall:", sum(cnt5$r), "of", sum(cnt5$n), sprintf("(%.1f%%)", 100 * sum(cnt5$r) / sum(cnt5$n)), "\n")
cat(sprintf("chi-square = %.2f, df = %d, p = %.2f\n", chi5$statistic, chi5$parameter, chi5$p.value))
lab5 <- setNames(sprintf("%s — %s of %s removed", cnt5$Stock,
                         formatC(cnt5$r, format = "d", big.mark = ","),
                         formatC(cnt5$n, format = "d", big.mark = ",")), cnt5$Stock)
pts5$Stock <- factor(pts5$Stock, levels = STOCK_ORDER)

## ---- the figure --------------------------------------------
xlab5 <- setNames(sprintf("%s\n%s of %s removed", cnt5$Stock,
                          formatC(cnt5$r, format = "d", big.mark = ","),
                          formatC(cnt5$n, format = "d", big.mark = ",")), cnt5$Stock)
cat("D2 range:", signif(range(pts5$D2), 3), "\n")

COL_OUT  <- COL_STOCK[["WBSS"]]   # orange
COL_KEPT <- "#BAC2CA"             # retained otoliths

set.seed(1)                        # fixed jitter, same picture every run
pts5$xj <- as.integer(pts5$Stock) + runif(nrow(pts5), -0.32, 0.32)
pts5$grp <- ifelse(pts5$removed, "out", "kept")

fig <- ggplot(pts5, aes(xj, D2)) +
  geom_hline(aes(yintercept = thr5, linetype = "thr"), colour = COL_INK, linewidth = 0.45) +
  geom_point(data = ~ filter(.x, !removed), aes(fill = grp),
             shape = 21, colour = COL_KEPT, size = 0.6, stroke = 0.15) +
  geom_point(data = ~ filter(.x, removed), aes(fill = grp),
             shape = 21, colour = "white", size = 1.5, stroke = 0.3) +
  scale_fill_manual(name = NULL, values = c(kept = COL_KEPT, out = COL_OUT),
    breaks = c("out", "kept"), labels = c(out = "removed", kept = "retained"),
    guide = guide_legend(order = 1,
                         override.aes = list(size = 1.8, colour = "white", stroke = 0.3))) +
  scale_linetype_manual(name = NULL, values = c(thr = "22"),
    labels = c(thr = sprintf("removal threshold (D\u00b2 = %.1f)", thr5)),
    guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = 1:5, labels = xlab5[STOCK_ORDER],
                     limits = c(0.5, 5.5), expand = expansion(0)) +
  ## log scale: retained and removed otoliths both get room, the threshold
  ## sits mid-height. (A square-root scale squeezed the retained otoliths
  ## into a thin strip at the bottom.)
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(x = NULL, y = "Mahalanobis distance D\u00b2 (PC1\u2013PC5, log scale)") +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey92"),
        panel.border       = element_rect(colour = "grey80", linewidth = 0.4),
        axis.text.x        = element_text(colour = COL_INK, size = 9.5, lineheight = 1.1),
        axis.text.y        = element_text(colour = COL_MUT, size = 8.5),
        axis.title         = element_text(colour = COL_INK, size = 10),
        axis.ticks.x       = element_blank(),
        legend.position    = "top", legend.justification = "right",
        legend.text        = element_text(colour = COL_INK, size = 9.5),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"),
        legend.key.width   = unit(7, "mm"))

name <- "fig5_outliers_d2"
ggsave(file.path(rev.pdf, paste0(name, ".pdf")), fig,
       width = 16, height = 10, units = "cm", device = cairo_pdf)
ggsave(file.path(rev.png, paste0(name, ".png")), fig,
       width = 16, height = 10, units = "cm", dpi = 600, device = ragg::agg_png)
cat("Done:", fig.dir, "\n")
