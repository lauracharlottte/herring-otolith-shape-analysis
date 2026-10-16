# ============================================================
# Appendix -- mean otolith outlines of the age subsets 0-3 and 1-10
# ------------------------------------------------------------
# Same design as Figures 16 and 17 (fig_composition_and_mean_shapes.R):
# the mean outline of every population against the mean outline of
# WBSS (filled), here for (a) otoliths aged 0-3 and (b) aged 1-10,
# one row per subset.
# Output: Shape Analysis/Appendix/{pdf,png}/figA_mean_shapes_age0_3_and_age1_10
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(Momocs) })
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- c(file.path(d, "config.R"), file.path(d, "Master-Thesis", "config.R"))
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) stop("config.R not found -- please setwd() into the repository folder.", call. = FALSE)
    d <- dirname(d)
  }
})
RESULTS <- DIR_RESULTS
out <- proj_path("Appendix")
e6 <- local({ e <- new.env(); load(file.path(RESULTS, "eOTO.RData"), envir = e); e$eOTO })
age <- suppressWarnings(as.numeric(as.character(e6$fac$Age)))
REF <- "WBSS"

shapes <- function(lo, hi, lab) {
  idx <- which(!is.na(age) & age >= lo & age <= hi)
  eo <- e6; eo$coe <- e6$coe[idx, , drop = FALSE]; eo$fac <- e6$fac[idx, , drop = FALSE]   # same rows as slice()
  ms <- MSHAPES(eo, fac = "Assigned.Stock")$shp
  lv <- intersect(POP_ORDER, names(ms))
  shp <- bind_rows(lapply(lv, function(g) { s <- as.data.frame(ms[[g]]); s <- rbind(s, s[1, ]); s$Group <- g; s }))
  oth <- setdiff(lv, REF)
  bind_rows(lapply(oth, function(g) bind_rows(
    filter(shp, Group == REF) %>% mutate(Panel = g),
    filter(shp, Group == g) %>% mutate(Panel = g)))) %>%
    mutate(Subset = lab)
}
pd <- bind_rows(shapes(0, 3, "(a) Age 0–3"), shapes(1, 10, "(b) Age 1–10"))
oth <- setdiff(POP_ORDER, REF)
pd <- pd %>% mutate(Panel = factor(Panel, levels = oth), Group = factor(Group, levels = POP_ORDER),
                    Subset = factor(Subset, levels = unique(Subset)))

f <- ggplot(pd, aes(x, y)) +
  geom_polygon(data = ~ filter(.x, Group == REF), fill = COL_STOCK[[REF]], alpha = 0.18, colour = NA) +
  geom_path(aes(colour = Group), linewidth = 0.6) +
  facet_grid(Subset ~ Panel, switch = "y",
             labeller = labeller(Panel = setNames(paste(REF, "and", oth), oth))) +
  scale_colour_manual(values = COL_STOCK[POP_ORDER], limits = POP_ORDER, name = "Population") +
  coord_equal() +
  theme_void(base_size = 11, base_family = FONT_PLOT) +
  theme(legend.position = "bottom", legend.key.width = unit(8, "mm"),
        legend.title = element_text(colour = COL_INK, size = 10),
        legend.text  = element_text(colour = COL_INK, size = 10),
        strip.text.x = element_text(colour = COL_INK, size = 10, margin = margin(b = 4)),
        strip.text.y.left = element_text(colour = COL_INK, size = 10, angle = 90, margin = margin(r = 6)),
        panel.spacing = unit(4, "mm"), strip.clip = "off", strip.placement = "outside",
        plot.margin = margin(4, 6, 4, 10))
ensure_dir(file.path(out, "pdf")); ensure_dir(file.path(out, "png"))
ggsave(file.path(out, "pdf", "figA_mean_shapes_age0_3_and_age1_10.pdf"), f, width = 17, height = 7.5, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(out, "png", "figA_mean_shapes_age0_3_and_age1_10.png"), f, width = 17, height = 7.5, units = "cm", dpi = 600, device = ragg::agg_png, bg = "white")
cat("done\n")
