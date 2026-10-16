# ============================================================
# Add Maturity column to Otolith_Shape_DATA_SET_with_Age.csv
# Main join key: Year + IdentNr
# ============================================================

# ------------------------------------------------------------
# Path setup
# ------------------------------------------------------------
library(tidyverse)

## ---- load project paths ------------------------------------
## Searches upwards from the working directory for config.R (also in a
## subfolder Master-Thesis/) and sources it. Afterwards SHAPE_DIR,
## DIR_IMAGES, DIR_FOURIER, ensure_dir() etc. are available.
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- c(file.path(d, "config.R"), file.path(d, "Master-Thesis", "config.R"))
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) {
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    }
    d <- dirname(d)
  }
})

main_dir <- DIR_DATASETS

main_file <- file.path(main_dir, "Otolith_Shape_DATA_SET_with_Age.csv")
file_2023 <- file.path(main_dir, "clean_Herring_genetic_SB827_FB.csv")
file_2024 <- file.path(main_dir, "GERAS2024_MSA.csv")
file_2025 <- file.path(main_dir, "GERAS_2025 assignments_partial.csv")
output_file <- file.path(main_dir, "Otolith_Shape_DATA_SET_with_Age_Maturity.csv")

cat("Checking files...\n")
cat("Main file exists: ", file.exists(main_file), "\n")
cat("2023 file exists: ", file.exists(file_2023), "\n")
cat("2024 file exists: ", file.exists(file_2024), "\n")
cat("2025 file exists: ", file.exists(file_2025), "\n\n")

stopifnot(file.exists(main_file))
stopifnot(file.exists(file_2023))
stopifnot(file.exists(file_2024))
stopifnot(file.exists(file_2025))


# ------------------------------------------------------------
# Load main dataset
# ------------------------------------------------------------
main_data <- read.csv(
  main_file,
  row.names = NULL,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

cat("Main dataset loaded.\n")
cat("Rows:", nrow(main_data), "\n")
cat("Columns:\n")
print(names(main_data))

main_data <- main_data %>%
  mutate(
    Year = as.integer(Year),
    IdentNr = as.integer(IdentNr)
  )

# ------------------------------------------------------------
# Load 2023 source
# ------------------------------------------------------------
d2023 <- read.csv(
  file_2023,
  row.names = NULL,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

d2023_sub <- d2023 %>%
  transmute(
    Source = "2023",
    Year = as.integer(Year),
    IdentNr = as.integer(IdentNr),
    Age_source = suppressWarnings(as.numeric(age)),
    Maturity_source = suppressWarnings(as.numeric(maturity))
  )

# ------------------------------------------------------------
# Load 2024 source
# ------------------------------------------------------------
d2024 <- read.csv(
  file_2024,
  row.names = NULL,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

d2024_sub <- d2024 %>%
  transmute(
    Source = "2024",
    Year = as.integer(year),
    IdentNr = as.integer(IDENTNR),
    Age_source = suppressWarnings(as.numeric(age)),
    Maturity_source = suppressWarnings(as.numeric(maturity))
  )

# ------------------------------------------------------------
# Load 2025 source
# This file is semicolon-separated
# ------------------------------------------------------------
d2025 <- read.csv2(
  file_2025,
  row.names = NULL,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

d2025_sub <- d2025 %>%
  transmute(
    Source = "2025",
    Year = as.integer(JAHR),
    IdentNr = as.integer(IDENTNR),
    Age_source = suppressWarnings(as.numeric(ALTER)),
    Maturity_source = suppressWarnings(as.numeric(REIFE))
  )

# ------------------------------------------------------------
# Combine all source tables
# ------------------------------------------------------------
source_all <- bind_rows(d2023_sub, d2024_sub, d2025_sub)

cat("\nCombined source rows:", nrow(source_all), "\n")

source_all <- source_all %>%
  filter(!(is.na(Age_source) & is.na(Maturity_source)))

cat("Rows after removing empty source rows:", nrow(source_all), "\n")

# ------------------------------------------------------------
# Check duplicate Year + IdentNr combinations
# ------------------------------------------------------------
dup_keys <- source_all %>%
  count(Year, IdentNr) %>%
  filter(n > 1)

cat("Duplicate Year + IdentNr combinations:", nrow(dup_keys), "\n")

# ------------------------------------------------------------
# Collapse duplicates
# Keep first non-NA value
# ------------------------------------------------------------
first_non_na_num <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA_real_)
  y[1]
}

source_unique <- source_all %>%
  group_by(Year, IdentNr) %>%
  summarise(
    Age_source = first_non_na_num(Age_source),
    Maturity_source = first_non_na_num(Maturity_source),
    n_source_rows = n(),
    .groups = "drop"
  )

cat("Unique source rows after collapsing duplicates:", nrow(source_unique), "\n")
cat("Rows with Maturity source:", sum(!is.na(source_unique$Maturity_source)), "\n")


# ------------------------------------------------------------
# Optional check: compare Age in main file vs source Age
# ------------------------------------------------------------
if ("Age" %in% names(main_data)) {
  age_check <- main_data %>%
    left_join(
      source_unique %>% select(Year, IdentNr, Age_source),
      by = c("Year", "IdentNr")
    ) %>%
    mutate(
      Age = suppressWarnings(as.numeric(Age))
    )
  
  matched_age <- age_check %>%
    filter(!is.na(Age), !is.na(Age_source))
  
  cat("\nAge consistency check:\n")
  cat("Rows with both Age and Age_source:", nrow(matched_age), "\n")
  
  if (nrow(matched_age) > 0) {
    cat("Age matches:", sum(matched_age$Age == matched_age$Age_source, na.rm = TRUE), "\n")
    cat("Age differs :", sum(matched_age$Age != matched_age$Age_source, na.rm = TRUE), "\n")
  }
}

# ------------------------------------------------------------
# Join Maturity into main dataset
# ------------------------------------------------------------
main_joined <- main_data %>%
  left_join(
    source_unique %>% select(Year, IdentNr, Maturity_source),
    by = c("Year", "IdentNr")
  ) %>%
  rename(Maturity = Maturity_source)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
cat("\nJoin summary:\n")
cat("Rows in main dataset:", nrow(main_joined), "\n")
cat("Rows with Maturity:", sum(!is.na(main_joined$Maturity)), "\n")
cat("Rows without Maturity:", sum(is.na(main_joined$Maturity)), "\n")

cat("\nMaturity distribution:\n")
print(table(main_joined$Maturity, useNA = "ifany"))

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
write.csv(
  main_joined,
  file = output_file,
  row.names = FALSE,
  na = ""
)

cat("\nFinished.\n")
cat("New file written:\n")
cat(output_file, "\n")

