############################################################################
## Age_Join.R
## Purpose: Join age data (from 3 separate yearly sources) onto the main
## otolith dataset, without altering the original CSV.
############################################################################

library(dplyr)
library(readxl)
library(readr)

## ---- 1. Load main dataset (from the Outline Extraction folder) ----
## ---- load project paths ------------------------------------
## Searches upwards from the working directory for config.R
## and sources it. Afterwards SHAPE_DIR,
## DIR_IMAGES, DIR_FOURIER, ensure_dir() etc. are available.
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- file.path(d, "config.R")
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) {
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    }
    d <- dirname(d)
  }
})

## Folder with the age source files; the outputs are written here as well.
in.dir <- DIR_DATASETS

data <- read.csv2(file.path(DIR_DATASETS, "Otolith_Shape_DATA_SET.csv"),
                  row.names = NULL, fileEncoding = "UTF-8-BOM")

## ---- 2. 2023 ages ----
age2023 <- read.csv(file.path(in.dir, "clean_Herring_genetic_SB827_FB.csv"), stringsAsFactors = FALSE) %>%
  select(indivCODE, age) %>%
  rename(Age = age)

## ---- 3. 2024 ages ----
age2024 <- read.csv(file.path(in.dir, "GERAS2024_MSA.csv"), stringsAsFactors = FALSE) %>%
  select(TUBE_ID, age) %>%
  rename(indivCODE = TUBE_ID, Age = age)

## ---- 4. 2025 ages ----
## (the 2025 file is found by its prefix; the full name is an internal file name)
f2025 <- list.files(in.dir, pattern = "^GERAS2025_SB855.*\\.xlsx$", full.names = TRUE)
stopifnot(length(f2025) == 1)
age2025 <- read_excel(f2025) %>%
  mutate(Haul = as.character(station), IdentNr = as.character(IDENTNR)) %>%
  select(Haul, IdentNr, Age = age)

## ---- 5. Join onto main dataset ----
## 2023 and 2024 are matched by the individual code, 2025 by haul and
## fish number; non-numeric age entries (e.g. "NA", "NB") become NA.
data <- data %>%
  mutate(Haul = as.character(Haul), IdentNr = as.character(IdentNr)) %>%
  left_join(
    bind_rows(
      age2023 %>% mutate(Age = as.character(Age)),
      age2024 %>% mutate(Age = as.character(Age))
    ) %>% mutate(Age = as.numeric(Age)),
    by = "indivCODE"
  ) %>%
  left_join(age2025, by = c("Haul", "IdentNr")) %>%
  mutate(Age = coalesce(Age.x, Age.y)) %>%
  select(-Age.x, -Age.y)

## ---- 6. DIAGNOSTIC: check match rate per year ----
match_summary <- data %>%
  group_by(Year) %>%
  summarise(
    n_total = n(),
    n_matched = sum(!is.na(Age)),
    n_missing = sum(is.na(Age)),
    pct_matched = round(100 * n_matched / n_total, 1)
  )
print(match_summary)

## ---- 7. See exactly which rows failed to match, for spot-checking ----
unmatched <- data %>% filter(is.na(Age))
write.csv(unmatched, file.path(in.dir, "unmatched_age.csv"), row.names = FALSE)
cat("Number of unmatched rows:", nrow(unmatched), "-- see unmatched_age.csv\n")

## ---- 8. Save the age-annotated dataset ----
write.csv(data, file.path(in.dir, "Otolith_Shape_DATA_SET_with_Age.csv"), row.names = FALSE)
save(data, file = file.path(in.dir, "dataBIO_with_Age.RData"))

message("Age_Join.R complete. Check match_summary and unmatched_age.csv before proceeding.")
