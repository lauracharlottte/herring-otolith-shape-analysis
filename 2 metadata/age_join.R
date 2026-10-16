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
age2025 <- read_excel(file.path(in.dir, "GERAS2025_SB855_2025_SnglIOR_WITHgenID_andPopassign_for Leo.xlsx")) %>%
  mutate(Haul = as.character(station), IdentNr = as.character(IDENTNR)) %>%
  select(Haul, IdentNr, Age = age)

## ---- 5. Combine all years into one age lookup table ----
age_all <- bind_rows(
  age2023 %>% mutate(Age = as.character(Age)),
  age2024 %>% mutate(Age = as.character(Age)),
  age2025 %>% mutate(Age = as.character(Age))
) %>%
  mutate(Age = as.numeric(Age))  # turns any non-numeric entries (e.g. "NA", "NB") into NA

## ---- 6. Join onto main dataset ----
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

## ---- 7. DIAGNOSTIC: check match rate per year ----
match_summary <- data %>%
  group_by(Year) %>%
  summarise(
    n_total = n(),
    n_matched = sum(!is.na(Age)),
    n_missing = sum(is.na(Age)),
    pct_matched = round(100 * n_matched / n_total, 1)
  )
print(match_summary)

## ---- 8. See exactly which rows failed to match, for spot-checking ----
unmatched <- data %>% filter(is.na(Age))
write.csv(unmatched, file.path(in.dir, "unmatched_age.csv"), row.names = FALSE)
cat("Number of unmatched rows:", nrow(unmatched), "-- see unmatched_age.csv\n")

## ---- 9. Save the age-annotated dataset ----
write.csv(data, file.path(in.dir, "Otolith_Shape_DATA_SET_with_Age.csv"), row.names = FALSE)
save(data, file = file.path(in.dir, "dataBIO_with_Age.RData"))

message("Age_Join.R complete. Check match_summary and unmatched_age.csv before proceeding.")
