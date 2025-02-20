
input_dir <- Sys.getenv("INPUT_DIR", unset = "/data/input")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "/data/output")

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)

options(chromote.chrome = "/usr/bin/chromium")

# Load the dataset and assign it to nsduh_2023
input_file <- file.path(input_dir, "NSDUH_2023_Tab.txt")
nsduh_2023 <- read_delim(input_file, show_col_types = FALSE)


# Filter rows where the age group is 18-64 and variables of interest

nsduh_2023_s <- nsduh_2023 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR, IRSUIPLANYR, MHNTENFCV, MHNTFFLKE, MHNTPROBS, MHNTTIME, MHNTINSCV, MHNTWHER, MHNTNOHLP, MHNTCOST, MHNTHNDL, MHNTFORCE, MHNTCONSQ, MHNTPRIV, MHNTPTHNK)


# Filter rows with suicidal ideation

nsduh_2023_ss <- nsduh_2023_s |>
  filter(
    ( IRSUICTHNK == 1 | IRSUIPLANYR == 1 | IRSUITRYYR == 1)
  )


#----------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2022

# Construct the file path for the dataset
input_file <- file.path(input_dir, "NSDUH_2022_Tab.txt")
nsduh_2022 <- read_delim(input_file, show_col_types = FALSE)



# Filter rows where the age group is 18-64 and variables of interest


nsduh_2022_s <- nsduh_2022 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR, IRSUIPLANYR, MHNTENFCV, MHNTFFLKE, MHNTPROBS, MHNTTIME, MHNTINSCV, MHNTWHER, MHNTNOHLP, MHNTCOST, MHNTHNDL, MHNTFORCE, MHNTCONSQ, MHNTPRIV, MHNTPTHNK)


# Filter rows with suicidal ideation and violent behaviours

nsduh_2022_ss <- nsduh_2022_s |>
  filter(
    ( IRSUICTHNK == 1 | IRSUIPLANYR == 1 | IRSUITRYYR == 1)
  )


#----------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2021


# Construct the file path for the dataset
input_file <- file.path(input_dir, "NSDUH_2021_Tab.txt")
nsduh_2021 <- read_delim(input_file, show_col_types = FALSE)



# Filter rows where the age group is 18-64 and variables of interest

nsduh_2021_s <- nsduh_2021 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR, IRSUIPLANYR, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2021_ss <- nsduh_2021_s |>
  filter(
    (IRSUICTHNK == 1 | IRSUIPLANYR == 1 | IRSUITRYYR == 1)   )


#----------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2020

# Construct the file path for the dataset
input_file <- file.path(input_dir, "NSDUH_2020_Tab.txt")
nsduh_2020 <- read_delim(input_file, show_col_types = FALSE)



# Filter rows where the age group is 18-64 and variables of interest

nsduh_2020_s <- nsduh_2020 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2020_ss <- nsduh_2020_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )


#----------------------------------------------------------------------------


# Load the dataset and assign it to nsduh_2019

# Construct the file path for the dataset
input_file <- file.path(input_dir, "NSDUH_2019_Tab.txt")
nsduh_2019 <- read_delim(input_file, show_col_types = FALSE)


# Filter rows where the age group is 18-64 and variables of interest

nsduh_2019_s <- nsduh_2019 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours


nsduh_2019_ss <- nsduh_2019_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )


#-----------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2018

input_file <- file.path(input_dir, "NSDUH_2018_Tab.tsv")
nsduh_2018 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)



# Filter rows where the age group is 18-64 and variables of interest

nsduh_2018_s <- nsduh_2018 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2018_ss <- nsduh_2018_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )


#-----------------------------------------------------------------------------


# Load the dataset and assign it to nsduh_2017

input_file <- file.path(input_dir, "NSDUH_2017_Tab.tsv")
nsduh_2017 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)


# Filter rows where the age group is 18-64 and variables of interest

nsduh_2017_s <- nsduh_2017 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2017_ss <- nsduh_2017_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )

#----------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2016

input_file <- file.path(input_dir, "NSDUH_2016_Tab.tsv")
nsduh_2016 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)


# Filter rows where the age group is 18-64 and variables of interest

nsduh_2016_s <- nsduh_2016 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2016_ss <- nsduh_2016_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )


#------------------------------------------------------------------------------

# Load the dataset and assign it to nsduh_2015
input_file <- file.path(input_dir, "NSDUH_2015_Tab.tsv")
nsduh_2015 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)


# Filter rows where the age group is 18-64 and variables of interest

nsduh_2015_s <- nsduh_2015 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN, MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )


# Filter rows with suicidal ideation and violent behaviours

nsduh_2015_ss <- nsduh_2015_s |>
  filter(
    MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1
  )




library(data.table)

# List of datasets from 2015 to 2021

dataset_list <- list(
  nsduh_2015_ss,
  nsduh_2016_ss,
  nsduh_2017_ss,
  nsduh_2018_ss,
  nsduh_2019_ss,
  nsduh_2020_ss,
  nsduh_2021_ss,
  nsduh_2022_ss,
  nsduh_2023_ss

)

years <- 2015:2023

# Convert datasets to data.tables and add YEAR column

dt_list <- mapply(function(df, year) {
  dt <- as.data.table(df)
  dt[, YEAR := year]
  return(dt)
}, dataset_list, years, SIMPLIFY = FALSE)

# Combine datasets from 2015 to 2023

combined_dt <- rbindlist(dt_list, use.names = FALSE)

# Get column names from the 2023 dataset and add 'YEAR'

colnames_2023 <- c(names(nsduh_2023_ss), "YEAR")

# Assign column names to the combined dataset

setnames(combined_dt, colnames_2023)

# Convert the combined data.table to a data.frame

combined_nsduh_ss <- as.data.frame(combined_dt)

#------------------------------------------------------------------------------

 #List of your datasets from 2015 to 2021

dataset_list <- list(
  nsduh_2015_s,
  nsduh_2016_s,
  nsduh_2017_s,
  nsduh_2018_s,
  nsduh_2019_s,
  nsduh_2020_s,
  nsduh_2021_s,
  nsduh_2022_s,
  nsduh_2023_s
  )
years <- 2015:2023

# Convert datasets to data.tables and add YEAR column

dt_list <- mapply(function(df, year) {
  dt <- as.data.table(df)
  dt[, YEAR := year]
  return(dt) }, dataset_list, years, SIMPLIFY = FALSE)

# Combine datasets

combined_dt <- rbindlist(dt_list, use.names = FALSE)

# Get column names from the 2023 dataset and add 'YEAR'

colnames_2023 <- c(names(nsduh_2023_s), "YEAR")

# Assign column names to the combined dataset

setnames(combined_dt, colnames_2023)

# Convert the combined data.table to a data.frame

combined_nsduh_s <- as.data.frame(combined_dt)

# Create the SI dataset
combined_nsduh_si <- combined_nsduh_s |>
  mutate(
    SI = ifelse(
      IRSUICTHNK == 1 |
      IRSUITRYYR == 1 |
      IRSUIPLANYR == 1,
      1,
      0
    )
  )

# Create the dataset for utilization of mental health services analysis
combined_nsduh_svhb <- combined_nsduh_ss |>
  mutate(
    SI = ifelse(
      IRSUICTHNK == 1 |
      IRSUITRYYR == 1 |
      IRSUIPLANYR == 1,
      1,
      0
    )
  )

# Export the combined_nsduh_svhb CSV
export_file <- file.path(output_dir, "combined_nsduh_svhb.csv")
write.csv(combined_nsduh_svhb, file = export_file, row.names = FALSE)

## ----warning=FALSE, message=FALSE----

library(gtsummary)
library(gt)

combined_nsduh_svhbd <- combined_nsduh_svhb |>
  mutate(
    CATAG6 = factor(
      CATAG6,
      levels = 2:5,
      labels = c(
        "18-25 Years Old",
        "26-34 Years Old",
        "35-49 Years Old",
        "50-64 Years Old"
      )
    ),
    IRSEX = factor(IRSEX, levels = 1:2, labels = c("Male", "Female")),
    NEWRACE2 = factor(
      NEWRACE2,
      levels = 1:7,
      labels = c(
        "White",
        "Black/African American",
        "Native American/Alaskan",
        "Hawaiian/Pacific Islander",
        "Asian",
        "Mixed",
        "Hispanic"
      )
    ),
    EDUHIGHCAT = factor(
      EDUHIGHCAT,
      levels = 1:4,
      labels = c(
        "< High school",
        "High school",
        "Associate degree",
        "College graduate"
      )
    )
  )

summary_table <- combined_nsduh_svhbd |>
  select(
    YEAR,
    CATAG6,
    IRSEX,
    NEWRACE2,
    EDUHIGHCAT
  ) |>
  tbl_summary(
    by = YEAR,
    label = list(
      CATAG6 ~ "Age",
      IRSEX ~ "Sex",
      NEWRACE2 ~ "Race/Hispanicity",
      EDUHIGHCAT ~ "Educational Attainment"
    ),
    type = all_categorical() ~ "categorical",
    statistic = all_categorical() ~ "{n} / {N} ({p}%)",
    missing = "no"
  ) |>
  modify_spanning_header(starts_with("stat_") ~ "**Year**") |>
  modify_caption("**Table 1: Summary of Sociodemographics by Year**") |>
  bold_labels() |>
  modify_header(
    label ~ "**Sociodemographics** (N = 23,996)"
  ) |>
  modify_footnote(
    everything() ~ NA_character_
  )

summary_table_gt <- summary_table |> as_gt()
summary_table_gt <- summary_table_gt |> cols_width(label ~ px(500))
summary_table_gt <- summary_table_gt |> tab_options(table.font.size = "small")

summary_table_gt

library(gt)
library(pagedown)
library(webshot2)

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "/data/output")
output_html <- file.path(output_dir, "tab1.html")
gtsave(summary_table_gt, output_html)




prevalence_sic <- combined_nsduh_si |>
  group_by(YEAR) |>
  summarize(
    total_number = n(),
    si_occurrences = sum(SI == 1, na.rm = TRUE),
    Prevalence = (si_occurrences / total_number) * 100
  ) |>
  select(Year = YEAR, Prevalence)


# Create the line graph for SI
plot_si <- ggplot(prevalence_sic, aes(x = Year, y = Prevalence)) +
  geom_line(color = "#ffdd99", linewidth = 1) +
  geom_point(color = "#ffdd99", size = 2) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(
    title = "Prevalence of Suicidal Ideation",
    x = "Year",
    y = "Prevalence (%)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 12)
  )


library(ggplot2)


output_dir <- Sys.getenv("OUTPUT_DIR", unset = "/data/output")

# Define the output file path
output_file <- file.path(output_dir, "prevalence_si_plot.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Export the plot as a PNG image
ggsave(
  filename = output_file,
  plot = plot_si,
  width = 8,
  height = 6,
  dpi = 300
)




library(scales)



# List of variables to test

variables_of_interest <- c(
  "MHNTENFCV", "MHNTFFLKE", "MHNTPROBS", "MHNTTIME", "MHNTINSCV",
  "MHNTWHER", "MHNTNOHLP", "MHNTCOST", "MHNTHNDL", "MHNTFORCE",
  "MHNTCONSQ", "MHNTPRIV", "MHNTPTHNK"
)

# Map of variable codes to descriptive labels

variable_labels <- c(
  MHNTENFCV = "Insurance coverage not enough",
  MHNTFFLKE = "Social disapproval",
  MHNTPROBS = "Logistics issue",
  MHNTTIME  = "Time constraint",
  MHNTINSCV = "No insurance coverage",
  MHNTWHER  = "Do not know where to get treatment",
  MHNTNOHLP = "Skeptical about treatment",
  MHNTCOST  = "Cost barrier",
  MHNTHNDL  = "Thought able to handle on their own",
  MHNTFORCE = "Worry about treatment commitment",
  MHNTCONSQ = "Concern about social security (loss of job/home/child)",
  MHNTPRIV  = "Confidentiality concern",
  MHNTPTHNK = "Worry about social stigma"
)


# Filter dataset to include SI and the YEAR column

data_filtered <- combined_nsduh_svhb |>
  filter(SI == 1) |>
  select(YEAR, all_of(variables_of_interest))

# Convert the data to long format

data_long <- data_filtered |>
  pivot_longer(
    cols = -YEAR, # All columns except YEAR
    names_to = "Variable",
    values_to = "Count"
  )

# Count the occurrences

summary_counts <- data_long |>
  filter(Count == 1) |>
  group_by(YEAR, Variable) |>
  summarise(Count = n(), .groups = "drop")

# Reverse the order of the years within each bar

summary_counts$YEAR <- factor(summary_counts$YEAR, levels = sort(unique(summary_counts$YEAR), decreasing = TRUE))

# Create the horizontal bar chart

rea_si <- ggplot(summary_counts, aes(x = Variable, y = Count, fill = YEAR)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Reasons for Not Receiving MHS Among U.S. Adults with Suicidal Ideation",
    x = "Reasons",
    y = "Number of Cases",
    fill = "Year"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(
    breaks = seq(0, 6000, by = 500),
    labels = seq(0, 6000, by = 500)
  ) +
  scale_x_discrete(labels = variable_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "left",
    # Key: make sure title is positioned at the plot level and left-aligned
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),
    # Optionally reduce left margin if needed:
    plot.margin = margin(t = 10, r = 10, b = 10, l = 5)
  ) +
  coord_flip()


library(ggplot2)


output_dir <- Sys.getenv("OUTPUT_DIR", unset = "/data/output")

# Define the output file path
output_file <- file.path(output_dir, "rea_si_plot.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Export the plot as a PNG image
ggsave(
  filename = output_file,
  plot = rea_si,
  width = 8,
  height = 6,
  dpi = 300
)


