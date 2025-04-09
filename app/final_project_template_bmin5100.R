

library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(aws.s3)
library(data.table)
library(gtsummary)
library(gt)
library(webshot2)
library(ggplot2)
library(aws.ec2metadata)

environment <- Sys.getenv("ENVIRONMENT", "LOCAL")
bucket      <- Sys.getenv("S3_BUCKET", "local-data-bucket")
region      <- Sys.getenv("AWS_DEFAULT_REGION", "us-east-1")
input_dir   <- Sys.getenv("INPUT_DIR", unset = "/data/input")
output_dir  <- Sys.getenv("OUTPUT_DIR", unset = "/data/output")

message(sprintf("Running with ENVIRONMENT=%s, S3_BUCKET=%s, REGION=%s", environment, bucket, region))

if (toupper(environment) == "ECS") {
  message("Fetching ECS metadata credentials...")
  creds <- tryCatch(
    {
      metadata_credentials()
    },
    error = function(e) {
      message(sprintf("Failed to fetch ECS metadata credentials: %s", e$message))
      return(NULL)
    }
  )
  if (!is.null(creds)) {
    Sys.setenv("AWS_ACCESS_KEY_ID" = creds$AccessKeyId,
               "AWS_SECRET_ACCESS_KEY" = creds$SecretAccessKey,
               "AWS_SESSION_TOKEN" = creds$SessionToken,
               "AWS_DEFAULT_REGION" = region)
    message(sprintf("ECS Credentials set: AccessKeyId=%s", creds$AccessKeyId))
  } else {
    message("ECS metadata credentials not available; falling back to environment variables")
  }
} else {
  message("Not in ECS environment; relying on default credentials")
}

if (!dir.exists(input_dir))  dir.create(input_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (toupper(environment) != "LOCAL") {
  message(sprintf("ENVIRONMENT=%s => Downloading input files from s3://%s/", environment, bucket))

  s3_files <- c(
    "NSDUH_2023_Tab.txt",
    "NSDUH_2022_Tab.txt",
    "NSDUH_2021_Tab.txt",
    "NSDUH_2020_Tab.txt",
    "NSDUH_2019_Tab.txt",
    "NSDUH_2018_Tab.tsv",
    "NSDUH_2017_Tab.tsv",
    "NSDUH_2016_Tab.tsv",
    "NSDUH_2015_Tab.tsv"
  )

  for (f in s3_files) {
    s3_object  <- f
    local_path <- file.path(input_dir, f)
    message(sprintf("  Downloading s3://%s/%s -> %s", bucket, s3_object, local_path))
    tryCatch(
      {
        save_object(
          object = s3_object,
          bucket = bucket,
          file   = local_path,
          region = region
        )
        message(sprintf("  Successfully downloaded %s", f))
      },
      error = function(e) {
        message(sprintf("  Error downloading %s: %s", f, e$message))
      }
    )
  }
} else {
  message("ENVIRONMENT=LOCAL => skipping S3 download.")
}



message("Beginning data import and processing...")

options(chromote.chrome = "/usr/bin/chromium")



# 2023
input_file <- file.path(input_dir, "NSDUH_2023_Tab.txt")
nsduh_2023 <- read_delim(input_file, show_col_types = FALSE)

nsduh_2023_s <- nsduh_2023 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR,
    IRSUIPLANYR, MHNTENFCV, MHNTFFLKE, MHNTPROBS, MHNTTIME, MHNTINSCV,
    MHNTWHER, MHNTNOHLP, MHNTCOST, MHNTHNDL, MHNTFORCE, MHNTCONSQ,
    MHNTPRIV, MHNTPTHNK
  )

nsduh_2023_ss <- nsduh_2023_s |>
  filter(IRSUICTHNK == 1 | IRSUITRYYR == 1 | IRSUIPLANYR == 1)

# 2022
input_file <- file.path(input_dir, "NSDUH_2022_Tab.txt")
nsduh_2022 <- read_delim(input_file, show_col_types = FALSE)

nsduh_2022_s <- nsduh_2022 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR,
    IRSUIPLANYR, MHNTENFCV, MHNTFFLKE, MHNTPROBS, MHNTTIME, MHNTINSCV,
    MHNTWHER, MHNTNOHLP, MHNTCOST, MHNTHNDL, MHNTFORCE, MHNTCONSQ,
    MHNTPRIV, MHNTPTHNK
  )

nsduh_2022_ss <- nsduh_2022_s |>
  filter(IRSUICTHNK == 1 | IRSUITRYYR == 1 | IRSUIPLANYR == 1)

# 2021
input_file <- file.path(input_dir, "NSDUH_2021_Tab.txt")
nsduh_2021 <- read_delim(input_file, show_col_types = FALSE)

nsduh_2021_s <- nsduh_2021 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, IRSUICTHNK, IRSUITRYYR, IRSUIPLANYR,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2021_ss <- nsduh_2021_s |>
  filter(IRSUICTHNK == 1 | IRSUITRYYR == 1 | IRSUIPLANYR == 1)

# 2020
input_file <- file.path(input_dir, "NSDUH_2020_Tab.txt")
nsduh_2020 <- read_delim(input_file, show_col_types = FALSE)

nsduh_2020_s <- nsduh_2020 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2020_ss <- nsduh_2020_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

# 2019
input_file <- file.path(input_dir, "NSDUH_2019_Tab.txt")
nsduh_2019 <- read_delim(input_file, show_col_types = FALSE)

nsduh_2019_s <- nsduh_2019 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2019_ss <- nsduh_2019_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

# 2018
input_file <- file.path(input_dir, "NSDUH_2018_Tab.tsv")
nsduh_2018 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

nsduh_2018_s <- nsduh_2018 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2018_ss <- nsduh_2018_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

# 2017
input_file <- file.path(input_dir, "NSDUH_2017_Tab.tsv")
nsduh_2017 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

nsduh_2017_s <- nsduh_2017 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2017_ss <- nsduh_2017_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

# 2016
input_file <- file.path(input_dir, "NSDUH_2016_Tab.tsv")
nsduh_2016 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

nsduh_2016_s <- nsduh_2016 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2016_ss <- nsduh_2016_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

# 2015
input_file <- file.path(input_dir, "NSDUH_2015_Tab.tsv")
nsduh_2015 <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

nsduh_2015_s <- nsduh_2015 |>
  filter(CATAG6 %in% c(2, 3, 4, 5)) |>
  select(
    CATAG6, IRSEX, NEWRACE2, EDUHIGHCAT, MHSUITHK, MHSUITRY, MHSUIPLN,
    MHRENUF2, MHRNBRS2, MHRTRAN2, MHRTIME2, MHRNCOV2, MHRWHER2,
    MHRNOHP2, MHRCOST2, MHRHAND2, MHRCMIT2, MHRJOBS2, MHRCFID2, MHRFOUT2
  )

nsduh_2015_ss <- nsduh_2015_s |>
  filter(MHSUITHK == 1 | MHSUITRY == 1 | MHSUIPLN == 1)

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

dt_list <- mapply(function(df, year) {
  dt <- as.data.table(df)
  dt[, YEAR := year]
  return(dt)
}, dataset_list, years, SIMPLIFY = FALSE)

combined_dt <- rbindlist(dt_list, use.names = FALSE)
colnames_2023 <- c(names(nsduh_2023_ss), "YEAR")
setnames(combined_dt, colnames_2023)
combined_nsduh_ss <- as.data.frame(combined_dt)

#------------------------------------------------------------------------------
# Combine 2015-2023 "S" data frames
#------------------------------------------------------------------------------
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

dt_list <- mapply(function(df, year) {
  dt <- as.data.table(df)
  dt[, YEAR := year]
  return(dt)
}, dataset_list, years, SIMPLIFY = FALSE)

combined_dt <- rbindlist(dt_list, use.names = FALSE)
colnames_2023 <- c(names(nsduh_2023_s), "YEAR")
setnames(combined_dt, colnames_2023)
combined_nsduh_s <- as.data.frame(combined_dt)

# Create an SI variable
combined_nsduh_si <- combined_nsduh_s |>
  mutate(
    SI = ifelse(
      IRSUICTHNK == 1 | IRSUITRYYR == 1 | IRSUIPLANYR == 1,
      1,
      0
    )
  )

combined_nsduh_svhb <- combined_nsduh_ss |>
  mutate(
    SI = ifelse(
      IRSUICTHNK == 1 | IRSUITRYYR == 1 | IRSUIPLANYR == 1,
      1,
      0
    )
  )

# Export CSV
export_file <- file.path(output_dir, "combined_nsduh_svhb.csv")
write.csv(combined_nsduh_svhb, file = export_file, row.names = FALSE)

# Summarize & create gtsummary table
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

output_html <- file.path(output_dir, "tab1.html")
gtsave(summary_table_gt, output_html)

# Create prevalence line graph
prevalence_sic <- combined_nsduh_si |>
  group_by(YEAR) |>
  summarize(
    total_number   = n(),
    si_occurrences = sum(SI == 1, na.rm = TRUE),
    Prevalence     = (si_occurrences / total_number) * 100
  ) |>
  select(Year = YEAR, Prevalence)

plot_si <- ggplot(prevalence_sic, aes(x = Year, y = Prevalence)) +
  geom_line(linewidth = 1, color = "#ffdd99") +
  geom_point(size = 2, color = "#ffdd99") +
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

output_file <- file.path(output_dir, "prevalence_si_plot.png")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = output_file,
  plot     = plot_si,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Create stacked bar chart
variables_of_interest <- c(
  "MHNTENFCV", "MHNTFFLKE", "MHNTPROBS", "MHNTTIME", "MHNTINSCV",
  "MHNTWHER", "MHNTNOHLP", "MHNTCOST", "MHNTHNDL", "MHNTFORCE",
  "MHNTCONSQ", "MHNTPRIV", "MHNTPTHNK"
)

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

data_filtered <- combined_nsduh_svhb |>
  filter(SI == 1) |>
  select(YEAR, all_of(variables_of_interest))

data_long <- data_filtered |>
  pivot_longer(
    cols      = -YEAR,
    names_to  = "Variable",
    values_to = "Count"
  )

summary_counts <- data_long |>
  filter(Count == 1) |>
  group_by(YEAR, Variable) |>
  summarise(Count = n(), .groups = "drop")

# Reverse the year factor so the stacked bars show 2023 at top, 2015 at bottom
summary_counts$YEAR <- factor(summary_counts$YEAR, levels = sort(unique(summary_counts$YEAR), decreasing = TRUE))

rea_si <- ggplot(summary_counts, aes(x = Variable, y = Count, fill = YEAR)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Reasons for Not Receiving MHS Among U.S. Adults with Suicidal Ideation",
    x = "Reasons",
    y = "Number of Cases",
    fill = "Year"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(breaks = seq(0, 6000, by = 500)) +
  scale_x_discrete(labels = variable_labels) +
  theme_bw() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "left",
    plot.title        = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.margin       = margin(t = 10, r = 10, b = 10, l = 5)
  ) +
  coord_flip()

output_file <- file.path(output_dir, "rea_si_plot.png")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = output_file,
  plot     = rea_si,
  width    = 8,
  height   = 6,
  dpi      = 300
)

message("Data processing and visualization complete!")

if (toupper(environment) != "LOCAL") {
  message(sprintf("ENVIRONMENT=%s => Uploading results to s3://%s/output/", environment, bucket))

  files_to_upload <- c(
    "combined_nsduh_svhb.csv",
    "tab1.html",
    "prevalence_si_plot.png",
    "rea_si_plot.png"
  )

  for (f in files_to_upload) {
    local_path <- file.path(output_dir, f)
    s3_object  <- paste0("output/", f)
    message(sprintf("  Uploading %s -> s3://%s/%s", local_path, bucket, s3_object))
    put_object(
      file   = local_path,
      object = s3_object,
      bucket = bucket,
      opts   = list(region = region)
    )
  }
  message("All output files uploaded to S3!")
} else {
  message("ENVIRONMENT=LOCAL => skipping S3 upload.")
}

message("Script complete!")

