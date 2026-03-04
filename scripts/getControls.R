library(tidyr)
library(lubridate)
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

smdb <- read_tsv(args[1], col_types = cols())  

datetime <- str_extract(basename(args[1]), "(?<=SMDB_).*(?=\\.tsv)")
name <- paste0("controls/control_", datetime, ".tsv")

# ----- Name controls 

# archive_sample_id or robot_sample_id should be smplntc (case insensitive) 

# this means sample negative control 

# if robot_sample_id and not SmplNTC,the first 6 chars can say 
# ExrNTC - extraction negative control
# ExrPTC - extraction positive control
# LibPTC - library positive control
# LibNTC - library negative control
# the numbers after are of the form YYMMDDNN where NN is a number that is incremented for each control sample on a given day

controls <- smdb %>%
  mutate(
    robot_sample_id = as.character(robot_sample_id),
    archive_sample_id = as.character(archive_sample_id),
    control_id = coalesce(robot_sample_id, archive_sample_id)
  ) %>%
  
  mutate(
    control_type = case_when(
      str_detect(robot_sample_id, regex("^ExrNTC", TRUE)) ~ "Extraction_Negative",
      str_detect(robot_sample_id, regex("^ExrPTC", TRUE)) ~ "Extraction_Positive",
      str_detect(robot_sample_id, regex("^LibNTC", TRUE)) ~ "Library_Negative",
      str_detect(robot_sample_id, regex("^LibPTC", TRUE)) ~ "Library_Positive",
      str_detect(control_id, regex("SmplNTC", TRUE)) ~ "Sample_Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(!is.na(control_type)) %>%
  
  mutate(
    # Only extract 8 digits for Exr/Lib controls
    yymmddnn = str_extract(
      control_id,
      "(?<=^(ExrNTC|ExrPTC|LibNTC|LibPTC))\\d{8}"
    ),
    
    # NA-safe date parsing
    extracted_date = lubridate::ymd(
      if_else(
        !is.na(yymmddnn),
        paste0("20", stringr::str_sub(yymmddnn, 1, 6)),
        NA_character_
      )
    ),
    
    # Fallback for SmplNTC
    control_date = coalesce(
      extracted_date,
      as.Date(robot_sample_sampling_date)
    )
  ) %>% select(library_id, control_id, control_type, control_date)

controls <- controls %>% filter(!is.na(library_id)) %>% distinct()


# ----- Find library results paths for controls -----

res_tsv <- system2(
  "scripts/findLibrary.sh",
  args = controls$library_id,
  stdout = TRUE
)
res_df <- read_tsv(
  paste(res_tsv, collapse = "\n"),
  show_col_types = FALSE
) %>% distinct()

df <- controls %>% left_join(res_df, by = c("library_id" = "library"))

df <- df %>% filter(!is.na(results_path) & results_wf != "ARCHIVE")

# ----- Get euk results  

df <- df %>%
  mutate(
    stat_path = file.path(
      results_path,
      "results",
      "metadmg",
      "aggregate",
      paste0("Lib_", library_id, "_collapsed.stat.gz")
    ),
    stat_exists = file.exists(stat_path)
  )

stat_data <- pmap_dfr(
  list(df$stat_path, df$library_id, df$results_path),
  function(stat_path, library_id, results_path) {
    
    read_tsv(
      stat_path,
      col_types = cols(.default = col_character()), 
      show_col_types = FALSE
    ) %>%
      mutate(
        library_id = library_id,
        results_path = results_path,
        stat_path = stat_path
      )
  }
)
stat_data$pipeline <- "EUKARYOTE"

# ----- Get prefilter results

df <- df %>%
  mutate(
    prefilter_path = file.path(
      results_path,
      "results",
      "prefilter_metadmg",
      "aggregate",
      paste0("Lib_", library_id, "_collapsed.stat.gz")
    ),
    stat_exists = file.exists(prefilter_path)
  )

prefilter_data <- pmap_dfr(
  list(df$prefilter_path, df$library_id, df$results_path),
  function(prefilter_path, library_id, results_path) {
    
    read_tsv(
      prefilter_path,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ) %>%
      mutate(
        library_id   = library_id,
        results_path = results_path,
        prefilter_path = prefilter_path
      )
  }
)

prefilter_data$pipeline <- "PREFILTER"

# ----- Join results with controls

final_data <- bind_rows(stat_data, prefilter_data)

df_final <- final_data %>% select(-prefilter_path, -stat_path) %>%
  left_join(df,by=c(library_id="library_id", results_path="results_path")) 

# remove columns called fw* or bw*
df_final <- df_final %>% select(-matches("^(fw|bw)")) 

df_final <- df_final %>% select(-nalign) %>%
  mutate(
    nreads = as.numeric(nreads),
    taxid  = as.character(taxid),
    A = as.numeric(A),
    mean_rlen = as.numeric(mean_rlen),
    mean_gc = as.numeric(mean_gc),
    date = as.Date(as.character(date),"%Y%m%d")
  )

write_tsv(df_final, name)


rmarkdown::render(
  "scripts/controlFreq.Rmd",
  params = list(
    input_file = name,
    datetime = datetime
  ),
  output_file = paste0("controlFreq_", datetime, ".html"),
  output_dir = "reports",
  knit_root_dir = getwd(),
  envir = new.env()
)