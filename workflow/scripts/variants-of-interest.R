#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(data.table)
library(glue)
library(openxlsx)

load_metadata <- function(metadata_path) {
  # Check the file extension and load metadata accordingly
  if (tools::file_ext(metadata_path) == "xlsx") {
    metadata <- readxl::read_excel(metadata_path)
  } else if (tools::file_ext(metadata_path) == "tsv") {
    metadata <- data.table::fread(metadata_path, sep = "\t")
  } else if (tools::file_ext(metadata_path) == "csv") {
    metadata <- data.table::fread(metadata_path, sep = ",")
  } else {
    stop("Metadata file must be .xlsx, .tsv, or .csv")
  }
  return(metadata)
}

sanitize_token <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "_", x)
  gsub("[^A-Za-z0-9._-]", "_", x)
}

make_mut_id <- function(names, locations) {
  raw_ids <- paste0(sanitize_token(names), "__", sanitize_token(locations))
  dup_idx <- ave(seq_along(raw_ids), raw_ids, FUN = seq_along)
  ifelse(dup_idx == 1, raw_ids, paste0(raw_ids, "__dup", dup_idx))
}

parse_location <- function(location) {
  location <- as.character(location)
  loc_parts <- strsplit(location, ":", fixed = TRUE)[[1]]
  chrom <- ifelse(length(loc_parts) >= 1, loc_parts[1], NA_character_)
  pos_block <- ifelse(length(loc_parts) >= 2, loc_parts[2], NA_character_)
  pos <- suppressWarnings(as.integer(strsplit(pos_block, "-", fixed = TRUE)[[1]][1]))
  list(chrom = chrom, pos = pos)
}

fallback_counts <- function(location) {
  loc <- parse_location(location)
  data.table(
    chrom = loc$chrom,
    pos = loc$pos,
    ref = "N",
    cov = 0,
    A = 0,
    C = 0,
    G = 0,
    T = 0
  )
}

safe_read_counts <- function(path, location) {
  if (!file.exists(path)) {
    warning(glue("Missing allele count file: {path}. Using zero-coverage fallback row."))
    return(fallback_counts(location))
  }

  dt <- tryCatch(
    fread(path, sep = "\t", fill = TRUE),
    error = function(e) data.table()
  )

  if (nrow(dt) == 0) {
    return(fallback_counts(location))
  }

  required_cols <- c("chrom", "pos", "ref", "cov", "A", "C", "G", "T")
  for (col_name in required_cols) {
    if (!col_name %in% names(dt)) {
      default_value <- if (col_name %in% c("chrom", "ref")) NA_character_ else 0
      dt[, (col_name) := default_value]
    }
  }

  dt <- dt[, ..required_cols]
  numeric_cols <- c("pos", "cov", "A", "C", "G", "T")
  dt[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
  dt
}

first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[[1]]
}

make_unique_sheet_names <- function(labels) {
  labels <- as.character(labels)
  labels <- gsub("[:\\\\/?*\\[\\]]", "_", labels)
  labels[is.na(labels) | labels == ""] <- "Mutation"
  output <- character(length(labels))
  seen <- integer(0)
  names(seen) <- character(0)

  for (i in seq_along(labels)) {
    base <- labels[[i]]
    if (!base %in% names(seen)) {
      seen[[base]] <- 1L
      output[[i]] <- substr(base, 1, 31)
    } else {
      seen[[base]] <- seen[[base]] + 1L
      suffix <- paste0("_", seen[[base]])
      keep <- max(1, 31 - nchar(suffix))
      output[[i]] <- paste0(substr(base, 1, keep), suffix)
    }
  }

  if (any(duplicated(output))) {
    output <- make.unique(output, sep = "_")
    output <- substr(output, 1, 31)
  }
  output
}

#### allele imbalance ####
metadata = load_metadata(snakemake@input[['metadata']])
samples = metadata$sampleID
# Read IR mutation data 
mutation_data = fread(snakemake@input[['mutations']], sep="\t")
mutation_data[, mutID := make_mut_id(Name, Location)]
mutation_data[, ALT2 := ifelse(is.na(ALT2), "", ALT2)]

all_list = vector("list", nrow(mutation_data))
mean_list = vector("list", nrow(mutation_data))
display_names = mutation_data$Name

# Loop through each mutation, finding allele coverage at each position
for (i in seq_len(nrow(mutation_data))){
  mut <- mutation_data[i]
  m <- mut$Name
  mut_id <- mut$mutID
  
  base <- toupper(mut$ALT)
  propstring <- glue("proportion{base}")
  base2 <- toupper(mut$ALT2)
  propstring2 <- glue("proportion{base2}")
  has_base <- base %in% c("A", "C", "G", "T")
  has_base2 <- nzchar(base2) && base2 %in% c("A", "C", "G", "T")
  
  #### load allele balance data ####
  allele_list = list()
  # read data for each sample and subset to what we want
  for (sample in samples){
    counts_path <- glue("results/variantAnalysis/variantsOfInterest/counts/{mut_id}/{sample}_allele_counts.tsv")
    allele_list[[sample]] <- safe_read_counts(counts_path, mut$Location)
    allele_list[[sample]]$sample = sample                                            #add sample column
    allele_list[[sample]]$treatment = metadata$treatment[samples == sample]         #add treatment column
    allele_list[[sample]]$mutation = m
    allele_list[[sample]]$gene = mut$Gene
    
    cover = allele_list[[sample]] %>% select(A,C,G,T) %>% rowSums()
    if (has_base) {
      allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring := ifelse(cover > 0, (!!sym(base))/cover, NA_real_))
    } else {
      warning(glue("Mutation '{m}' has invalid ALT '{base}'. Setting frequency to NA."))
      allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring := NA_real_)
    }
    
    if (has_base2){
      allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring2 := ifelse(cover > 0, (!!sym(base2))/cover, NA_real_))
    }
  }
  
  # We have 24 separate dataframes in a list, bind them together into one big dataframe
  alleles = rbindlist(allele_list, fill = TRUE)
  
  # now lets calculate population means and lower and upper CIs
  alleles_per_pop_list = list()
  for (pop in unique(metadata$treatment)){
    
    alleles_per_pop = alleles %>% filter(treatment == pop) 
    
    sum_cov <- sum(alleles_per_pop$cov, na.rm = TRUE)
    n_obs <- sum(alleles_per_pop$cov > 0, na.rm = TRUE)
    pop_prop <- if (has_base && sum_cov > 0) sum(alleles_per_pop[, ..base], na.rm = TRUE) / sum_cov else NA_real_
    error <- if (!is.na(pop_prop) && n_obs > 0) sqrt((pop_prop * (1 - pop_prop)) / n_obs) * 1.96 else NA_real_
    lower <- if (!is.na(error)) pmax(pop_prop - error, 0) else NA_real_
    upper <- if (!is.na(error)) pmin(pop_prop + error, 1) else NA_real_
    
    alleles_per_pop_list[[pop]] = alleles_per_pop %>% mutate(!!propstring := pop_prop, lowerCI = lower, upperCI = upper)
    
    if (has_base2){
      pop_prop2 <- if (sum_cov > 0) sum(alleles_per_pop[, ..base2], na.rm = TRUE) / sum_cov else NA_real_
      error2 <- if (!is.na(pop_prop2) && n_obs > 0) sqrt((pop_prop2 * (1 - pop_prop2)) / n_obs) * 1.96 else NA_real_
      lower2 <- if (!is.na(error2)) pmax(pop_prop2 - error2, 0) else NA_real_
      upper2 <- if (!is.na(error2)) pmin(pop_prop2 + error2, 1) else NA_real_
      
      alleles_per_pop_list[[pop]] = alleles_per_pop_list[[pop]] %>%
        mutate(!!propstring2 := pop_prop2, lowerCI_2 = lower2, upperCI_2 = upper2)
    }
  }
  
  mean_alleles = rbindlist(alleles_per_pop_list, fill=TRUE)
  
  if (has_base2){
    mean_alleles = mean_alleles %>%
      group_by(chrom, pos, ref, mutation, treatment, gene) %>%
      summarise(
        cov_mean = mean(cov, na.rm = TRUE),
        A_mean = mean(A, na.rm = TRUE),
        C_mean = mean(C, na.rm = TRUE),
        G_mean = mean(G, na.rm = TRUE),
        T_mean = mean(T, na.rm = TRUE),
        lowerCI = first_or_na(lowerCI),
        upperCI = first_or_na(upperCI),
        prop_value = first_or_na(!!sym(propstring)),
        prop_value_2 = first_or_na(!!sym(propstring2)),
        lowerCI_2 = first_or_na(lowerCI_2),
        upperCI_2 = first_or_na(upperCI_2),
        .groups = "drop"
      )
    mean_alleles <- mean_alleles %>% rename(!!propstring := prop_value, !!propstring2 := prop_value_2)
  } else {
    mean_alleles = mean_alleles %>%
      group_by(chrom, pos, ref, mutation, treatment, gene) %>%
      summarise(
        cov_mean = mean(cov, na.rm = TRUE),
        A_mean = mean(A, na.rm = TRUE),
        C_mean = mean(C, na.rm = TRUE),
        G_mean = mean(G, na.rm = TRUE),
        T_mean = mean(T, na.rm = TRUE),
        lowerCI = first_or_na(lowerCI),
        upperCI = first_or_na(upperCI),
        prop_value = first_or_na(!!sym(propstring)),
        .groups = "drop"
      )
    mean_alleles <- mean_alleles %>% rename(!!propstring := prop_value)
  }
  
  fwrite(alleles, glue("results/variantAnalysis/variantsOfInterest/csvs/{mut_id}_alleleBalance.csv"))
  fwrite(mean_alleles, glue("results/variantAnalysis/variantsOfInterest/csvs/mean_{mut_id}_alleleBalance.csv"))
  
  all_list[[i]] = alleles
  mean_list[[i]] = mean_alleles
}

#### write to excel file on diff sheets #### 
results_list = all_list
sheets = make_unique_sheet_names(display_names)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = TRUE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file=snakemake@output[['alleleBalance']], overwrite = TRUE)


### mean balance ####
#### write to excel file on diff sheets #### 
results_list = mean_list
sheets = make_unique_sheet_names(display_names)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = TRUE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file=snakemake@output[['mean_alleleBalance']], overwrite = TRUE)

sessionInfo()
