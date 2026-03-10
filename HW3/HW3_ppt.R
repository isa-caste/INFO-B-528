###### PART 1: Predict operons using ptt files for 4 genomes ######

# load library
library(dplyr)
# set working directory & load in files
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW3")
e_coli <- read.table("E_coli.ptt", header = TRUE, sep = "\t")
b_subtilis <- read.delim("B_subtilis.ptt", header = TRUE, sep = "\t")
halobac <- read.delim("Halobacterium.ptt", header = TRUE, sep = "\t")
synechoc <- read.delim("Synechocystis.ptt", header = TRUE, sep = "\t")

# create function for detecting operons
operon_pred <- function(genome, filename) {
  # split location column to get start & end position for gene
  split_data <- strsplit(genome$Location, "\\.\\.")
  matrix_data <- do.call(rbind, split_data)
  genome$Start <- matrix_data[,1]
  genome$End <- matrix_data[,2]
  # drop location column
  genome$Location <- NULL
  # select necessary columns for prediction
  genome <- genome[, c("Start", "End", "Strand", "Gene", "Product")]
  # compare adjacent genes 
  genome_with_operons <- genome %>%
    arrange(Start) %>%
    mutate(
      Start = as.numeric(Start),
      End = as.numeric(End),
      next_Start = lead(Start),
      next_Strand = lead(Strand),
      same_operon = (Strand == next_Strand) & ((next_Start - End) < 50),
      same_operon = ifelse(is.na(same_operon), FALSE, same_operon),
      operon_group = cumsum(!same_operon),
      Operon_ID = dense_rank(operon_group)
    ) %>%
    group_by(Operon_ID) %>%
    filter(n() > 1) %>%  # Keep only operons with 2+ genes
    ungroup() %>%
    select(Operon_ID, Gene, Product)
  num_operons <- max(genome_with_operons$Operon_ID)
  # save results as csv
  write.csv(genome_with_operons, paste0(filename, "_operons.csv"), row.names = FALSE)
  print(paste("Number of operons found in", filename, ":", num_operons))
  # return the results
  return(genome_with_operons)
}

# run function for all 4 files
operon_pred(e_coli, "E_coli")
operon_pred(b_subtilis, "B_subtilis")
operon_pred(halobac, "Halobacterium")
operon_pred(synechoc, "Synechocystis")

