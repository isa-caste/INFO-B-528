###### PART 2: Predict the operons in a Crop microbiome ######

# load library
library(dplyr)
library(stringr)
# load in file & add headers
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW3")
metagenome <- read.delim("metagenome.gff", header = TRUE, sep = "\t")
colnames(metagenome) <- c("Seq_ID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")

# create function for  detecting operons
operon_pred <- function(genome, filename) {
  # Parse Attributes column to extract locus_tag and product
  genome <- genome %>%
    mutate(
      Gene = str_extract(Attributes, "locus_tag=([^;]+)", group = 1),
      Product = str_extract(Attributes, "product=(.+)$", group = 1),
      Location = paste0(Start, "..", End)
    ) %>%
    select(Start, End, Strand, Gene, Product) 
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
    select(Operon_ID, Gene, Product)
  num_operons <- max(genome_with_operons$Operon_ID)
  # save results as csv
  write.csv(genome_with_operons, paste0(filename, "_operons.csv"), row.names = FALSE)
  print(paste("Number of operons found in", filename, ":", num_operons))
  # return the results
  return(genome_with_operons)
}

# run function for given file
operon_pred(metagenome, "Metagenome")