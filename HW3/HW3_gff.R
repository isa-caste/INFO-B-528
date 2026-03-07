## PART 2: Predict the operons in a Crop microbiome from Hoatzin
# load in file & add headers
metagenome <- read.delim("metagenome.gff", header = TRUE, sep = "\t")
colnames(metagenome) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# View(metagenome)

