## PART 1: Predict operons using ptt files for all 4 genomes
# set working directory & load in files
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW3")
e_coli <- read.table("E_coli.ptt", header = TRUE, sep = "\t")
b_subtilis <- read.delim("B_subtilis.ptt", header = TRUE, sep = "\t")
halobac <- read.delim("Halobacterium.ptt", header = TRUE, sep = "\t")
synechoc <- read.delim("Synechocystis.ptt", header = TRUE, sep = "\t")
# View(e_coli)
# View(b_subtilis)
# View(halobac)
# View(synechoc)

#######################################################
# Operons:  groups of genes that work together as a unit.
# genes belong in the same operon if they meet these criteria:
# They're adjacent to each other (next to one another on the genome)
# They're co-directional (pointing in the same direction — both reading left-to-right or both right-to-left)
# The gap between them is less than 50 base pairs
#################################################################
