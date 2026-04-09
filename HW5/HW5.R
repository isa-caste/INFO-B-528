##### Part 1 #####
- use libraries for computing degree
- compute avg closeting coeff
- see if its scale free
- if # nodes with given degree vs degree itself has a slope of gamma where gamma is between 2-3
- do log and log on both lines and should be a straight line 
^ this is for whole network

setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW5")
prot_list_1 <- read.table("protein-list1.txt")
prot_list_2 <- read.table("protein-list2.txt")
human_ppi <- read.table("Human-PPI.txt", header = TRUE)
View(prot_list_2)
View(prot_list_1)
View(human_ppi)

##### Part 2 #####
- shortest path length btwn all nodes
^ for protein lists