##### Part 1 #####
# install libraries & load in file
# install.packages("igraph")
library(igraph)
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW5")
human_ppi <- read.table("Human-PPI.txt", header = TRUE, sep = "\t", row.names = NULL)
human_ppi <- na.omit(human_ppi)
# convert it to graph object
ppi_graph <- graph_from_data_frame(human_ppi, directed = FALSE)
# calculate degree for every node
node_degrees <- degree(ppi_graph)
degree_df <- data.frame(
  protein = names(node_degrees),
  degree = as.numeric(node_degrees))
# calculate clustering coeff for each node
node_cc <- transitivity(ppi_graph, type = "local")
cc_df <- data.frame(
  protein = V(ppi_graph)$name,
  clustering_coefficient = node_cc)
print(cc_df)                         
write.csv(cc_df, "clustering_coefficients.csv", row.names = FALSE)
# get average of clustering coeff
avg_cc <- mean(node_cc, na.rm = TRUE)
print(paste("Average Clustering Coefficient of the network:", round(avg_cc, 4)))
# test for scale-free structure
degree_counts <- table(node_degrees)
degree_vals <- as.numeric(names(degree_counts))
counts <- as.numeric(degree_counts)
# normal plot 
png(filename = "Degree_Distribution.png", width = 800, height = 600, res = 100)
plot(degree_vals, counts,
     xlab = "Degree", ylab = "Number of Nodes",
     main = "Degree Distribution",
     pch = 16)
dev.off()
# fit model to log-log data
fit <- lm(log10(counts) ~ log10(degree_vals))         
# extract gamma and r-squared value
gamma <- -coef(fit)[2]                                  
r_squared <- summary(fit)$r.squared                     
print(paste("Power-law exponent (gamma):", round(gamma, 3)))       
print(paste("R-squared:", round(r_squared, 3)))  
# get conclusion if network is scale free or not
if (gamma >= 2 & gamma <= 3) {                                      
  print("Network is scale-free (gamma between 2 and 3)")
} else {                                                            
  print(paste("Network is not be scale-free") )
}                                                                   
# plot log-log degree distribution with fitted line
png(filename = "loglog_Degree_Distribution.png", width = 800, height = 600, res = 100)
plot(log10(degree_vals), log10(counts),
     xlab = "log10(Degree)", ylab = "log10(Count)",
     main = "Degree Distribution (Log-Log Scale)",
     pch = 16)
abline(fit, col = "purple")                            
legend("topright", legend = paste("gamma =", round(gamma, 3), "\nR² =", round(r_squared, 3)),
       bty = "n")                                     
dev.off()
##### Part 2 #####
# load protein lists & filter
prot_1 <- read.table("protein-list1.txt")
prot_2 <- read.table("protein-list2.txt")
# make first column character
proteins1 <- as.character(prot_1$V1)
proteins2 <- as.character(prot_2$V1)
# filter & see how many passed through
proteins1_in_network <- proteins1[proteins1 %in% V(ppi_graph)$name]
proteins2_in_network <- proteins2[proteins2 %in% V(ppi_graph)$name]
length(proteins1_in_network)
length(proteins2_in_network)
# calculate all pairwise shortest paths
paths1 <- distances(ppi_graph, v = proteins1_in_network, to = proteins1_in_network)
paths2 <- distances(ppi_graph, v = proteins2_in_network, to = proteins2_in_network)
# take only unique pairs
paths1_vec <- paths1[upper.tri(paths1)]
paths2_vec <- paths2[upper.tri(paths2)]
# remove infinite values 
paths1_vec <- paths1_vec[is.finite(paths1_vec)]
paths2_vec <- paths2_vec[is.finite(paths2_vec)]
# see how many pairs there are 
length(paths1_vec)
length(paths2_vec)
# plot distributions
png(filename = "Path_Length_Distribution.png", width = 800, height = 600, res = 100)
par(mfrow = c(1, 2))
hist(paths1_vec,
     main = "Path Lengths - List 1",
     xlab = "Shortest Path Length",
     col = "purple")
hist(paths2_vec,
     main = "Path Lengths - List 2",
     xlab = "Shortest Path Length",
     col = "blue")
dev.off()
# run wilcox test
wilcox_result <- wilcox.test(paths1_vec, paths2_vec)
print(wilcox_result)
if (wilcox_result$p.value < 0.05) {
  print("Conclusion: The two protein sets have significantly different path length distributions (p < 0.05)")
} else {
  print(paste("Conclusion: No significant difference in path length distributions (p =", 
              round(wilcox_result$p.value, 4), ")"))}