# set working directory & load in file
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW2")
data <- read.delim("DecayTimecourse.txt")
head(data)
# extract time points & gene names
times <- as.numeric(data[1, 2:10])
genes <- data[-1, 1]
# define columns for each timecourse
tc1_cols <- 2:10
tc2_cols <- 12:20
tc3_cols <- 20:28
# make function to calculate half life for one gene's expression data
calc_half_life <- function(expr, times) {
  expr <- as.numeric(expr)
  valid <- !is.na(expr) & expr > 0  
  if (sum(valid) < 3) return(NA)    
  fit <- lm(log(expr[valid]) ~ times[valid])
  k <- -coef(fit)[2]
  if (k <= 0) return(NA)    # only want to see when decay is positive
  return(log(2) / k)
}
# calculate half lives for each gene and timecourse
n_genes <- nrow(data) - 1  # exclude header row
hl1 <- hl2 <- hl3 <- numeric(n_genes)
for (i in 1:n_genes) {
  row <- i + 1  # first row 1 is time points
  hl1[i] <- calc_half_life(data[row, tc1_cols], times)
  hl2[i] <- calc_half_life(data[row, tc2_cols], times)
  hl3[i] <- calc_half_life(data[row, tc3_cols], times)
}
# get average for the three timecourses
avg_hl <- rowMeans(cbind(hl1, hl2, hl3), na.rm=TRUE)
# put all information in one dataframe
results <- data.frame(
  gene = genes,
  hl_tc1 = hl1,
  hl_tc2 = hl2,
  hl_tc3 = hl3,
  avg_half_life = avg_hl
)
# get rid of NAs for clean results
results_clean <- results[!is.na(results$avg_half_life), ]
top10 <- results_clean[results_clean$avg_half_life >= quantile(results_clean$avg_half_life, 0.9), ]
bottom10 <- results_clean[results_clean$avg_half_life <= quantile(results_clean$avg_half_life, 0.1), ]
# save the gene names
write.table(bottom10$gene, "bottom_10_gene_list.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(top10$gene, "top_10_gene_list.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
