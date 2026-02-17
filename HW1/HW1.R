# load in files
mat1 <- read.delim("C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW1/matrix1.txt")
mat2 <- read.delim("C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW1/matrix2.txt")
# Remove first &last column
Mat1 <- mat1 <- mat1[, -c(1, 14)]
Mat2 <- mat2 <- mat2[, -c(1, 14)]
# get correlation
cor_mat1 <- cor(Mat1, method = "pearson")
cor_mat2 <- cor(Mat2, method = "pearson")
# make heatmaps
# install.packages("corrplot")
library(corrplot)
png("heatmap_mat1.png", width = 800, height = 800, res = 120)
col <- colorRampPalette(c("blue", "white", "red"))(40)
heatmap(x = cor_mat1, 
        col = col, 
        symm = TRUE,
        margins = c(6, 6),
        cexRow = 1,         
        cexCol = 1,         
        main = "Cancer Correlation Matrix 1 Heatmap")
dev.off()
png("heatmap_mat2.png", width = 800, height = 800, res = 120)
heatmap(x = cor_mat2, 
        col = col, 
        symm = TRUE,
        margins = c(6, 6),  
        cexRow = 1,         
        cexCol = 1,         
        main = "Cancer Correlation Matrix 2 Heatmap")
dev.off()
# compute a pearson correlation between 2 matrices
cor_mat <- cor(t(cor_mat1),t(cor_mat2), method = "pearson")
png("final_heatmap.png", width = 800, height = 800, res = 120)
col <- colorRampPalette(c("blue", "white", "red"))(40)
heatmap(x = cor_mat, 
        col = col, 
        symm = TRUE,
        margins = c(6, 6),  
        cexRow = 1,         
        cexCol = 1,         
        main = "Cancer Correlation Matrices Heatmap")  
dev.off()
# export final matrix
write.csv(cor_mat, file = "C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW1/final_matrix.csv")