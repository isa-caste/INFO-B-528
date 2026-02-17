## Homework 1 README file

Title: Homework 1 <br>
Class: INFO B528, Computational Analysis of High-throughput Biomedical Data <br>
Author: Isabella Castellano <br>
Due Date: 2/6/26 <br>
Description: HW1.R is a R script that calculates the correlation of two matricies individally and then with eachother. It computes a pearson correlation for the matricies & creates heatmaps for vizualizations. <br>

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Language, tools & packages
### Language
- R v4.5.2
### Packages
- corrplot v0.95


Install the corrplot library in RStudio's global environment:
```r
install.packages("corrplot")
```
- The R package was installed globally in RStudio rather than in a project-specific library
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Inputs
- 'matrix1.txt' - first matrix of RPKM values showing the miRNA expression across cancers for a particular patient group as .txt file
- 'matrix2.txt' - second matrix of RPKM values showing the miRNA expression across cancers for a particular patient group as .txt file
### Files Location
- 'C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework?HW'
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Output 
### Files
- 'mat1_corr.csv' - correlation matrix from matrix 1
- 'mat2_corr.csv' - correlation matrix from matrix 2
- 'final_matrix.csv' - correlation matrix from matrix 1 & 2
- 'heatmap_mat1.png' - heatmap for from correlation matrix 1
- 'heatmap_mat2.png' - heatmap for from correlation matrix 1
- 'final_heatmap.png' - heatmap for from correlation matrix 1 & 2
### Files Location
- 'C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW1'
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## Generative AI Disclosure
Generative AI tools were used for debugging assistance & exporting files.