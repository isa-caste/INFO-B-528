# Homework 5 README file
Title: Homework 5<br>
Class: INFO B528, Computational Analysis of High-throughput Biomedical Data <br>
Author: Isabella Castellano <br>
Due Date: 4/10/26 <br>
Description: HW5.R is a R script that calculates the degree and average clustering coefficient of a given network and tests if the network is scale-free. <br>
Then the script calculates shortest path lengths and compares the distributions between two given protein lists with a wilcox test.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Language, tools & packages
### Language
- R v4.5.2
### Packages
- igraph v2.2.2

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Inputs
- 'protein-list1.txt' - txt file containing first list of proteins 
- 'protein-list2.txt' - txt file containing second list of proteins 
- 'Human-PPI.txt' - txt file containing human protein interaction network
### Files Location
- 'C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW5'
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Output 
<<<<<<< HEAD
- 'clustering_coefficients.csv' - csv file containing the local clustering coefficient for every node in the network
- 'Degree_Distribution.png' - scatter plot of the degree distribution of the network
- 'loglog_Degree_Distribution.png' - log-log scale degree distribution plot with power law fit line, gamma value, and R-squared
- 'Path_Length_Distribution.png' - histograms comparing shortest path length distributions for both protein lists
### Key results
- Wilcox Test Results:
=======
- "Degree_Distribution.png" - png containing plot of degree distribution of the human protein interaction network
- "loglog_Degree_Distribution.png" - png containing plot of log degree distribution of the human protein interaction network
- "Path_Length_Distribution.png" - png containing histogram comparing plot distributions of both protein lists
- Wilcox Test:
    - Test Statistic: 259326
    - P-value: 0.4879
    - 1176 pairs from protein list 1 and 435 pairs from protein list 2
    - This means that the distribution of lengths for paths for each protein pair in both lists are the same.
- Average Clustering Coefficient: 0.210
- Gamma: 1.548
- R-squared: 0.857
- The degree distribution fits a power law with gamma being 1.548 and R-squared value being 0.857. The high R-squared confirms a good power law fit, but the gamme value is suggests that the network has scale-free-like properties but is not strictly scale-free.
### Files Location
- 'C:/Users/icwin/OneDrive - Indiana University/INFO-B 528 High-TP/Homework/HW5'
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Generative AI Disclosure
Generative AI tools were used for debugging.
