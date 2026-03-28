##### Part 1: compute the frequency matrix ##### 
# load counts matrix
setwd("C:/Users/icwin/OneDrive - Indiana University/Documents/INFO-B-528/HW4")
count_mat <- read.table("counts-matrix.txt", header = FALSE)
# remove uneeded columns
count_mat<- count_mat[-c(1,2)]
View(count_mat)
# define variables
N_SITES   <- 27      # real binding sites
PSEUDO    <- 1       # pseudocount added to each base at each position
N_PSEUDO  <- N_SITES + (PSEUDO * nrow(count_mat))
BG_FREQ   <- 0.25    # uniform background frequency
MOTIF_LEN <- ncol(count_mat)  
# compute frequency matrix
freq <- count_mat / N_SITES
print(round(freq, 4))
# compute adjusted frequency matrix 
freq_adj <- (count_mat + PSEUDO) / N_PSEUDO
print(round(freq_adj, 4))
# compute weight Mmtrix
pwm <- log2(freq_adj / BG_FREQ)
print(round(pwm, 4))

#### Part 2: Use matrix to search for genes ##### 
# load in fasta file - needed custom parser
parse_fasta <- function(filepath) {
  lines <- suppressWarnings(readLines(filepath))
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  parts    <- strsplit(lines, "\\s+")
  gene_ids <- sapply(parts, `[`, 1)
  seqs     <- sapply(parts, `[`, 3)       
  data.frame(gene_id = gene_ids, sequence = tolower(seqs), stringsAsFactors = FALSE)
}
fasta_df <- parse_fasta("E_coli.400_50")
# compute reverse complement of sequence string
rev_comp <- function(seq) {
  comp <- c(a = "t", t = "a", c = "g", g = "c", n = "n")
  bases <- strsplit(seq, "")[[1]]
  paste(rev(comp[bases]), collapse = "")}
# create score for a single window 
score_window <- function(bases) {
  total <- 0
  for (j in seq_along(bases)) {
    b <- bases[j]
    if (!(b %in% rownames(pwm))) return(-Inf)   # ambiguous base
    total <- total + pwm[b, j]
  }
  total
}
# scan sequences
scan_sequences <- function(fasta_df) {
  results <- vector("list", nrow(fasta_df))
  for (i in seq_len(nrow(fasta_df))) {
    gene_id <- fasta_df$gene_id[i]
    seq     <- fasta_df$sequence[i]
    rc_seq  <- rev_comp(seq)
    seq_len <- nchar(seq)
    best_score  <- -Inf
    # scan forward strand
    if (seq_len >= MOTIF_LEN) {
      for (start in 1:(seq_len - MOTIF_LEN + 1)) {
        window <- strsplit(substring(seq, start, start + MOTIF_LEN - 1), "")[[1]]
        s      <- score_window(window)
        if (s > best_score) best_score <- s
      }
    }
    # scan reverse complement strand
    if (nchar(rc_seq) >= MOTIF_LEN) {
      for (start in 1:(nchar(rc_seq) - MOTIF_LEN + 1)) {
        window <- strsplit(substring(rc_seq, start, start + MOTIF_LEN - 1), "")[[1]]
        s      <- score_window(window)
        if (s > best_score) best_score <- s
      }
    }
    # create results
    results[[i]] <- data.frame(gene_id = gene_id, score = best_score,
                               stringsAsFactors = FALSE)
  }
  do.call(rbind, results)
}
gene_scores <- scan_sequences(fasta_df)
# sort by score, take top 30, add rank
top30 <- head(gene_scores[order(-gene_scores$score), ], 30)
top30 <- data.frame(Rank = 1:nrow(top30), gene_id = top30$gene_id)
# write results to txt file
write.table(top30, "top30genes.txt", row.names = FALSE, quote = FALSE, append = FALSE)