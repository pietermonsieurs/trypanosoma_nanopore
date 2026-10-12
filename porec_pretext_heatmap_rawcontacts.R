matrix_path <- '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/contactmap_mm05_matrix.tsv'
out_prefix  <- '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/contactmap_mm05'
fai_path <- '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/purged.fa.fai'

#!/usr/bin/env Rscript
#
# contig_contact_cluster.R
#
# Reads the contig x contig raw contact-count matrix produced by
# contig_contact_cluster.py (the "*_matrix.tsv" output) and clusters
# contigs by shared Pore-C contact frequency.
#
# Confirmed on this project (2026-09-16): the raw/log1p clustering below
# sorts contigs almost entirely by LENGTH (Spearman rho = 0.87 between
# heatmap position and contig length on the MM05 purged Flye set) rather
# than by genuine proximity - longer contigs simply pick up more total
# Pore-C hits with everything, which swamps the specific pairwise signal.
# The length-normalized section (part 3 below) corrects for this by
# dividing each cell by sqrt(len_i * len_j) before clustering/plotting -
# run that version FIRST if you have a .fai available; only fall back to
# the unnormalized raw/log1p versions if you don't.
#
# Produces, in order:
#   1. base R heatmap(), clustered on RAW counts,   colored on RAW counts
#   2. base R heatmap(), clustered on LOG1P counts, colored on LOG1P counts
#   3. pheatmap(),        clustered on RAW counts,   colored on RAW counts
#   4. pheatmap(),        clustered on LOG1P counts, colored on LOG1P counts
#   5. [only if a .fai is given] pheatmap(), clustered + colored on
#      LENGTH-NORMALIZED counts (raw normalized values)
#   6. [only if a .fai is given] pheatmap(), clustered + colored on
#      LOG1P of the length-normalized counts
#
# Requires the "pheatmap" package (install once: install.packages("pheatmap"),
# or on Debian/Ubuntu: apt-get install r-cran-pheatmap). Everything else is
# base R (stats/grDevices/graphics/utils).
#
# Usage:
#   Rscript contig_contact_cluster.R mm05_flye_matrix.tsv mm05_flye_R
#   Rscript contig_contact_cluster.R mm05_flye_matrix.tsv mm05_flye_R purged.fa.fai
#
# The .fai is the standard samtools faidx index of the SAME purged contig
# set the Pore-C reads were mapped against (samtools faidx purged.fa) -
# only the contig-name/length columns are used.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2, 3))) {
  stop("Usage: Rscript contig_contact_cluster.R <matrix.tsv> <out_prefix> [contigs.fai]")
}
matrix_path <- args[1]
out_prefix  <- args[2]
fai_path    <- if (length(args) == 3) args[3] else NA

library(pheatmap)

palette100 <- colorRampPalette(c("white", "orange", "red", "darkred"))(100)

to_dist <- function(x) {
  d <- 1 / (1 + x)
  diag(d) <- 0
  as.dist(d)
}

# ---------------------------------------------------------------
# 1. Load the contig x contig count matrix
#    (tab-separated, first column = row names, header = column names,
#    exactly as written by contig_contact_cluster.py's
#    pd.DataFrame(...).to_csv(..., sep="\t"))
# ---------------------------------------------------------------
cat(sprintf("[1] Reading %s ...\n", matrix_path))
df  <- read.delim(matrix_path, row.names = 1, check.names = FALSE)
mat <- as.matrix(df)
mode(mat) <- "numeric"

mat = mat[rownames(mat) != "ontig_42_1", colnames(mat) != "ontig_42_1"]  # remove the single contig that is a known artifact

n <- nrow(mat)
cat(sprintf("    %d contigs\n", n))

mat_log <- log1p(mat)  # log(1 + counts); log1p handles zeros cleanly, plain log(0) would not

# ---------------------------------------------------------------
# 2. Raw and log1p clustering + plots (unnormalized - kept for
#    comparison against the length-normalized version below; on this
#    project these two are dominated by contig length, see note above)
# ---------------------------------------------------------------
cat("[2] Clustering + plotting RAW and LOG1P (unnormalized) ...\n")

hc_raw <- hclust(to_dist(mat),     method = "average")
hc_log <- hclust(to_dist(mat_log), method = "average")

writeLines(rownames(mat)[hc_raw$order], con = paste0(out_prefix, "_clustered_order_raw.txt"))
writeLines(rownames(mat)[hc_log$order], con = paste0(out_prefix, "_clustered_order_log.txt"))

png(paste0(out_prefix, "_heatmap_raw_R.png"), width = 2000, height = 2000, res = 200)
heatmap(mat, Rowv = as.dendrogram(hc_raw), Colv = as.dendrogram(hc_raw), scale = "none",
        col = palette100, margins = c(8, 8), cexRow = 0.4, cexCol = 0.4,
        main = "Contig-level Pore-C contacts\n(raw counts, raw-based clustering - NOT length-normalized)")
dev.off()

png(paste0(out_prefix, "_heatmap_log_R.png"), width = 2000, height = 2000, res = 200)
heatmap(mat_log, Rowv = as.dendrogram(hc_log), Colv = as.dendrogram(hc_log), scale = "none",
        col = palette100, margins = c(8, 8), cexRow = 0.4, cexCol = 0.4,
        main = "Contig-level Pore-C contacts\n(log1p counts, log1p-based clustering - NOT length-normalized)")
dev.off()

pheatmap(mat, cluster_rows = hc_raw, cluster_cols = hc_raw, color = palette100,
         fontsize_row = 6, fontsize_col = 6,
         main = "Contig-level Pore-C contacts (raw, NOT length-normalized)",
         filename = paste0(out_prefix, "_pheatmap_raw.png"), width = 10, height = 10)

pheatmap(mat_log, cluster_rows = hc_log, cluster_cols = hc_log, color = palette100,
         fontsize_row = 6, fontsize_col = 6,
         main = "Contig-level Pore-C contacts (log1p, NOT length-normalized)",
         filename = paste0(out_prefix, "_pheatmap_log.png"), width = 10, height = 10)

# ---------------------------------------------------------------
# 3. Length-normalized version (only runs if a .fai was given)
#
#    Divides every cell by sqrt(len_i * len_j) - the geometric mean of
#    the two contigs' lengths - so two big contigs no longer look
#    "close" just because both are big. Rescaled by the median contig
#    length afterwards purely so the resulting numbers land in a
#    readable range; the rescaling constant does not affect clustering
#    or relative color contrast, only the numbers on the color-scale
#    legend.
# ---------------------------------------------------------------
fai <- read.delim(fai_path, header = FALSE, col.names = c("name","length","offset","linebases","linewidth"))
len <- setNames(fai$length, fai$name)

## filter out the contigs with small size (< 10000)
min_size <- 20000
small_contigs = names(len[len < min_size])

mat_filtered = mat[! rownames(mat) %in% small_contigs, ! colnames(mat) %in% small_contigs] 
len_filtered = len[! names(len) %in% small_contigs]

missing <- setdiff(rownames(mat_filtered), names(len_filtered))
if (length(missing) > 0) {
  stop("These contigs are in the matrix but not in the .fai (name mismatch?): ",
       paste(missing, collapse = ", "))
}
len_vec <- len[rownames(mat_filtered)]  # reorder to match the matrix

cat("    Length range in this matrix: ", format(min(len_vec), big.mark=","),
    " - ", format(max(len_vec), big.mark=","), " bp\n", sep = "")

denom <- sqrt(outer(len_vec, len_vec))
mat_norm <- (mat_filtered / denom) * median(len_vec)  # rescaled for readability only
mat_norm_log <- log1p(mat_norm)

cat("[3] Clustering + plotting LENGTH-NORMALIZED (raw and log1p) ...\n")

hc_norm     <- hclust(to_dist(mat_norm),     method = "average")
hc_norm_log <- hclust(to_dist(mat_norm_log), method = "average")

writeLines(rownames(mat)[hc_norm$order],     con = paste0(out_prefix, "_clustered_order_normalized.txt"))
writeLines(rownames(mat)[hc_norm_log$order], con = paste0(out_prefix, "_clustered_order_normalized_log.txt"))

pheatmap(mat_norm, cluster_rows = hc_norm, cluster_cols = hc_norm, color = palette100,
         fontsize_row = 6, fontsize_col = 6,
         main = "Contig-level Pore-C contacts\n(LENGTH-NORMALIZED counts, normalized-based clustering)",
         filename = paste0(out_prefix, "_pheatmap_normalized.png"), width = 10, height = 10)

pheatmap(mat_norm_log, cluster_rows = hc_norm_log, cluster_cols = hc_norm_log, color = palette100,
         fontsize_row = 6, fontsize_col = 6,
         main = "Contig-level Pore-C contacts\n(log1p of LENGTH-NORMALIZED counts, normalized-log-based clustering)",
         filename = paste0(out_prefix, "_pheatmap_normalized_log.png"), width = 10, height = 10)

# Quick sanity check: how strongly does the (still unnormalized) raw
# clustering order correlate with contig length? Printed for reference -
# a high value here reproduces the length-bias diagnosis; after
# normalizing, that correlation should drop substantially.
ord_raw  <- match(rownames(mat)[hc_raw$order],  names(len_vec))
ord_norm <- match(rownames(mat)[hc_norm$order], names(len_vec))
cat(sprintf("\n    Spearman(order position, contig length) - raw clustering:        %.3f\n",
            suppressWarnings(cor(seq_along(ord_raw),  len_vec[rownames(mat)[hc_raw$order]],  method = "spearman"))))
cat(sprintf(  "    Spearman(order position, contig length) - normalized clustering: %.3f\n",
              suppressWarnings(cor(seq_along(ord_norm), len_vec[rownames(mat)[hc_norm$order]], method = "spearman"))))
cat("    (normalized should be much closer to 0 than raw, if length bias was the issue)\n")


cat("\nDone. Files written with prefix:", out_prefix, "\n")

