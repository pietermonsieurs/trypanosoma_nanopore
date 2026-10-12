#!/usr/bin/env Rscript
#
# plot_porec_heatmap.R
#
# Plot a genome-wide Pore-C contact heatmap from the binned contact table
# produced by bin_porec_contacts.py -- a misassembly-QC check for a "final"
# assembly. Bins are kept in assembly/coordinate order (no re-clustering),
# with chromosome/scaffold boundaries drawn in, so the reading is the same
# as PretextView/Juicebox: a clean, continuous diagonal band within each
# sequence is consistent with the physical contact data; a broken diagonal,
# or an off-diagonal block linking two spots that shouldn't be near each
# other, flags a misjoin/mis-orientation candidate.
#
# Usage:
#   Rscript plot_porec_heatmap.R <prefix> [log1p|raw|oe] [out.png] [n_largest]
#
# scale = "oe" plots log2(observed/expected) instead of raw counts, with a
# blue-white-red diverging palette centered on 0 (= exactly as expected,
# i.e. no enrichment). "Expected" for cell (i,j) is what you'd see if all
# of bin i's and bin j's contacts were distributed with no proximity
# preference at all: rowSum(i) * rowSum(j) / total. This cancels out each
# bin's overall coverage/mappability level, which is what a raw or log1p
# map cannot do -- use this when the background trans level is so high
# (as here) that it swamps raw counts even after coarsening the bins.
#
# n_largest (optional, e.g. 11): keep only the N longest sequences in the
# assembly -- for this project, the 11 megabase chromosomes
# (Tbg972_01-11_RagTag), dropping the ~52 smaller unplaced scaffolds from
# both axes. Length is measured directly from the bins table (max bin end
# per contig), so nothing needs to be re-run upstream -- this only
# subsets the already-binned table. Omit it (or pass 0/NA) to keep every
# sequence, as before.
#
# Expects <prefix>.bins.tsv and <prefix>.contacts.tsv (from
# bin_porec_contacts.py) in the same directory/prefix.
#
# Requires: pheatmap, Matrix
#   install.packages(c("pheatmap", "Matrix"), repos = "https://cloud.r-project.org")
#
# Always prints two diagnostics before plotting (cheap, and the first thing
# to check if the map looks like a uniform wash instead of a clean diagonal):
#   - cis (within-sequence) vs. trans (between-sequence) contact fraction
#   - the 10 "hub" bins with the most total contacts, which is what a
#     repeat-driven multi-mapping artifact (VSG arrays, satellite repeats)
#     looks like -- a handful of bins linked to everything, everywhere.
#
# normalize (optional 5th arg, TRUE/FALSE, default FALSE): apply ICE-style
# iterative row/column balancing before plotting. Raw Hi-C/Pore-C counts are
# dominated by each bin's overall "visibility" (mappability/coverage), which
# is the usual reason a raw matrix looks like a uniform plaid rather than a
# clean diagonal -- this is the same normalization cooler/Juicer/HiCExplorer
# apply by default, done here by hand since we're not going through cooler.

suppressMessages({
  library(pheatmap)
  library(Matrix)
})

ice_normalize <- function(m, max_iter = 100, tol = 1e-4) {
  # Sinkhorn-Knopp / ICE-style balancing: rescale rows and columns
  # iteratively until every bin has (approximately) equal total contact
  # mass. Removes per-bin coverage/mappability bias without touching the
  # genuine proximity signal on the diagonal.
  mc <- m
  for (iter in seq_len(max_iter)) {
    rs <- rowSums(mc)
    nz <- rs > 0
    target <- mean(rs[nz])
    scale <- rep(1, length(rs))
    scale[nz] <- rs[nz] / target
    scale[scale == 0] <- 1
    mc <- sweep(mc, 1, scale, "/")
    mc <- sweep(mc, 2, scale, "/")
    if (max(abs(scale[nz] - 1)) < tol) break
  }
  mc
}

compute_oe <- function(m, pseudocount = 1) {
  # log2(observed/expected) under a no-proximity-preference null: expected
  # contact count for (i,j) is proportional to how much total contact mass
  # bins i and j carry overall, regardless of distance or chromosome. A
  # genuinely enriched cell (real 3D proximity, e.g. the true diagonal)
  # comes out positive; plain background trans should hover near 0.
  #
  # Pseudocount (standard Hi-C practice, e.g. cooltools): without it, a
  # bin with essentially zero real coverage (an assembly gap, a masked/
  # low-complexity patch) hits log2(0) = -Inf, which is indistinguishable
  # from -- and gets rendered the same extreme blue as -- a genuinely
  # depleted, well-covered bin. Adding it to both observed and expected
  # shrinks near-empty bins toward 0 (white = "no information") instead,
  # so a true low-coverage gap reads differently from a real hub/repeat
  # artifact (which has plenty of coverage, just spread thin -- that one
  # still shows up blue).
  rs <- rowSums(m)
  total <- sum(m)
  expected <- outer(rs, rs) / total
  oe <- (m + pseudocount) / (expected + pseudocount)
  oe[!is.finite(oe) | expected == 0] <- NA
  log2(oe)
}

args      <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_porec_heatmap.R <prefix> [log1p|raw|oe] [out.png] [n_largest] [normalize] [chrom]")
}
prefix    <- args[1]
scale     <- if (length(args) >= 2) args[2] else "log1p"
out_png   <- if (length(args) >= 3) args[3] else paste0(prefix, "_heatmap.png")
n_largest <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
if (!is.na(n_largest) && n_largest <= 0) n_largest <- NA_integer_
normalize <- if (length(args) >= 5) as.logical(args[5]) else FALSE
if (is.na(normalize)) normalize <- FALSE
# chrom (optional 6th arg, exact contig name e.g. Tbg972_05_RagTag): zoom
# the PLOT to just this chromosome's own bins, at full resolution with bp
# position labels. Zooming happens AFTER O/E is computed on the full
# genome-wide matrix (below), not by recomputing it on a subset -- doing
# that would change what "expected" means and make the values wrong.

# prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/porec_200k_bins/MM05_v2_porec"
# prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_v2_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_medaka_noragtag_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_medaka_ragtag_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_nomedaka_noragtag_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_nomedaka_ragtag_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/final_niki/MM05_medaka_noragtag.joined_porec"

prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/DAL972/TriTrypDB-68_TbruceigambienseDAL972_Genome_porec"
prefix = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_mm05/TREU927/TriTrypDB-67_TbruceiTREU927_Genome_nuclear_porec"



n_largest = 11
# n_largest = NA_integer_
normalize = FALSE
scale = "oe"

out_png   <- if (length(args) >= 3) args[3] else paste0(prefix, "_heatmap.png")


chrom     <- if (length(args) >= 6) args[6] else ""

bins     <- read.delim(paste0(prefix, ".bins.tsv"),     stringsAsFactors = FALSE)
contacts <- read.delim(paste0(prefix, ".contacts.tsv"), stringsAsFactors = FALSE)

# --- Optional: keep only the N longest sequences (e.g. the 11 chromosomes,
# dropping unplaced scaffolds from both axes of the heatmap) -------------
if (!is.na(n_largest)) {
  contig_len <- aggregate(end ~ contig, data = bins, FUN = max)
  contig_len <- contig_len[order(-contig_len$end), ]
  keep_contigs <- contig_len$contig[seq_len(min(n_largest, nrow(contig_len)))]

  keep_mask    <- bins$contig %in% keep_contigs
  keep_bin_ids <- bins$bin_id[keep_mask]                 # original 0-based ids
  remap        <- setNames(seq_along(keep_bin_ids) - 1L, # new 0-based ids,
                            as.character(keep_bin_ids))  # same relative order

  bins     <- bins[keep_mask, , drop = FALSE]
  contacts <- contacts[contacts$bin_i %in% keep_bin_ids &
                          contacts$bin_j %in% keep_bin_ids, , drop = FALSE]
  contacts$bin_i <- unname(remap[as.character(contacts$bin_i)])
  contacts$bin_j <- unname(remap[as.character(contacts$bin_j)])

  cat(sprintf("Filtered to %d largest sequences: %s\n",
              length(keep_contigs), paste(keep_contigs, collapse = ", ")))
}

n <- nrow(bins)

# Build the upper-triangular sparse matrix from the (bin_i, bin_j, count)
# table, then symmetrize (mirror below the diagonal) without double-adding
# the diagonal itself.
m <- sparseMatrix(
  i    = contacts$bin_i + 1L,
  j    = contacts$bin_j + 1L,
  x    = contacts$count,
  dims = c(n, n)
)
m <- m + t(m) - Diagonal(x = diag(m))
m <- as.matrix(m)

# --- Diagnostics: run BEFORE any normalization, on the raw counts -------
cis_mask   <- contacts$cis <- (bins$contig[contacts$bin_i + 1L] == bins$contig[contacts$bin_j + 1L])
cis_total  <- sum(contacts$count[cis_mask])
trans_total<- sum(contacts$count[!cis_mask])
cat(sprintf("cis (within-sequence): %d contacts (%.1f%%)\n",
            cis_total, 100 * cis_total / (cis_total + trans_total)))
cat(sprintf("trans (between-sequence): %d contacts (%.1f%%)\n",
            trans_total, 100 * trans_total / (cis_total + trans_total)))
cat("(a well-behaved map is heavily cis-dominated, typically >70-80%;",
    "a low cis fraction points to a mapping/specificity problem upstream",
    "of this plot, not something normalization here can fix)\n")

hub_order <- order(-rowSums(m))[seq_len(min(10, n))]
hub_tbl   <- data.frame(bins[hub_order, c("contig", "start", "end")],
                         total_contacts = rowSums(m)[hub_order])
cat("Top 10 hub bins (most total contacts -- cross-check against\n")
cat("08_repeats/main_assembly_repeats.tsv for repeat overlap):\n")
print(hub_tbl, row.names = FALSE)

if (normalize && scale != "oe") {
  # oe's own row/column-sum correction supersedes ICE balancing -- applying
  # both would double-correct the same coverage bias.
  cat("Applying ICE-style row/column balancing before plotting...\n")
  m <- ice_normalize(m)
} else if (normalize && scale == "oe") {
  cat("normalize=TRUE ignored for scale=oe (oe already corrects for per-bin coverage)\n")
}

m_plot <- switch(scale,
  "log1p" = log1p(m),
  "oe"    = compute_oe(m),
  m  # "raw" or anything else: plot as-is
)

if (scale == "oe" && !nzchar(chrom)) {
  # Genome-wide survey diagnostics -- skipped in zoom mode (chrom set),
  # since they're meant to scan across all chromosomes at once; the whole
  # point of zoom mode is to already know which chromosome you're
  # inspecting.
  # Two distinct candidate-artifact lists, matched to the two patterns
  # visible in an oe map: low-coverage bins (assembly gaps, masked
  # regions, or a genuinely foreign inserted stretch -- these are what
  # can produce a thin white/near-zero break cutting through an
  # otherwise-enriched diagonal), and uniformly-depleted "hub" bins
  # (repeat-driven ambiguous multi-mapping -- the solid blue stripes).
  raw_cov <- rowSums(m)
  low_cov_order <- order(raw_cov)[seq_len(min(15, n))]
  low_cov_tbl <- data.frame(bins[low_cov_order, c("contig", "start", "end")],
                             raw_contacts = raw_cov[low_cov_order])
  cat("\n15 lowest-coverage bins (candidates for the 'white line' pattern --\n")
  cat("assembly gap, masked/low-complexity patch, or a mis-placed small\n")
  cat("contig; check N-content and provenance -- e.g. was this stretch a\n")
  cat("TGS-GapCloser fill? -- before deciding to remove anything):\n")
  print(low_cov_tbl, row.names = FALSE)

  diag_idx <- seq_len(n)
  row_mean_oe <- sapply(diag_idx, function(i) mean(m_plot[i, -i], na.rm = TRUE))
  depleted_order <- order(row_mean_oe)[seq_len(min(15, n))]
  depleted_tbl <- data.frame(bins[depleted_order, c("contig", "start", "end")],
                              mean_log2_oe = round(row_mean_oe[depleted_order], 2),
                              raw_contacts = raw_cov[depleted_order])
  cat("\n15 most uniformly-depleted bins (candidates for the 'blue line'\n")
  cat("pattern -- repeat-driven hub bins; cross-check against\n")
  cat("08_repeats/main_assembly_repeats.tsv for repeat family overlap):\n")
  print(depleted_tbl, row.names = FALSE)
  cat("\n")
}

if (nzchar(chrom)) {
  chrom_mask <- bins$contig == chrom
  if (!any(chrom_mask)) {
    stop(sprintf("chrom '%s' not found. Contigs in %s.bins.tsv: %s",
                  chrom, prefix, paste(unique(bins$contig), collapse = ", ")))
  }
  m_plot <- m_plot[chrom_mask, chrom_mask, drop = FALSE]
  bins   <- bins[chrom_mask, , drop = FALSE]
  n      <- nrow(bins)
  cat(sprintf("Zoomed to %s: %d bins (%.2f - %.2f Mb)\n",
              chrom, n, min(bins$start) / 1e6, max(bins$end) / 1e6))
}

# Chromosome/scaffold boundaries in genome order (bins.tsv is already in
# .fai / assembly order) -- these become the grid lines pheatmap draws
# between blocks, same purpose as the thin lines in a Juicebox/Pretext view.
# In zoom mode there's only one contig left, so this naturally comes out
# as no internal boundaries -- correct, nothing to add here.
contig_rle <- rle(bins$contig)
boundaries <- cumsum(contig_rle$lengths)
boundaries <- boundaries[-length(boundaries)]

if (scale == "oe") {
  # Diverging palette centered on 0 (= exactly as expected). Clip
  # symmetrically at the 99th percentile of |log2(O/E)| so a few extreme
  # cells don't wash out contrast in the bulk, same rationale as the
  # 99th-percentile cap on the other scales, just symmetric here.
  cap <- as.numeric(quantile(abs(m_plot[upper.tri(m_plot)]), 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(abs(m_plot), na.rm = TRUE)
  breaks <- seq(-cap, cap, length.out = 101)
  palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
} else {
  # Clip the color scale at the 99th percentile rather than the raw max: a
  # handful of repeat-driven hub bins can otherwise stretch the whole ramp
  # and wash out the contrast between diagonal and background that you
  # actually came here to see.
  cap <- as.numeric(quantile(m_plot[upper.tri(m_plot)], 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(m_plot, na.rm = TRUE)
  breaks <- seq(0, cap, length.out = 101)
  palette <- colorRampPalette(c("white", "#b30000"))(100)
}

if (nzchar(chrom)) {
  # Zoomed to one chromosome: label bp position (Mb) at regular intervals
  # instead of the chromosome name -- position is what you need to read
  # off here, e.g. to find exactly where a break sits.
  label_every <- max(1, round(n / 20))
  label_idx   <- seq(1, n, by = label_every)
  axis_labels <- rep("", n)
  axis_labels[label_idx] <- sprintf("%.2fMb", bins$start[label_idx] / 1e6)
} else {
  # Label each chromosome/scaffold block at its midpoint bin only (not
  # every bin -- with ~15-20 bins per chromosome that would be unreadably
  # dense). Reading chromosome identity off raw grid position is exactly
  # the kind of thing that's easy to miscount by one in an 11x11 image --
  # this removes that ambiguity directly in the figure.
  mid_idx <- tapply(seq_len(n), bins$contig, function(idx) idx[ceiling(length(idx) / 2)])
  mid_idx <- mid_idx[unique(bins$contig)]  # keep genome order, not alphabetical
  axis_labels <- rep("", n)
  axis_labels[mid_idx] <- sub("_RagTag$", "", names(mid_idx))
}

png(out_png, width = 2400, height = 2400, res = 300)
pheatmap(
  m_plot,
  cluster_rows   = FALSE,   # keep genome/coordinate order -- do NOT re-cluster
  cluster_cols   = FALSE,   # (clustering would scramble the diagonal you're checking for)
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  labels_row     = axis_labels,
  labels_col     = axis_labels,
  fontsize_row   = 9,
  fontsize_col   = 9,
  angle_col      = "45",
  gaps_row       = boundaries,
  gaps_col       = boundaries,
  color          = palette,
  breaks         = breaks,
  na_col         = "grey80",
  main           = paste0("Pore-C contact map: ", basename(prefix),
                           if (scale == "oe") " (log2 obs/exp)" else ""),
  border_color   = NA
)
dev.off()

cat(sprintf("Wrote %s (%d x %d bins, %d sequences)\n",
            out_png, n, n, length(contig_rle$lengths)))


## focus only on scaffold_7
pheatmap(
  m_plot[1280:nrow(m_plot),1280:ncol(m_plot)],
  cluster_rows   = FALSE,   # keep genome/coordinate order -- do NOT re-cluster
  cluster_cols   = FALSE,   # (clustering would scramble the diagonal you're checking for)
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  # labels_row     = axis_labels,
  # labels_col     = axis_labels,
  # fontsize_row   = 9,
  # fontsize_col   = 9,
  # angle_col      = "45",
  # gaps_row       = boundaries,
  # gaps_col       = boundaries,
  color          = palette,
  breaks         = breaks,
  na_col         = "grey80",
  main           = paste0("Pore-C contact map: ", basename(prefix),
                          if (scale == "oe") " (log2 obs/exp)" else ""),
  border_color   = NA
)


