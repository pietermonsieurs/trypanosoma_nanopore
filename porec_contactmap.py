#!/usr/bin/env python3
"""
contig_contact_cluster.py

Build a contig-level Pore-C contact matrix directly from the
`*_contacts.bed` file produced by porec_to_pairs.py (Step 5 of
02a_porec_scaffold_flye.slurm / 02b_porec_scaffold_hifiasm.slurm),
and reorder it by hierarchical clustering so that frequently-linked
contigs sit next to each other -> visible "blocks" without running
YaHS or building a full binned .hic/.mcool matrix.

This is a quick sanity-check / exploration tool, run on the purged
(purge_dups output) but NOT-yet-scaffolded contig set. It is NOT a
replacement for YaHS/HapHiC: it ignores enzyme cut-site density,
contig length normalisation subtleties, and orientation - it only
answers "which contigs cluster together", which is useful to eyeball
before committing to a full scaffolding run (e.g. does the number of
clusters look like ~11 chromosomes + some unplaced junk, or something
clearly wrong).

Usage:
    python contig_contact_cluster.py MM05_flye_contacts.bed out_prefix

Outputs:
    out_prefix_matrix.tsv         raw contig x contig contact counts
    out_prefix_heatmap.png        naive (alphabetical) order
    out_prefix_heatmap_clustered.png   reordered by clustering - this
                                        is the one where blocks show up
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform


def load_contact_counts(bed_path):
    """
    Parse a paired contacts BED (two lines per contact, sharing the
    same read+pair-index name in column 4) into a
    {(contig_a, contig_b): count} dict, contig_a < contig_b.
    """
    cols = ["chrom", "start", "end", "pair_id", "mapq", "strand"]
    df = pd.read_csv(bed_path, sep="\t", names=cols, dtype={"chrom": str, "pair_id": str})

    counts = defaultdict(int)
    n_skipped_malformed = 0
    n_intra = 0

    for pair_id, grp in df.groupby("pair_id", sort=False):
        contigs = grp["chrom"].tolist()
        if len(contigs) != 2:
            # a well-formed paired BED should always have exactly 2 lines
            # per pair_id; anything else means the file was filtered/
            # truncated/concatenated oddly - skip rather than guess
            n_skipped_malformed += 1
            continue
        a, b = sorted(contigs)
        if a == b:
            # shouldn't normally appear (--intra-only-skip in the
            # upstream script), but skip defensively if it does -
            # not informative for a contig-level matrix
            n_intra += 1
            continue
        counts[(a, b)] += 1

    if n_skipped_malformed:
        print(f"  note: skipped {n_skipped_malformed} pair_id groups with != 2 lines", file=sys.stderr)
    if n_intra:
        print(f"  note: skipped {n_intra} intra-contig pairs", file=sys.stderr)

    return counts


def build_matrix(counts):
    contigs = sorted({c for pair in counts for c in pair})
    idx = {c: i for i, c in enumerate(contigs)}
    n = len(contigs)
    mat = np.zeros((n, n), dtype=float)
    for (a, b), c in counts.items():
        i, j = idx[a], idx[b]
        mat[i, j] = c
        mat[j, i] = c
    return contigs, mat


def cluster_order(mat):
    """
    Hierarchical clustering on a 'distance' derived from contact counts:
    more contacts between two contigs -> smaller distance -> more likely
    to end up adjacent after reordering. This is the step that actually
    produces the visible blocks.
    """
    with np.errstate(divide="ignore"):
        dist = 1.0 / (1.0 + mat)
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    return leaves_list(Z)


def plot_heatmap(mat, labels, out_png, title, fontsize=4):
    n = len(labels)
    plt.figure(figsize=(max(8, n * 0.12), max(8, n * 0.12)))
    plt.imshow(np.log1p(mat), cmap="YlOrRd")
    plt.xticks(range(n), labels, rotation=90, fontsize=fontsize)
    plt.yticks(range(n), labels, fontsize=fontsize)
    plt.colorbar(label="log(1 + contacts)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <contacts.bed> <out_prefix>", file=sys.stderr)
        sys.exit(1)

    bed_path, out_prefix = sys.argv[1], sys.argv[2]

    print(f"[1/4] Parsing {bed_path} ...", file=sys.stderr)
    counts = load_contact_counts(bed_path)
    print(f"      {len(counts)} distinct contig pairs with >=1 contact", file=sys.stderr)

    print("[2/4] Building contig x contig matrix ...", file=sys.stderr)
    contigs, mat = build_matrix(counts)
    n = len(contigs)
    print(f"      {n} contigs involved in at least one inter-contig contact", file=sys.stderr)

    # raw counts table, useful on its own for sanity-checking specific pairs
    pd.DataFrame(mat, index=contigs, columns=contigs).to_csv(f"{out_prefix}_matrix.tsv", sep="\t")

    print("[3/4] Plotting naive (alphabetical) order heatmap ...", file=sys.stderr)
    plot_heatmap(mat, contigs, f"{out_prefix}_heatmap.png",
                 "Contig-level Pore-C contacts (alphabetical order - no structure expected)")

    print("[4/4] Clustering and plotting reordered heatmap ...", file=sys.stderr)
    order = cluster_order(mat)
    mat_ordered = mat[np.ix_(order, order)]
    contigs_ordered = [contigs[i] for i in order]
    plot_heatmap(mat_ordered, contigs_ordered, f"{out_prefix}_heatmap_clustered.png",
                 "Contig-level Pore-C contacts (clustered order - blocks = candidate groups)")

    print("Done.", file=sys.stderr)
    print(f"  {out_prefix}_matrix.tsv", file=sys.stderr)
    print(f"  {out_prefix}_heatmap.png", file=sys.stderr)
    print(f"  {out_prefix}_heatmap_clustered.png   <- look at this one", file=sys.stderr)


if __name__ == "__main__":
    main()