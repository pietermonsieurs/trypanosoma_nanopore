#!/usr/bin/env python3
"""
bin_porec_contacts.py

Bin Pore-C contacts (from a name-sorted BAM of Pore-C reads aligned to an
assembly) into a genome-wide bin x bin contact matrix, for a whole-assembly
misassembly-QC heatmap -- the same underlying signal PretextMap/PretextView
show, aggregated into a table that plots easily in R.

Logic: group alignments by read name (works because the BAM is name-sorted).
A Pore-C read with N monomer alignments (primary + supplementary) contributes
N*(N-1)/2 pairwise contacts -- same decomposition this project's
porec_to_pairs.py uses for YaHS input. Each contact is assigned to a pair of
fixed-size genomic bins and tallied.

Usage:
    python3 bin_porec_contacts.py \
        --bam porec_vs_final.namesorted.bam \
        --fai MM05_v2.fasta.fai \
        --bin-size 20000 \
        --out-prefix 19_porec_qc_final/MM05_v2_porec

Outputs:
    <out-prefix>.bins.tsv      bin_id  contig  start  end
    <out-prefix>.contacts.tsv  bin_i   bin_j   count   (bin_i <= bin_j)

Requires: pysam (pip install --user pysam)
"""
import argparse
import itertools
import statistics
import sys
from collections import defaultdict

import pysam


def build_bins(fai_path, bin_size):
    """Read a .fai and lay out fixed-size bins in .fai (i.e. assembly) order.

    Using .fai order rather than re-sorting means the resulting axis matches
    however the final assembly's sequences are already ordered -- megabase
    chromosomes first, then unplaced scaffolds, matching MM05_v2.fasta.
    """
    bins = []               # list of (contig, start, end)
    bin_offset = {}         # contig -> index of its first global bin
    with open(fai_path) as fh:
        for line in fh:
            fields = line.split()
            contig, length = fields[0], int(fields[1])
            bin_offset[contig] = len(bins)
            n_bins = (length + bin_size - 1) // bin_size
            for i in range(n_bins):
                start = i * bin_size
                end = min(start + bin_size, length)
                bins.append((contig, start, end))
    return bins, bin_offset


def global_bin(contig, pos, bin_offset, bin_size):
    return bin_offset[contig] + (pos // bin_size)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True,
                     help="name-sorted BAM, Pore-C reads vs assembly (samtools sort -n)")
    ap.add_argument("--fai", required=True,
                     help="assembly .fai (samtools faidx <assembly.fasta>)")
    ap.add_argument("--bin-size", type=int, default=20000,
                     help="genomic bin size in bp (default 20000)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-mapq", type=int, default=20,
                     help="drop alignments below this MAPQ (matches the YaHS -q 20 "
                          "convention used elsewhere in this pipeline; default 20)")
    ap.add_argument("--min-aln-len", type=int, default=0,
                     help="drop alignments with aligned (query) length below this many bp "
                          "(default 0 = off). A short monomer alignment is more likely to "
                          "be an ambiguous/repeat-driven multi-mapper than a real, "
                          "confidently-placed Pore-C fragment -- try e.g. 150-300.")
    ap.add_argument("--max-monomers", type=int, default=250,
                     help="drop reads with more than this many filtered alignments "
                          "(default 250, matching wf-pore-c's own --max_monomers "
                          "default). A read with a runaway monomer count contributes "
                          "N*(N-1)/2 contacts -- a handful of such reads can dominate "
                          "the whole contact table, disproportionately as trans, if "
                          "the extra splits are alignment noise rather than real "
                          "ligation junctions (a real risk here: this pipeline has no "
                          "enzyme-aware virtual digestion step to validate a split "
                          "before counting it -- see pore-c-methods-overview.md).")
    args = ap.parse_args()

    bins, bin_offset = build_bins(args.fai, args.bin_size)

    with open(f"{args.out_prefix}.bins.tsv", "w") as fh:
        fh.write("bin_id\tcontig\tstart\tend\n")
        for i, (contig, start, end) in enumerate(bins):
            fh.write(f"{i}\t{contig}\t{start}\t{end}\n")

    counts = defaultdict(int)
    bam = pysam.AlignmentFile(args.bam, "rb")

    n_reads = 0
    n_multi = 0
    n_contacts = 0
    n_dropped_runaway = 0
    dropped_contacts_est = 0
    monomer_counts = []
    cis_kept, trans_kept = 0, 0
    cis_dropped_est, trans_dropped_est = 0, 0

    for qname, group in itertools.groupby(bam, key=lambda r: r.query_name):
        alns = [r for r in group
                if not r.is_unmapped
                and r.mapping_quality >= args.min_mapq
                and (args.min_aln_len == 0
                     or (r.query_alignment_length or 0) >= args.min_aln_len)]
        n_reads += 1
        if len(alns) < 2:
            continue
        n_multi += 1
        monomer_counts.append(len(alns))

        if len(alns) > args.max_monomers:
            n_dropped_runaway += 1
            for a, b in itertools.combinations(alns, 2):
                dropped_contacts_est += 1
                if a.reference_name == b.reference_name:
                    cis_dropped_est += 1
                else:
                    trans_dropped_est += 1
            continue

        # N monomer alignments -> N*(N-1)/2 pairwise contacts (same
        # decomposition as porec_to_pairs.py)
        for a, b in itertools.combinations(alns, 2):
            if a.reference_name == b.reference_name:
                cis_kept += 1
            else:
                trans_kept += 1
            gi = global_bin(a.reference_name, a.reference_start, bin_offset, args.bin_size)
            gj = global_bin(b.reference_name, b.reference_start, bin_offset, args.bin_size)
            if gi > gj:
                gi, gj = gj, gi
            counts[(gi, gj)] += 1
            n_contacts += 1

    with open(f"{args.out_prefix}.contacts.tsv", "w") as fh:
        fh.write("bin_i\tbin_j\tcount\n")
        for (i, j), c in sorted(counts.items()):
            fh.write(f"{i}\t{j}\t{c}\n")

    sys.stderr.write(
        f"[bin_porec_contacts] {n_reads} reads seen, {n_multi} with >=2 "
        f"filtered alignments (MAPQ>={args.min_mapq}), {n_contacts} pairwise "
        f"contacts binned into a {len(bins)} x {len(bins)} matrix "
        f"({args.bin_size} bp bins)\n"
    )

    if monomer_counts:
        sm = sorted(monomer_counts)
        def pct(p):
            return sm[min(len(sm) - 1, int(p / 100 * len(sm)))]
        sys.stderr.write(
            f"[bin_porec_contacts] monomers/read (reads with >=2 filtered "
            f"alignments): median={statistics.median(sm)}, p95={pct(95)}, "
            f"p99={pct(99)}, max={sm[-1]}\n"
        )

    sys.stderr.write(
        f"[bin_porec_contacts] cis/trans among KEPT contacts: "
        f"cis={cis_kept} ({100*cis_kept/max(1,cis_kept+trans_kept):.1f}%), "
        f"trans={trans_kept} ({100*trans_kept/max(1,cis_kept+trans_kept):.1f}%)\n"
    )

    sys.stderr.write(
        f"[bin_porec_contacts] dropped {n_dropped_runaway} reads with "
        f">{args.max_monomers} monomers (--max-monomers); these would have "
        f"contributed an estimated {dropped_contacts_est} contacts if kept "
        f"(cis={cis_dropped_est}, trans={trans_dropped_est}) -- compare "
        f"against the {n_contacts} contacts actually counted above. If "
        f"dropped_contacts_est is large and heavily trans, a small number "
        f"of runaway reads were the main driver of the low cis fraction.\n"
    )


if __name__ == "__main__":
    main()

    