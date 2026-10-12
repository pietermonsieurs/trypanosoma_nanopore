#!/usr/bin/env python3
"""
porec_to_pairs.py

Convert Pore-C alignments (BED, sorted by read name) into pairwise
contacts for Hi-C scaffolding tools such as YaHS.

WHY THIS IS NEEDED
------------------
Hi-C produces read PAIRS: two ends, one contact per read pair.
Pore-C produces multi-way CONCATEMERS: a single ONT read can contain
3, 5, 10+ ligated fragments from different parts of the genome.

Scaffolders like YaHS expect pairwise contacts. So each Pore-C read
with N mapped fragments is decomposed into its N*(N-1)/2 pairwise
combinations. This is the standard "pairwise decomposition" approach
(the same idea used by the pore_c / pore-c-py toolkits).

INPUT
-----
BED file from `bamToBed`, sorted by read name (column 4), e.g.:
    contig_1  1000  1150  read_abc  60  +
    contig_7  8000  8140  read_abc  60  -
    contig_7  9200  9350  read_abc  60  +

Fragments belonging to one Pore-C read share the same name, which is
why supplementary alignments must be KEPT when filtering the BAM.

OUTPUT (default: BED)  <- recommended for Pore-C
-------------------------------------------------
Paired BED, two consecutive lines per contact, sharing a base name with
/1 and /2 suffixes (the suffixes are optional in the YaHS spec, but make
the pairing explicit):

    contig_1  1000  1150  read_abc_p1/1  60  +
    contig_7  8000  8140  read_abc_p1/2  60  -

BED is preferred because it carries the real alignment start and end.

OUTPUT (--format pa5)
---------------------
Seven-column pairwise contact format:
    read_name  contig1  pos1  contig2  pos2  mapq1  mapq2

PA5 stores only a single position per fragment, so YaHS infers the
alignment span from --read-length (default 150), which is wrong for
Pore-C fragments. Use BED unless you have a reason not to.

USAGE
-----
    python3 porec_to_pairs.py in.bed contacts.bed
    python3 porec_to_pairs.py in.bed contacts.bed --max-order 30 --min-mapq 20
    python3 porec_to_pairs.py in.bed contacts.pa5 --format pa5

NOTE
----
YaHS accepts BED / BAM / PA5 / BIN. Both formats emitted here are valid;
BED is the default for the reason given above.
"""

import argparse
import sys
from itertools import combinations


def parse_args():
    p = argparse.ArgumentParser(
        description="Decompose Pore-C concatemers into pairwise contacts."
    )
    p.add_argument("infile", help="Input BED, sorted by read name (col 4)")
    p.add_argument("outfile", help="Output contact file")
    p.add_argument(
        "--format",
        choices=["bed", "pa5"],
        default="bed",
        help=(
            "Output format (default: bed). BED is recommended for Pore-C: "
            "it carries real alignment start/end, whereas PA5 gives only a "
            "position and makes YaHS assume a fixed --read-length."
        ),
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Discard alignment fragments below this MAPQ (default: 20)",
    )
    p.add_argument(
        "--max-order",
        type=int,
        default=30,
        help=(
            "Skip reads with more than this many mapped fragments. "
            "Very high-order concatemers explode combinatorially "
            "(N=100 gives 4950 pairs) and are often chimeric. "
            "Set to 0 to disable the cap. (default: 30)"
        ),
    )
    p.add_argument(
        "--pair-naming",
        choices=["identical", "suffix"],
        default="identical",
        help=(
            "BED output only. How to name the two mates of each contact. "
            "'identical' (default) gives both mates the SAME name, which is "
            "how a real Hi-C paired-end BAM looks and what YaHS pairs on. "
            "'suffix' appends /1 and /2. The YaHS spec calls the suffix "
            "optional, but if it compares raw names, differing names mean "
            "ZERO pairs are found and no scaffolding happens."
        ),
    )
    p.add_argument(
        "--intra-only-skip",
        action="store_true",
        help=(
            "Skip pairs where both fragments map to the SAME contig. "
            "NOT RECOMMENDED for YaHS. Intra-contig contacts look "
            "uninformative for scaffolding, but YaHS uses them to build "
            "its distance-decay model, which is how it separates real "
            "inter-contig links from background. Removing them yields "
            "'noise level 0.000' in the YaHS log and zero joins. Only "
            "use this for tools that ignore intra-contig signal."
        ),
    )
    return p.parse_args()


def emit_pa5(out, read_name, frags, skip_intra):
    """Write pairwise contacts in PA5 format.

    Columns: name  c1  pos1  c2  pos2  mapq1  mapq2

    Columns 6 and 7 are optional in the YaHS spec but worth including:
    without them YaHS cannot apply its own -q mapping-quality filter to
    PA5 input.

    NOTE: PA5 carries only a POSITION, not a span, so YaHS infers the
    alignment length from its --read-length option (default 150). That
    default is wrong for Pore-C fragments. Prefer --format bed, which
    carries real start/end coordinates and needs no such assumption.
    """
    n = 0
    for a, b in combinations(frags, 2):
        if skip_intra and a[0] == b[0]:
            continue
        # use fragment midpoint as the contact position
        pos_a = (a[1] + a[2]) // 2
        pos_b = (b[1] + b[2]) // 2
        out.write(
            f"{read_name}\t{a[0]}\t{pos_a}\t{b[0]}\t{pos_b}"
            f"\t{a[3]}\t{b[3]}\n"
        )
        n += 1
    return n


def emit_bed(out, read_name, frags, skip_intra, naming="identical"):
    """Write pairwise contacts as paired BED, two consecutive lines each.

    naming="identical": both mates share one name (mimics a real
        paired-end Hi-C BAM). This is the safe default — YaHS groups
        consecutive lines by name, and identical names always pair.
    naming="suffix": append /1 and /2. Valid per the YaHS spec, but if
        YaHS compares raw name strings the mates never match and NO
        links are built, silently producing zero joins.
    """
    n = 0
    for i, (a, b) in enumerate(combinations(frags, 2), start=1):
        if skip_intra and a[0] == b[0]:
            continue
        pair_id = f"{read_name}_p{i}"
        if naming == "suffix":
            name_a, name_b = f"{pair_id}/1", f"{pair_id}/2"
        else:
            name_a = name_b = pair_id
        out.write(f"{a[0]}\t{a[1]}\t{a[2]}\t{name_a}\t{a[3]}\t{a[4]}\n")
        out.write(f"{b[0]}\t{b[1]}\t{b[2]}\t{name_b}\t{b[3]}\t{b[4]}\n")
        n += 1
    return n


def main():
    args = parse_args()
    if args.format == "bed":
        def emit(out, name, frags, skip_intra):
            return emit_bed(out, name, frags, skip_intra, args.pair_naming)
    else:
        emit = emit_pa5

    stats = {
        "lines": 0,
        "low_mapq": 0,
        "reads": 0,
        "reads_single": 0,
        "reads_skipped_order": 0,
        "pairs": 0,
    }

    current_name = None
    frags = []

    def flush():
        """Process the accumulated fragments for one read."""
        if current_name is None:
            return
        stats["reads"] += 1
        if len(frags) < 2:
            stats["reads_single"] += 1
            return
        if args.max_order and len(frags) > args.max_order:
            stats["reads_skipped_order"] += 1
            return
        stats["pairs"] += emit(out, current_name, frags, args.intra_only_skip)

    with open(args.infile) as fh, open(args.outfile, "w") as out:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            stats["lines"] += 1

            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            name = f[3]
            score = f[4] if len(f) > 4 else "60"
            strand = f[5] if len(f) > 5 else "+"

            # MAPQ filter (bamToBed puts MAPQ in the score column)
            try:
                if int(float(score)) < args.min_mapq:
                    stats["low_mapq"] += 1
                    continue
            except ValueError:
                pass

            if name != current_name:
                flush()
                current_name = name
                frags = []

            frags.append((chrom, start, end, score, strand))

        flush()

    sys.stderr.write(
        "porec_to_pairs.py summary\n"
        f"  alignment lines read     : {stats['lines']}\n"
        f"  dropped (low MAPQ)       : {stats['low_mapq']}\n"
        f"  reads seen               : {stats['reads']}\n"
        f"  reads with <2 fragments  : {stats['reads_single']}\n"
        f"  reads over --max-order   : {stats['reads_skipped_order']}\n"
        f"  pairwise contacts written: {stats['pairs']}\n"
        f"  output format            : {args.format}\n"
        + (f"  pair naming              : {args.pair_naming}\n"
           if args.format == "bed" else "")
    )

    if stats["pairs"] == 0:
        sys.stderr.write(
            "\nWARNING: zero contacts written. Common causes:\n"
            "  - BED was not sorted by read name (column 4)\n"
            "  - supplementary alignments were filtered out of the BAM\n"
            "    (keep flag 0x800; filter only 0x4 and 0x100)\n"
            "  - --min-mapq set too high\n"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
