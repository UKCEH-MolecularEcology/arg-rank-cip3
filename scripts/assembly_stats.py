#!/usr/bin/env python3
"""
assembly_stats.py — Metagenome assembly statistics calculator

Usage:
    python assembly_stats.py *.fasta
    python assembly_stats.py assembly1.fa assembly2.fasta assembly3.fna
    python assembly_stats.py --dir /path/to/assemblies --ext .fasta
    python assembly_stats.py *.fa --out results.tsv

Outputs a summary table to stdout (TSV) and optionally saves to a file.
"""

import sys
import os
import argparse
import gzip
from pathlib import Path


def parse_fasta(filepath):
    """Parse a FASTA file (plain or gzipped), yield (header, sequence) tuples."""
    open_fn = gzip.open if str(filepath).endswith(".gz") else open
    mode = "rt"
    header = None
    seq_parts = []
    with open_fn(filepath, mode) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - seq.count("N")
    return (gc / total * 100) if total > 0 else 0.0


def compute_nx(lengths, x=50):
    """Compute Nx (e.g. N50) for a sorted (descending) list of contig lengths."""
    total = sum(lengths)
    target = total * (x / 100)
    cumsum = 0
    for L in lengths:
        cumsum += L
        if cumsum >= target:
            return L
    return 0


def assembly_stats(filepath, min_len=0):
    lengths = []
    total_gc_bases = 0
    total_non_n = 0
    total_n = 0

    for header, seq in parse_fasta(filepath):
        L = len(seq)
        if L < min_len:
            continue
        lengths.append(L)
        seq_up = seq.upper()
        total_gc_bases += seq_up.count("G") + seq_up.count("C")
        n_count = seq_up.count("N")
        total_n += n_count
        total_non_n += L - n_count

    if not lengths:
        return None

    lengths.sort(reverse=True)
    total_bp = sum(lengths)

    stats = {
        "file": os.path.basename(filepath),
        "num_contigs": len(lengths),
        "total_bp": total_bp,
        "longest": lengths[0],
        "shortest": lengths[-1],
        "mean_len": round(total_bp / len(lengths), 1),
        "median_len": lengths[len(lengths) // 2],
        "n50": compute_nx(lengths, 50),
        "n90": compute_nx(lengths, 90),
        "l50": next(i + 1 for i, _ in enumerate(lengths) if sum(lengths[:i+1]) >= total_bp * 0.5),
        "l90": next(i + 1 for i, _ in enumerate(lengths) if sum(lengths[:i+1]) >= total_bp * 0.9),
        "contigs_ge_1kb": sum(1 for L in lengths if L >= 1000),
        "contigs_ge_10kb": sum(1 for L in lengths if L >= 10000),
        "contigs_ge_100kb": sum(1 for L in lengths if L >= 100000),
        "total_n_bases": total_n,
        "gc_percent": round((total_gc_bases / total_non_n * 100) if total_non_n > 0 else 0, 2),
    }
    return stats


def human_bp(n):
    if n >= 1_000_000_000:
        return f"{n/1e9:.2f} Gbp"
    elif n >= 1_000_000:
        return f"{n/1e6:.2f} Mbp"
    elif n >= 1_000:
        return f"{n/1e3:.2f} Kbp"
    return f"{n} bp"


def print_pretty(all_stats):
    """Print a pretty human-readable summary for each assembly."""
    for s in all_stats:
        print(f"\n{'='*60}")
        print(f"  {s['file']}")
        print(f"{'='*60}")
        print(f"  Contigs            : {s['num_contigs']:,}")
        print(f"  Total length       : {human_bp(s['total_bp'])} ({s['total_bp']:,} bp)")
        print(f"  Longest contig     : {human_bp(s['longest'])}")
        print(f"  Shortest contig    : {human_bp(s['shortest'])}")
        print(f"  Mean length        : {human_bp(int(s['mean_len']))}")
        print(f"  Median length      : {human_bp(s['median_len'])}")
        print(f"  N50                : {human_bp(s['n50'])}")
        print(f"  N90                : {human_bp(s['n90'])}")
        print(f"  L50                : {s['l50']:,} contigs")
        print(f"  L90                : {s['l90']:,} contigs")
        print(f"  Contigs ≥ 1 kbp    : {s['contigs_ge_1kb']:,}")
        print(f"  Contigs ≥ 10 kbp   : {s['contigs_ge_10kb']:,}")
        print(f"  Contigs ≥ 100 kbp  : {s['contigs_ge_100kb']:,}")
        print(f"  GC content         : {s['gc_percent']}%")
        print(f"  N bases            : {s['total_n_bases']:,}")
    print()


def write_tsv(all_stats, outfile):
    cols = [
        "file", "num_contigs", "total_bp", "longest", "shortest",
        "mean_len", "median_len", "n50", "n90", "l50", "l90",
        "contigs_ge_1kb", "contigs_ge_10kb", "contigs_ge_100kb",
        "gc_percent", "total_n_bases",
    ]
    with open(outfile, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in all_stats:
            fh.write("\t".join(str(s[c]) for c in cols) + "\n")
    print(f"TSV saved to: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute assembly statistics for metagenome FASTA files."
    )
    parser.add_argument(
        "files", nargs="*", help="FASTA file(s) (.fa, .fasta, .fna, .gz)"
    )
    parser.add_argument(
        "--dir", help="Directory containing FASTA files"
    )
    parser.add_argument(
        "--ext", default=".fasta",
        help="Extension to glob when using --dir (default: .fasta)"
    )
    parser.add_argument(
        "--min-len", type=int, default=0,
        help="Minimum contig length to include (default: 0)"
    )
    parser.add_argument(
        "--out", help="Output TSV file path (optional)"
    )
    parser.add_argument(
        "--tsv-only", action="store_true",
        help="Skip pretty printing, only write TSV"
    )

    args = parser.parse_args()

    filepaths = list(args.files)
    if args.dir:
        d = Path(args.dir)
        filepaths += [str(p) for p in d.glob(f"*{args.ext}")]
        filepaths += [str(p) for p in d.glob(f"*{args.ext}.gz")]

    if not filepaths:
        parser.print_help()
        sys.exit(1)

    all_stats = []
    for fp in sorted(filepaths):
        print(f"Processing: {fp} ...", file=sys.stderr)
        s = assembly_stats(fp, min_len=args.min_len)
        if s is None:
            print(f"  WARNING: No sequences found in {fp}", file=sys.stderr)
        else:
            all_stats.append(s)

    if not all_stats:
        print("No valid assemblies found.", file=sys.stderr)
        sys.exit(1)

    if not args.tsv_only:
        print_pretty(all_stats)

    if args.out:
        write_tsv(all_stats, args.out)
    else:
        # Always print TSV to stdout as well (for piping/redirection)
        cols = [
            "file", "num_contigs", "total_bp", "longest", "shortest",
            "mean_len", "median_len", "n50", "n90", "l50", "l90",
            "contigs_ge_1kb", "contigs_ge_10kb", "contigs_ge_100kb",
            "gc_percent", "total_n_bases",
        ]
        print("\t".join(cols))
        for s in all_stats:
            print("\t".join(str(s[c]) for c in cols))


if __name__ == "__main__":
    main()
