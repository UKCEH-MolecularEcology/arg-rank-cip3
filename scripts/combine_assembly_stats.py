#!/usr/bin/env python3
"""
combine_stats.py — Merge multiple assembly_stats.py TSV outputs into one table.

Usage:
    python combine_stats.py *.tsv
    python combine_stats.py --dir /path/to/results --out combined.tsv
    python combine_stats.py results/*.tsv --out combined.tsv --sort n50
"""

import sys
import os
import argparse
import csv
from pathlib import Path


def read_tsv(filepath):
    rows = []
    with open(filepath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Combine multiple assembly_stats TSV files into one table."
    )
    parser.add_argument("files", nargs="*", help="TSV file(s) to combine")
    parser.add_argument("--dir", help="Directory to search for TSV files")
    parser.add_argument("--out", default="combined_stats.tsv", help="Output file (default: combined_stats.tsv)")
    parser.add_argument(
        "--sort", default=None,
        help="Column to sort by, e.g. n50, total_bp, num_contigs (descending)"
    )
    parser.add_argument(
        "--add-source", action="store_true",
        help="Add a 'source_file' column showing which TSV each row came from"
    )
    args = parser.parse_args()

    filepaths = list(args.files)
    if args.dir:
        filepaths += [str(p) for p in Path(args.dir).glob("*.tsv")]

    if not filepaths:
        parser.print_help()
        sys.exit(1)

    all_rows = []
    for fp in sorted(filepaths):
        if not os.path.isfile(fp):
            print(f"WARNING: {fp} not found, skipping.", file=sys.stderr)
            continue
        rows = read_tsv(fp)
        if not rows:
            print(f"WARNING: {fp} is empty, skipping.", file=sys.stderr)
            continue
        if args.add_source:
            for row in rows:
                row["source_file"] = os.path.basename(fp)
        all_rows.extend(rows)
        print(f"  Loaded {len(rows):>4} row(s) from {fp}", file=sys.stderr)

    if not all_rows:
        print("No data found.", file=sys.stderr)
        sys.exit(1)

    # Sort if requested
    if args.sort:
        col = args.sort.lower()
        if col not in all_rows[0]:
            print(f"WARNING: sort column '{col}' not found. Available: {', '.join(all_rows[0].keys())}", file=sys.stderr)
        else:
            try:
                all_rows.sort(key=lambda r: float(r[col]), reverse=True)
            except ValueError:
                all_rows.sort(key=lambda r: r[col], reverse=True)

    cols = list(all_rows[0].keys())
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nCombined {len(all_rows)} assemblies into: {args.out}", file=sys.stderr)

    # Print a quick summary to stdout
    numeric_cols = ["num_contigs", "total_bp", "n50", "longest", "gc_percent"]
    print(f"\n{'Sample':<40} {'Contigs':>10} {'Total bp':>14} {'N50':>12} {'Longest':>12} {'GC%':>6}")
    print("-" * 96)
    for r in all_rows:
        name = r.get("file", "?")[:39]
        print(
            f"{name:<40}"
            f"{int(float(r.get('num_contigs',0))):>10,}"
            f"{int(float(r.get('total_bp',0))):>14,}"
            f"{int(float(r.get('n50',0))):>12,}"
            f"{int(float(r.get('longest',0))):>12,}"
            f"{float(r.get('gc_percent',0)):>6.2f}"
        )


if __name__ == "__main__":
    main()
