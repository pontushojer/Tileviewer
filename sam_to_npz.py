#!/usr/bin/env python3
"""Convert a paired-end SAM file into the sparse NPZ used by the genome explorer.

Usage:
    python sam_to_npz.py --sam aln.sam --fai ref.fa.fai --out matrix.npz

Matrix cell (x, y) holds the 1-based absolute genome coordinate (1..total
reference length) of the pair midpoint, or an explicitly stored 0 when the
cluster is unmapped or fails the filters. Empty grid cells are implicit zeros.
"""

import argparse
import sys
from array import array

import numpy as np
import scipy.sparse as sp


def load_fai_offsets(fai_path):
    """Return {chrom: cumulative offset} and the total reference length."""
    offsets = {}
    current = 0
    with open(fai_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                raise ValueError(f"Invalid FAI record on line {line_number}")
            if cols[0] in offsets:
                raise ValueError(f"Duplicate reference {cols[0]!r} in FAI")
            offsets[cols[0]] = current
            current += int(cols[1])
    return offsets, current


def pair_to_position(r1, r2, offsets, max_insert_size):
    """Return the 1-based absolute pair midpoint, or 0 if the pair fails."""
    if not (r1["valid"] and r2["valid"]) or r1["chr"] != r2["chr"]:
        return 0
    if abs(r1["pos"] - r2["pos"]) > max_insert_size:
        return 0
    # SAM POS and the cumulative offset make this 1-based already.
    return offsets[r1["chr"]] + (r1["pos"] + r2["pos"]) // 2


def parse_positions(sam_path, offsets, max_insert_size, min_mapq, stats):
    """Yield (x, y, absolute_position) once per read pair; 0 means unmapped."""
    pending = {}
    with open(sam_path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qname, flag = cols[0], int(cols[1])
            # Secondary/supplementary records must not overwrite the primary one.
            if flag & 0x900:
                continue
            if flag & 0x40:
                read_num = 1
            elif flag & 0x80:
                read_num = 2
            else:
                continue

            rname, pos, mapq = cols[2], int(cols[3]), int(cols[4])
            valid = rname != "*" and not flag & 0x4 and mapq >= min_mapq
            if valid and rname not in offsets:
                stats["unknown_chrom"] += 1
                valid = False

            reads = pending.setdefault(qname, {})
            reads[read_num] = {"chr": rname, "pos": pos, "valid": valid}
            if len(reads) == 2:
                del pending[qname]
                x, y = parse_xy(qname)
                yield x, y, pair_to_position(reads[1], reads[2], offsets, max_insert_size)

    # A cluster with only one record in the SAM cannot form a valid pair.
    stats["incomplete"] = len(pending)
    for qname in pending:
        x, y = parse_xy(qname)
        yield x, y, 0


def parse_xy(qname):
    *_, x, y = qname.split(":")
    return int(x), int(y)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sam", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--out", required=True, help="Output .npz path")
    parser.add_argument("--max-insert-size", type=int, default=2000)
    parser.add_argument("--min-mapq", type=int, default=40)
    args = parser.parse_args()

    offsets, total = load_fai_offsets(args.fai)
    stats = {"unknown_chrom": 0, "incomplete": 0}
    xs, ys, values = array("q"), array("q"), array("q")
    for x, y, position in parse_positions(
        args.sam, offsets, args.max_insert_size, args.min_mapq, stats
    ):
        xs.append(x)
        ys.append(y)
        values.append(position)

    if not xs:
        sys.exit("No read pairs found.")
    x = np.frombuffer(xs, dtype=np.int64)
    y = np.frombuffer(ys, dtype=np.int64)
    data = np.frombuffer(values, dtype=np.int64)
    if x.min() < 0 or y.min() < 0:
        sys.exit("Negative grid coordinates found.")
    shape = (int(x.max()) + 1, int(y.max()) + 1)

    # Duplicate (x, y) keys would be summed by scipy and corrupt coordinates
    # (e.g. clusters from different tiles). Keep one entry per cell, preferring mapped.
    order = np.lexsort((-data, y, x))
    x, y, data = x[order], y[order], data[order]
    keep = np.ones(len(x), dtype=bool)
    keep[1:] = (x[1:] != x[:-1]) | (y[1:] != y[:-1])
    duplicates = int(len(x) - keep.sum())
    x, y, data = x[keep], y[keep], data[keep]

    if data.max() > total:
        sys.exit("Computed coordinate exceeds the total reference length.")

    # CSR conversion preserves explicit zeros, which mark unmapped clusters.
    matrix = sp.coo_matrix((data, (x, y)), shape=shape, dtype=np.int64).tocsr()
    sp.save_npz(args.out, matrix)

    print(f"Grid shape: {shape}")
    print(f"Clusters: {len(data):,} (mapped {np.count_nonzero(data):,}, "
          f"unmapped {int((data == 0).sum()):,})")
    print(f"Incomplete pairs treated as unmapped: {stats['incomplete']:,}")
    print(f"Reads on references missing from FAI: {stats['unknown_chrom']:,}")
    print(f"Duplicate grid coordinates dropped: {duplicates:,}")


if __name__ == "__main__":
    main()