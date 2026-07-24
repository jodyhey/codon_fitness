#!/usr/bin/env python3
"""Convert SINGER text ARG samples to tskit with valid initial mutation rows."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tskit


def load_2d(path: Path) -> np.ndarray:
    values = np.loadtxt(path)
    return np.atleast_2d(values)


def read_arg(node_file: Path, branch_file: Path, mutation_file: Path) -> tskit.TreeSequence:
    node_times = np.atleast_1d(np.loadtxt(node_file))
    edges = load_2d(branch_file)
    edges = edges[edges[:, 2] >= 0]
    if len(edges) == 0:
        raise ValueError(f"No valid edges in {branch_file}")

    tables = tskit.TableCollection(sequence_length=float(np.max(edges[:, 1])))

    previous_internal_time = -1.0
    for raw_time in node_times:
        time = float(raw_time)
        if time == 0:
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        else:
            time = max(previous_internal_time + 1e-4, time)
            tables.nodes.add_row(time=time)
            previous_internal_time = time

    for left, right, parent, child, *_ in edges:
        if left >= right:
            raise ValueError(
                f"Invalid edge interval [{left}, {right}) in {branch_file}"
            )
        tables.edges.add_row(
            left=float(left),
            right=float(right),
            parent=int(parent),
            child=int(child),
        )

    # Sort the tree topology before adding mutations.
    tables.sort()

    mutations = load_2d(mutation_file)
    current_position: float | None = None
    site_id = tskit.NULL
    for row in mutations:
        position = float(row[0])
        if current_position is None or position != current_position:
            site_id = tables.sites.add_row(position=position, ancestral_state="0")
            current_position = position
        tables.mutations.add_row(
            site=site_id,
            node=int(row[1]),
            derived_state=str(int(row[3])),
            parent=tskit.NULL,
            time=tskit.UNKNOWN_TIME,
        )

    # tskit 1.0 validates recurrent mutations strictly. SINGER can emit
    # multiple mutations at a site, including back-mutations on descendant
    # branches. Compute the mutation.parent links before materializing the
    # TreeSequence, then sort so parent mutation IDs precede child IDs.
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.sort()

    return tables.tree_sequence()


def convert(
    input_prefix: str,
    output_prefix: str,
    start: int,
    end: int,
    step: int,
    fast: bool,
) -> None:
    infix = "_fast" if fast else ""
    for index in range(start, end, step):
        input_base = f"{input_prefix}{infix}"
        ts = read_arg(
            Path(f"{input_base}_nodes_{index}.txt"),
            Path(f"{input_base}_branches_{index}.txt"),
            Path(f"{input_base}_muts_{index}.txt"),
        )
        output = Path(f"{output_prefix}_{index}.trees")
        ts.dump(output)
        print(f"Wrote {output}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert SINGER text ARG output to tskit."
    )
    parser.add_argument("-input", required=True, help="Prefix of ARG files")
    parser.add_argument("-output", required=True, help="Prefix of output trees")
    parser.add_argument("-start", required=True, type=int)
    parser.add_argument("-end", required=True, type=int)
    parser.add_argument("-step", type=int, default=1)
    parser.add_argument("-fast", action="store_true")
    args = parser.parse_args()

    convert(
        args.input,
        args.output,
        args.start,
        args.end,
        args.step,
        args.fast,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
