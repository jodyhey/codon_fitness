#!/usr/bin/env python3
"""
Sum all intron SFSs in the specified SFS file and write a single summed SFS.

Input (hardcoded):
  /mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/
    ZIsuccess/Rootp0.9/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026_atrandom_unfolded_SFSs.txt

Output (same folder):
  /mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/
    ZIsuccess/Rootp0.9/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026_atrandom_unfolded_SFSs_summed_intron_SFS.txt

The input is expected to have repeating blocks like:
  Synonymous <CODONPAIR>
  <160 numbers>
  Intron for <CODONPAIR>
  <160 numbers>

We sum element-wise across every "Intron for ..." line that follows such a header.
"""
from __future__ import annotations
import os
import sys

IN_PATH = (
    "/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/"
    "ZIsuccess/Rootp0.9/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026_atrandom_unfolded_SFSs.txt"
)

OUT_PATH = (
    "/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/"
    "ZIsuccess/Rootp0.9/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026_atrandom_unfolded_SFSs_summed_intron_SFS.txt"
)


def parse_numbers_line(s: str):
    parts = s.strip().replace('\t', ' ').split()
    vals = []
    for p in parts:
        try:
            v = float(p)
        except ValueError:
            continue
        vals.append(v)
    return vals


def main():
    if not os.path.exists(IN_PATH):
        print(f"Input not found: {IN_PATH}")
        sys.exit(1)

    with open(IN_PATH, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f]

    n = len(lines)
    i = 0
    summed = None
    count_intron = 0
    bad_blocks = 0

    while i < n:
        line = lines[i].strip()
        if line.startswith('Intron for '):
            # Next non-empty line should contain the SFS values
            j = i + 1
            while j < n and not lines[j].strip():
                j += 1
            if j >= n:
                break
            vals = parse_numbers_line(lines[j])
            if not vals:
                bad_blocks += 1
                i = j + 1
                continue
            if summed is None:
                summed = [0.0] * len(vals)
            # Pad if inconsistent length (unlikely); keep the min common length
            L = min(len(summed), len(vals))
            for k in range(L):
                summed[k] += vals[k]
            count_intron += 1
            i = j + 1
        else:
            i += 1

    if summed is None:
        print("No intron SFS blocks found.")
        sys.exit(1)

    # Write output
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w') as out:
        out.write("Summed_Intron_SFS\n")
        out.write(' '.join(f"{v:.6g}" for v in summed) + "\n")

    print(f"Summed {count_intron} intron SFS blocks; wrote: {OUT_PATH}")
    if bad_blocks:
        print(f"Skipped {bad_blocks} malformed intron SFS blocks.")


if __name__ == '__main__':
    main()

