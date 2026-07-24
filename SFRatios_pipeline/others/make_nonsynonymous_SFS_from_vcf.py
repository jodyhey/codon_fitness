#!/usr/bin/env python3
"""
Build a Site Frequency Spectrum (SFS) for NONSYNONYMOUS coding SNPs from a VCF.

Requirements implemented:
- Use only biallelic SNPs (single-nucleotide REF and single ALT; no multi-allelic).
- Count ALT alleles using only the first allele of each genotype (haploid lines: use the first allele before '/' or '|').
- Ignore fixed sites in the final SFS (i.e., do not include bin = target_n).
- Downsample without replacement to a target haploid sample size (default 160) when a site has more than that many called haploid genotypes.
- Skip sites with fewer than the target haploid sample size.
- Filter to nonsynonymous coding using either old SnpEff EFF tag (NONSYNONYMOUS_CODING) or newer ANN terms (e.g., missense_variant, stop_gained, stop_lost, start_lost).

Output: one line with counts for bins 0..(target_n-1) separated by spaces.
Optional stats summary can be printed or written to a file.

Example:
  python make_nonsynonymous_SFS_from_vcf.py \
    -i ../vcf_files/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026.vcf \
    -o ZI_nonsynonymous_SFS_nc160.txt -n 160 -s 1
"""

import sys
import os
import gzip
import argparse
import numpy as np
from typing import Optional, Tuple, Dict


NONSYN_ANN_TERMS = {
    'missense_variant',
    'stop_gained',
    'stop_lost',
    'start_lost',
    'rare_amino_acid_variant',
}


def is_biallelic_snp(ref: str, alt: str) -> bool:
    if ',' in alt:
        return False
    return len(ref) == 1 and len(alt) == 1 and ref != '.' and alt != '.'


def is_nonsynonymous(info: str) -> bool:
    """Return True if INFO suggests a nonsynonymous coding SNP.

    Accept if either:
    - 'NONSYNONYMOUS_CODING' appears (older SnpEff EFF tag), or
    - ANN contains a typical nonsynonymous term.
    """
    L = info.lower()
    if 'nonsynonymous_coding' in L:
        return True
    # Quick ANN check without full CSV parsing
    for term in NONSYN_ANN_TERMS:
        if term in L:
            return True
    return False


def parse_gt_first_allele(sample_field: str, gt_index: int) -> Optional[str]:
    """Extract the first allele of the GT for a sample field ('0' or '1').

    Returns '0' or '1' or None if missing/invalid.
    """
    if sample_field == '.' or sample_field == './.':
        return None
    parts = sample_field.split(':')
    if gt_index >= len(parts):
        return None
    gt = parts[gt_index]
    if gt == '.' or gt == './.' or gt == '.|.':
        return None
    # Split on '/' or '|'
    sep = '/' if '/' in gt else ('|' if '|' in gt else None)
    if sep is None:
        # Haploid like '0' or '1'
        allele = gt.strip()
    else:
        allele = gt.split(sep)[0].strip()
    return allele if allele in ('0', '1') else None


def downsample_alt_count(alt_count: int, n_called: int, target_n: int, rng: np.random.Generator) -> int:
    """Downsample ALT count to target_n using hypergeometric (without replacement)."""
    if n_called == target_n:
        return alt_count
    # numpy hypergeometric: ngood, nbad, nsample
    ngood = alt_count
    nbad = n_called - alt_count
    draws = target_n
    if ngood < 0 or nbad < 0 or draws < 0 or (ngood + nbad) < draws:
        # Fallback guard; shouldn't happen if inputs are correct
        p = 0.0 if n_called == 0 else (alt_count / n_called)
        return int(round(p * target_n))
    return int(rng.hypergeometric(ngood, nbad, draws))


def build_sfs(vcf_path: str, target_n: int = 160, seed: int = 1) -> Tuple[list[int], Dict[str, int]]:
    """Stream the VCF and accumulate an SFS for nonsynonymous biallelic SNPs.

    Returns a list of length target_n with counts for bins 0..(target_n-1).
    Fixed-alt (count == target_n) is excluded by construction.
    Sites with fewer than target_n called haploid genotypes are skipped.
    """
    sfs = [0] * target_n
    rng = np.random.default_rng(seed)

    stats = {
        'total_records': 0,
        'skipped_not_biallelic': 0,
        'skipped_not_nonsyn': 0,
        'skipped_no_GT': 0,
        'skipped_lt_target': 0,
        'skipped_fixed_ALT': 0,
        'processed_sites': 0,
    }

    # Open plain or gz
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt', encoding='utf-8', errors='ignore') as f:
        gt_index = None
        for line in f:
            if not line or line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                # FORMAT is column 9 (0-based index 8), samples start at 10
                # We'll determine GT index per-variant since FORMAT can change
                continue

            fields = line.rstrip('\n').split('\t')
            stats['total_records'] += 1
            if len(fields) < 10:
                continue
            chrom, pos, var_id, ref, alt, qual, flt, info, fmt = fields[:9]
            samples = fields[9:]

            # biallelic SNP only
            if not is_biallelic_snp(ref, alt):
                stats['skipped_not_biallelic'] += 1
                continue

            # nonsynonymous coding filter
            if not is_nonsynonymous(info):
                stats['skipped_not_nonsyn'] += 1
                continue

            # Determine GT index in FORMAT
            fmt_keys = fmt.split(':')
            try:
                gt_index = fmt_keys.index('GT')
            except ValueError:
                # No GT in this record
                stats['skipped_no_GT'] += 1
                continue

            # Count ALT using first allele per genotype
            alt_count = 0
            n_called = 0
            for sf in samples:
                allele = parse_gt_first_allele(sf, gt_index)
                if allele is None:
                    continue
                n_called += 1
                if allele == '1':
                    alt_count += 1

            # Skip if not enough called haploid genotypes
            if n_called < target_n:
                stats['skipped_lt_target'] += 1
                continue

            # Downsample if necessary
            a160 = downsample_alt_count(alt_count, n_called, target_n, rng)

            # Ignore fixed ALT
            if a160 >= target_n:
                stats['skipped_fixed_ALT'] += 1
                continue

            # Increment SFS bin (0..target_n-1)
            sfs[a160] += 1
            stats['processed_sites'] += 1

    return sfs, stats


def main():
    ap = argparse.ArgumentParser(description='Build SFS for nonsynonymous coding SNPs from a VCF (ALT counts, haploid-first allele).')
    ap.add_argument('-i', '--input', dest='vcf_path', required=True, help='Input VCF path')
    ap.add_argument('-o', '--output', dest='out_path', required=False, help='Output file path (default: stdout)')
    ap.add_argument('-n', '--target-n', dest='target_n', type=int, default=160, help='Target haploid sample size (default: 160)')
    ap.add_argument('-s', '--seed', dest='seed', type=int, default=1, help='Random seed for downsampling (default: 1)')
    ap.add_argument('--stats', action='store_true', help='Print a summary of processing stats to stderr')
    ap.add_argument('--stats-file', dest='stats_file', help='Write processing stats to this file (JSON-like text)')
    args = ap.parse_args()

    if not os.path.exists(args.vcf_path):
        print(f"Error: VCF not found: {args.vcf_path}", file=sys.stderr)
        sys.exit(1)

    sfs, stats = build_sfs(args.vcf_path, target_n=args.target_n, seed=args.seed)
    out_line = ' '.join(str(x) for x in sfs)

    if args.out_path:
        with open(args.out_path, 'w') as fo:
            fo.write(out_line + '\n')
    else:
        print(out_line)

    if args.stats or args.stats_file:
        # Minimal text block (readable JSON-like)
        lines = [
            f"total_records: {stats['total_records']}",
            f"processed_sites: {stats['processed_sites']}",
            f"skipped_not_biallelic: {stats['skipped_not_biallelic']}",
            f"skipped_not_nonsyn: {stats['skipped_not_nonsyn']}",
            f"skipped_no_GT: {stats['skipped_no_GT']}",
            f"skipped_lt_target: {stats['skipped_lt_target']}",
            f"skipped_fixed_ALT: {stats['skipped_fixed_ALT']}",
        ]
        text = "\n".join(lines) + "\n"
        if args.stats:
            print(text, file=sys.stderr, end='')
        if args.stats_file:
            with open(args.stats_file, 'w') as sf:
                sf.write(text)


if __name__ == '__main__':
    main()
