#!/usr/bin/env python3
"""
codex script

Estimate the fraction of D. melanogaster autosomes that align (in a MAF) to valid bases in:
  1) both Anc0 and Anc2
  2) D_simulans

Inputs:
  - Reference FASTA: Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa (for autosome lengths)
  - MAF (gz): Dmelanogaster.maf.gz

Assumptions:
  - MAF is anchored on D. melanogaster.
  - Species names in MAF 's' lines contain identifiers like 'melanogaster', 'Anc0', 'Anc2', 'simulans'.
  - Chromosome names use some combination of 2L/2R/3L/3R with optional 'chr' prefix.

Output: prints overall fractions and per-arm breakdown; optionally writes to a file.
"""

import sys
import os
import argparse
import gzip
from collections import defaultdict

AUTOSOMES = { '2L', '2R', '3L', '3R' }
# Map D. melanogaster RefSeq chromosome sizes (dm6) to arms
MEL_ARM_BY_SIZE = {
    23513712: '2L',
    25286936: '2R',
    28110227: '3L',
    32079331: '3R',
}
VALID_BASE = set('ACGTacgt')

def norm_chrom(name: str) -> str:
    # Accept '2L' or 'chr2L'; return plain arm (e.g., '2L') or '' if not autosome
    x = name
    if x.lower().startswith('chr'):
        x = x[3:]
    return x if x in AUTOSOMES else ''

def fasta_autosome_lengths(fasta_path: str):
    lens = {a: 0 for a in AUTOSOMES}
    cur = None
    with open(fasta_path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                # Parse contig name up to first whitespace
                name = line[1:].strip().split()[0]
                arm = norm_chrom(name)
                cur = arm if arm else None
            else:
                if cur:
                    lens[cur] += sum(c in 'ACGTNacgtn' for c in line.strip())
    total = sum(lens.values())
    return lens, total

def parse_s_line(line: str):
    # MAF s-line: s src start size strand srcSize text
    parts = line.strip().split()
    # guard
    if len(parts) < 7:
        return None
    _, src, start, size, strand, srcSize, text = parts[:7]
    # Identify species and chrom
    sp = src.split('.')[0]
    chrom = src.split('.')[-1] if '.' in src else src
    return {
        'src': src,
        'species': sp,
        'chrom': chrom,
        'start': int(start),
        'size': int(size),
        'strand': strand,
        'srcSize': int(srcSize),
        'text': text,
    }

def label_of_species(species: str):
    s = species.lower()
    if 'melanogaster' in s or s in ('dm6','dmel','d_melanogaster','dmelanogaster'):
        return 'mel'
    if s.startswith('anc0'):
        return 'anc0'
    if s.startswith('anc2'):
        return 'anc2'
    if 'simulans' in s:
        return 'dsim'
    return ''

def estimate(fasta_path: str, maf_gz_path: str):
    lens, total_len = fasta_autosome_lengths(fasta_path)

    # counters overall
    count_anc_both = 0
    count_dsim = 0
    # per arm counters
    per_arm_anc_both = {a: 0 for a in AUTOSOMES}
    per_arm_dsim = {a: 0 for a in AUTOSOMES}

    with gzip.open(maf_gz_path, 'rt', encoding='utf-8', errors='replace') as fh:
        block = []
        for raw in fh:
            line = raw.rstrip('\n')
            if not line.strip():
                if block:
                    process_block(block, per_arm_anc_both, per_arm_dsim)
                    block = []
                continue
            first = line.split(None, 1)[0]
            if first == 'a':
                if block:
                    process_block(block, per_arm_anc_both, per_arm_dsim)
                    block = []
                continue
            if first == 's':
                block.append(line)
        if block:
            process_block(block, per_arm_anc_both, per_arm_dsim)

    # sum totals
    count_anc_both = sum(per_arm_anc_both.values())
    count_dsim = sum(per_arm_dsim.values())

    return {
        'lens': lens,
        'total_len': total_len,
        'anc_both': count_anc_both,
        'dsim': count_dsim,
        'per_arm_anc_both': per_arm_anc_both,
        'per_arm_dsim': per_arm_dsim,
    }

def process_block(s_lines, per_arm_anc_both, per_arm_dsim):
    # Collect sequences by label
    recs = defaultdict(list)
    mel_arm = ''
    for line in s_lines:
        rec = parse_s_line(line)
        if not rec:
            continue
        sp_lab = label_of_species(rec['species'])
        # attempt to detect mel by species label first
        if sp_lab == 'mel':
            arm = norm_chrom(rec['chrom'])
            if not arm:
                arm = MEL_ARM_BY_SIZE.get(rec['srcSize'], '')
            if arm in AUTOSOMES:
                mel_arm = arm
        recs[sp_lab].append(rec)

    if not mel_arm:
        # fallback: try any 's' whose chrom looks autosomal and species unknown
        for rec in [parse_s_line(x) for x in s_lines]:
            if not rec:
                continue
            arm = norm_chrom(rec['chrom'])
            if not arm:
                arm = MEL_ARM_BY_SIZE.get(rec['srcSize'], '')
            sp_lab = label_of_species(rec['species'])
            if arm in AUTOSOMES and sp_lab == '':
                mel_arm = arm
                recs['mel'].append(rec)
                break

    if not mel_arm:
        return

    # use the first occurrence per label (MAF typically has one per block)
    mel = recs['mel'][0] if recs.get('mel') else None
    a0  = recs['anc0'][0] if recs.get('anc0') else None
    a2  = recs['anc2'][0] if recs.get('anc2') else None
    ds  = recs['dsim'][0] if recs.get('dsim') else None
    if not mel:
        return

    t_m = mel['text']
    aln_len = len(t_m)
    # Check all texts have same length if present
    for t in (a0['text'] if a0 else '', a2['text'] if a2 else '', ds['text'] if ds else ''):
        if t and len(t) != aln_len:
            return

    # Iterate columns
    for i in range(aln_len):
        cm = t_m[i]
        if cm == '-':
            continue  # mel gap
        # Anc both
        if a0 and a2:
            c0 = a0['text'][i]
            c2 = a2['text'][i]
            if (c0 in VALID_BASE) and (c2 in VALID_BASE):
                per_arm_anc_both[mel_arm] += 1
        # D. simulans
        if ds:
            cds = ds['text'][i]
            if cds in VALID_BASE:
                per_arm_dsim[mel_arm] += 1

def main():
    ap = argparse.ArgumentParser(description='Estimate autosome fractions with valid bases in Anc0&Anc2 and D_simulans from MAF')
    ap.add_argument('-f', '--fasta', required=True, help='D. melanogaster FASTA path (BDGP6.54)')
    ap.add_argument('-m', '--maf', required=True, help='MAF gz path')
    ap.add_argument('-o', '--out', help='Output file to write summary')
    args = ap.parse_args()

    res = estimate(args.fasta, args.maf)

    total = res['total_len']
    frac_anc = res['anc_both'] / total if total else 0.0
    frac_dsim = res['dsim'] / total if total else 0.0

    lines = []
    lines.append('Autosome lengths (bp): ' + ', '.join(f"{k}={v}" for k,v in sorted(res['lens'].items())))
    lines.append(f"Total autosome length: {total}")
    lines.append(f"Anc0&Anc2 valid base count: {res['anc_both']}  fraction: {frac_anc:.6f}")
    lines.append(f"D_simulans valid base count: {res['dsim']}  fraction: {frac_dsim:.6f}")
    lines.append('Per-arm counts:')
    for arm in sorted(AUTOSOMES):
        la = res['per_arm_anc_both'][arm]
        ld = res['per_arm_dsim'][arm]
        L = res['lens'][arm]
        fa = la / L if L else 0.0
        fd = ld / L if L else 0.0
        lines.append(f"  {arm}: Anc0&Anc2={la} ({fa:.6f})  D_simulans={ld} ({fd:.6f})")

    out_text = '\n'.join(lines) + '\n'
    if args.out:
        with open(args.out, 'w') as f:
            f.write(out_text)
    else:
        sys.stdout.write(out_text)

if __name__ == '__main__':
    main()
