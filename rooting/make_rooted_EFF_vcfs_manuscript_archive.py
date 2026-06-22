"""
Make rooted VCFs from legacy CSV 'root' estimates and keep SnpEff EFF annotations consistent.

Changes vs original make_rooted_vcfs_from_old_csv_files.py:
- Robust GT flipping: flips only the GT field per sample (0<->1), preserves other FORMAT fields.
- When REF/ALT are swapped (because ALT equals the inferred ancestor), update EFF codon pairs
  for SYNONYMOUS_CODING and NON(S)_SYNONYMOUS_CODING by reversing Codon_Change (e.g., aaA/aaG -> aaG/aaA).
  This works regardless of coding strand because SnpEff codons are in coding orientation.
- Optionally updates INFO AC/AN/AF to remain consistent with flipped genotypes.
"""

import sys
import argparse
import pandas as pd
import os.path as op
import gzip
import re
from collections import Counter
from types import SimpleNamespace


EFF_RE = re.compile(r"(?:^|;)EFF=([^;]+)")


def swap_eff_codon_pairs(info: str) -> str:
    """Reverse Codon_Change (field 3) for SYNONYMOUS_CODING and NON-SYNONYMOUS_CODING tokens.

    Works on EFF=... only; returns original info if no EFF present.
    """
    m = EFF_RE.search(info or '')
    if not m:
        return info
    eff_val = m.group(1)
    # Split top-level tokens by commas while respecting parentheses
    tokens = []
    depth = 0
    buf = []
    for ch in eff_val:
        if ch == ',' and depth == 0:
            tokens.append(''.join(buf))
            buf = []
            continue
        buf.append(ch)
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
    if buf:
        tokens.append(''.join(buf))

    def should_swap(name: str) -> bool:
        # Handle multiple spellings used by snpEff variants
        n = name.upper()
        return (
            n == 'SYNONYMOUS_CODING' or
            n == 'NON_SYNONYMOUS_CODING' or
            n == 'NONSYNONYMOUS_CODING'
        )

    new_tokens = []
    for t in tokens:
        ts = t.strip()
        if '(' in ts and ts.endswith(')'):
            name, rest = ts.split('(', 1)
            name = name.strip()
            inner = rest[:-1]  # drop trailing ')'
            if should_swap(name):
                fields = inner.split('|')
                if len(fields) >= 3 and '/' in fields[2]:
                    a, b = fields[2].split('/', 1)
                    fields[2] = f"{b}/{a}"
                ts = f"{name}({ '|'.join(fields) })"
        new_tokens.append(ts)

    new_eff_val = ','.join(new_tokens)
    start, end = m.span(1)
    return info[:start] + new_eff_val + info[end:]


def is_biallelic_snp(ref: str, alt: str) -> bool:
    return (len(ref) == 1 and len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT' and ',' not in alt)


def flip_gt_field(gt: str) -> str:
    # Translate only 0/1 digits; keep separators and '.'
    return gt.translate(str.maketrans({'0': '1', '1': '0'}))


def flip_genotypes_and_recount(cols):
    """Flip GT (0<->1) per sample, recompute AC/AN, and return (cols, AC, AN).
    cols is a list of VCF columns for one record (split by '\t').
    """
    ac = 0
    an = 0
    if len(cols) <= 8 or not cols[8] or cols[8] == '.':
        return cols, ac, an
    fmt_keys = cols[8].split(':')
    key_to_idx = {k: i for i, k in enumerate(fmt_keys)}
    gt_idx = key_to_idx.get('GT')
    if gt_idx is None:
        return cols, ac, an
    # Process samples
    for si in range(9, len(cols)):
        if cols[si] in ('.', ''):
            continue
        sp = cols[si].split(':')
        if gt_idx < len(sp):
            g = sp[gt_idx]
            sp[gt_idx] = flip_gt_field(g)
            ac += sp[gt_idx].count('1')
            an += (sp[gt_idx].count('0') + sp[gt_idx].count('1'))
        cols[si] = ':'.join(sp)
    return cols, ac, an


def update_info_ac_af(info: str, ac: int, an: int) -> str:
    if not info:
        info = ''
    fields = [] if info == '' else info.split(';')
    kv = {}
    others = []
    for fld in fields:
        if '=' in fld:
            k, v = fld.split('=', 1)
            kv[k] = v
        elif fld:
            others.append(fld)
    kv['AC'] = str(ac)
    if an > 0:
        kv['AN'] = str(an)
        kv['AF'] = f"{(ac / an):.6g}"
    new_info = ';'.join([f"{k}={v}" for k, v in kv.items()] + others)
    return new_info


def load_ancestor_tables():
    ancdir = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/zimbabwe_allele_ancestor/ancestral"
    chrstrs = ["chr2L", "chr2R", "chr3L", "chr3R"]
    tables = {}
    for chrstr in chrstrs:
        csvfn = op.join(ancdir, "D_melanogaster_ancbase_" + chrstr[3:] + ".csv")
        df = pd.read_csv(csvfn, index_col=False)
        tables[chrstr] = df
    return tables


def main(args):
    anc_tables = load_ancestor_tables()
    bases = set('ACGT')
    counts = Counter()

    opener = gzip.open if args.vcf.endswith('.gz') else open
    with opener(args.vcf, 'rt') as vcf_in, open(args.out, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith('#'):
                vcf_out.write(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = cols[:8]
            if chrom not in anc_tables:
                continue
            if not is_biallelic_snp(ref, alt):
                continue
            pos1 = int(pos)
            df = anc_tables[chrom]
            # Dmel6pos1based column used in original script
            hit = df[df['Dmel6pos1based'] == pos1]
            if hit.empty:
                counts['not_found'] += 1
                continue
            row = hit.iloc[0]
            # Columns: ... probA probC probG probT EstimatedAncest ...
            # Using the same max prob >= 0.9 rule
            probs = [row.get('probA', 0), row.get('probC', 0), row.get('probG', 0), row.get('probT', 0)]
            maxprob = max(probs)
            estanc = str(row.get('EstimatedAncest', '')).upper()
            if ref in bases and alt in bases:
                if maxprob >= args.probcutoff :
                    if estanc != ref and estanc == alt:
                        # Flip REF/ALT
                        counts['flipped'] += 1
                        old_ref, old_alt = cols[3], cols[4]
                        cols[3], cols[4] = cols[4], cols[3]
                        # Update INFO EFF codon pairs for coding effects
                        cols[7] = swap_eff_codon_pairs(cols[7])
                        # Flip GT and recalc AC/AN/AF
                        cols, ac_after, an_after = flip_genotypes_and_recount(cols)
                        cols[7] = update_info_ac_af(cols[7], ac_after, an_after)
                        # Mark ID
                        cols[2] = (cols[2] + "_rooted_alt") if cols[2] != '.' else "rooted_alt"
                        vcf_out.write('\t'.join(cols) + '\n')
                    else:
                        counts['kept'] += 1
                        cols[2] = (cols[2] + "_rooted_ref") if cols[2] != '.' else "rooted_ref"
                        vcf_out.write('\t'.join(cols) + '\n')
                else:
                    counts['uncertain'] += 1
                    # write unchanged
                    vcf_out.write(line)
            else:
                # non-SNP or ambiguous, pass through
                vcf_out.write(line)

    with open(args.summary, 'w') as sf:
        sf.write(str(vars(args)) + '\n')
        for k in ['flipped', 'kept', 'uncertain', 'not_found']:
            sf.write(f"{k}: {counts.get(k, 0)}\n")


def parse_args():
    ap = argparse.ArgumentParser(description='Root VCF using legacy CSV roots and keep EFF codons consistent when flipping')
    ap.add_argument('-v', '--vcf', required=True, help='Input VCF')
    ap.add_argument('-o', '--out', required=True, help='Output VCF')
    ap.add_argument('-s', '--summary', required=True, help='Summary file')
    return ap.parse_args()


if __name__ == '__main__':
    # args = parse_args()
    args = SimpleNamespace()
    # args.vcf = "../vcf_files/ZI_2L2R3L3RX_remade_dm6_snpeff.vcf"
    # args.out = "../vcf_files/ZI_2L2R3L3RX_remade_dm6_snpeff.rooted_2_8_2026_p0.99.vcf"
    # args.summary = "../ZIoldrootmethod/ZI_2L2R3L3R_remade.dm6.snpeff.rooting_summary_p0.99.txt"

    # args.vcf = "../vcf_files/ZI_2L2R3L3R_dm6_imputed_remade_snpeff.vcf"
    # args.out = "../vcf_files/ZI_2L2R3L3R_dm6_imputed_remade_snpeff_rootedp0.99.vcf"
    # args.summary = "../ZIoldrootmethod/ZI_2L2R3L3R_dm6_imputed_remade_snpeff_rootedp0.99_summary.txt"
    # args.probcutoff = 0.99

    args.vcf = "../vcf_files/ZI_2L2R3L3R_dm6_imputed_remade_snpeff.vcf"
    args.out = "../vcf_files/ZI_2L2R3L3R_dm6_imputed_remade_snpeff_rooted0.9.vcf"
    args.summary = "../ZIimputed/ZI_2L2R3L3R_dm6_imputed_remade_snpeff_rooted0.9_summary.txt"
    args.probcutoff = 0.9
    main(args)

