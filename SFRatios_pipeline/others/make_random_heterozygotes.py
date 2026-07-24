#!/usr/bin/env python3
import sys
import os
import random
from typing import List, Tuple

USAGE = """
Usage:
  make_random_heterozygotes.py <input.vcf> <output.vcf> <summary.txt>

Description:
  - Reads a haploid-sample VCF and constructs diploid individuals by pairing the
    input samples at random (without replacement).
  - If sample ID 'ZI382' is present, it is excluded prior to pairing. After that, if
    the remaining number of samples is odd, one additional sample is excluded at random.
  - New sample names are "<ID1>_<ID2>" as requested.
  - For each variant, the new genotype for a pair is GT = a|b where:
      a = first allele (0 or 1) of the first sample in the pair
      b = first allele (0 or 1) of the second sample in the pair
    If either allele is missing or not 0/1, outputs ".|." for that pair at that site.
  - Writes a summary with per-new-sample counts: hom_ref(0|0), hom_alt(1|1), het_0_1(0|1), het_1_0(1|0),
    lists the random pairings used, and any excluded samples.

Notes:
  - Input VCF must have a GT field in FORMAT. Other FORMAT fields are ignored in the output; output uses FORMAT=GT only.
  - New sample names use '_' to join IDs.
"""


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def parse_args(argv: List[str]) -> Tuple[str, str, str]:
    if len(argv) != 4 or argv[1] in {"-h", "--help"}:
        print(USAGE.strip())
        sys.exit(0 if len(argv) != 4 else 0)
    invcf, outvcf, summary = argv[1], argv[2], argv[3]
    return invcf, outvcf, summary


def parse_header_and_samples(fp) -> Tuple[List[str], List[str]]:
    meta = []
    for line in fp:
        if line.startswith('##'):
            meta.append(line.rstrip('\n'))
            continue
        if line.startswith('#CHROM'):
            hdr = line.rstrip('\n').split('\t')
            if len(hdr) < 10:
                die("VCF header has fewer than 10 columns; no samples found?")
            samples = hdr[9:]
            return meta, samples
        die("Malformed VCF: no #CHROM header before data")
    die("Empty VCF or missing #CHROM header")


def ensure_gt_format(meta_lines: List[str]) -> List[str]:
    # Keep meta as-is, but add GT FORMAT definition if missing
    has_gt = any(l.startswith('##FORMAT=') and 'ID=GT,' in l for l in meta_lines)
    if not has_gt:
        meta_lines = list(meta_lines) + [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        ]
    return meta_lines


def first_allele_from_gt(gt: str) -> str:
    if not gt or gt == '.' or gt == './.' or gt == '.|.':
        return '.'
    # Accept patterns like 0, 1, 0/0, 1|0, etc. Return first allele digit if 0/1, else '.'
    for ch in gt:
        if ch in ('0', '1'):
            return ch
        if ch in ('/', '|'):
            continue
    return '.'


def build_pairs(samples: List[str]) -> List[Tuple[str, str]]:
    n = len(samples)
    if n % 2 != 0:
        die(f"Number of samples is odd ({n}); cannot pair without replacement")
    idx = list(range(n))
    random.shuffle(idx)
    pairs = []
    for i in range(0, n, 2):
        s1 = samples[idx[i]]
        s2 = samples[idx[i+1]]
        pairs.append((s1, s2))
    return pairs


def adjust_samples(samples: List[str]) -> Tuple[List[str], List[str]]:
    """
    Exclude 'ZI382' if present. If count remains odd, exclude one random sample.
    Returns (filtered_samples, excluded_list).
    """
    excluded = []
    work = list(samples)
    if 'ZI382' in work:
        work.remove('ZI382')
        excluded.append('ZI382')
    if len(work) % 2 == 1:
        victim = random.choice(work)
        work.remove(victim)
        excluded.append(victim)
    return work, excluded


def main(argv: List[str]) -> None:
    invcf, outvcf, sumfile = parse_args(argv)
    if not os.path.isfile(invcf):
        die(f"Input VCF not found: {invcf}")

    with open(invcf, 'r') as fin:
        meta, orig_samples = parse_header_and_samples(fin)
        if len(orig_samples) == 0:
            die("No samples in VCF")

        samples, excluded_samples = adjust_samples(orig_samples)
        if len(samples) == 0:
            die("No samples left after exclusions")

        pairs = build_pairs(samples)
        # New sample IDs using '_'
        new_samples = [f"{a}_{b}" for a, b in pairs]

        meta = ensure_gt_format(meta)

        # Precompute sample index mapping for speed (relative to original sample order)
        name_to_index = {s: i for i, s in enumerate(orig_samples)}

        # Prepare summary counters
        hom0 = {ns: 0 for ns in new_samples}
        hom1 = {ns: 0 for ns in new_samples}
        het01 = {ns: 0 for ns in new_samples}
        het10 = {ns: 0 for ns in new_samples}

        # Write output VCF
        with open(outvcf, 'w') as fout:
            # Meta headers
            for m in meta:
                fout.write(m + '\n')
            # Optional provenance line
            fout.write('##source=make_random_heterozygotes.py\n')
            # Header line with new samples
            hdr_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
            fout.write('\t'.join(hdr_cols + new_samples) + '\n')

            # Now process remaining lines (fin is positioned after #CHROM)
            for line in fin:
                if not line or line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 10:
                    continue  # malformed
                chrom,pos,vid,ref,alt,qual,fil,info,fmt = cols[:9]
                samp_fields = cols[9:]
                fmt_keys = fmt.split(':') if fmt else []
                try:
                    gt_idx = fmt_keys.index('GT')
                except ValueError:
                    gt_idx = -1

                # Extract first allele for each original sample (0/1 or '.')
                first_alleles = []
                for sf in samp_fields:
                    if gt_idx == -1:
                        first_alleles.append('.')
                        continue
                    parts = sf.split(':')
                    gt = parts[gt_idx] if gt_idx < len(parts) else '.'
                    a = first_allele_from_gt(gt)
                    if a not in ('0', '1'):
                        a = '.'
                    first_alleles.append(a)
                if len(first_alleles) != len(orig_samples):
                    continue  # malformed row

                # Build new diploid genotypes
                new_gts = []
                for (a,b) in pairs:
                    i = name_to_index[a]
                    j = name_to_index[b]
                    fa = first_alleles[i]
                    fb = first_alleles[j]
                    if fa in ('0','1') and fb in ('0','1'):
                        gt = f"{fa}|{fb}"
                        new_gts.append(gt)
                        ns = f"{a}_{b}"
                        if fa == '0' and fb == '0':
                            hom0[ns] += 1
                        elif fa == '1' and fb == '1':
                            hom1[ns] += 1
                        elif fa == '0' and fb == '1':
                            het01[ns] += 1
                        elif fa == '1' and fb == '0':
                            het10[ns] += 1
                    else:
                        new_gts.append('.|.')
                # Write line with FORMAT=GT only
                out_cols = [chrom,pos,vid,ref,alt,qual,fil,info,'GT'] + new_gts
                fout.write('\t'.join(out_cols) + '\n')

        # Write summary file
        with open(sumfile, 'w') as sf:
            sf.write("# Excluded samples\n")
            if excluded_samples:
                for ex in excluded_samples:
                    sf.write(f"{ex}\n")
            else:
                sf.write("(none)\n")
            sf.write("\n# Random pairs used (NewSample\tOrig1\tOrig2)\n")
            for a,b in pairs:
                ns = f"{a}_{b}"
                sf.write(f"{ns}\t{a}\t{b}\n")
            sf.write("\n# Per-new-sample genotype counts (hom0=0|0, hom1=1|1, het01=0|1, het10=1|0)\n")
            sf.write("Sample\thom0\thom1\thet01\thet10\n")
            for ns in new_samples:
                sf.write(f"{ns}\t{hom0[ns]}\t{hom1[ns]}\t{het01[ns]}\t{het10[ns]}\n")


if __name__ == '__main__':
    main(sys.argv)
