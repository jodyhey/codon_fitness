"""
Fixed version of make_codon_pair_SFS_from_SNP_paired_allele_counts.py

Changes:
- Robust header detection in input parsing (no accidental drop of first data row).
- Correct ordering when writing SFSs: Synonymous first, then Intron, matching labels.

Input format (tab-separated):
    codonpair  Synnc  SynCount  Intronnc  IntronCount

Output format:
    Synonymous <codonpair>\n
    <SFS line>\n
    Intron for <codonpair>\n
    <SFS line>\n
Usage:
    python make_codon_pair_SFS_from_SNP_paired_allele_counts_fixed.py \
        -i <input_pair_counts> -o <output_sfs_file> -n <target_nc> [-f] [-e seed]
"""

import sys
import argparse
import numpy as np
import random


def readfile(infilename):
    """Parse paired SNP allele counts file.

    Accepts either:
    1) A header line containing these required column names in any order:
        codonpair  Synnc  SynCount  Intronnc  IntronCount
    2) Headerless 5-column rows in the positional order above.

    Returns:
        Nsnpcounts: dict[codonpair][nc][count] -> multiplicity (intron)
        Ssnpcounts: dict[codonpair][nc][count] -> multiplicity (synonymous)
    """
    Nsnpcounts = {}
    Ssnpcounts = {}
    required_cols = ["codonpair", "synnc", "syncount", "intronnc", "introncount"]

    def add_row(parts, col_idx):
        if len(parts) <= max(col_idx.values()):
            return

        codonpair = parts[col_idx["codonpair"]]
        try:
            syn_nc = int(parts[col_idx["synnc"]])
            syn_cnt = int(parts[col_idx["syncount"]])
            intron_nc = int(parts[col_idx["intronnc"]])
            intron_cnt = int(parts[col_idx["introncount"]])
        except ValueError:
            # Header-like or malformed row
            return

        if codonpair not in Nsnpcounts:
            Nsnpcounts[codonpair] = {}
        if codonpair not in Ssnpcounts:
            Ssnpcounts[codonpair] = {}

        # Synonymous bucket
        if syn_nc not in Ssnpcounts[codonpair]:
            Ssnpcounts[codonpair][syn_nc] = {}
        if syn_cnt not in Ssnpcounts[codonpair][syn_nc]:
            Ssnpcounts[codonpair][syn_nc][syn_cnt] = 0
        Ssnpcounts[codonpair][syn_nc][syn_cnt] += 1

        # Intron bucket
        if intron_nc not in Nsnpcounts[codonpair]:
            Nsnpcounts[codonpair][intron_nc] = {}
        if intron_cnt not in Nsnpcounts[codonpair][intron_nc]:
            Nsnpcounts[codonpair][intron_nc][intron_cnt] = 0
        Nsnpcounts[codonpair][intron_nc][intron_cnt] += 1

    with open(infilename, "r") as f:
        first_parts = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            first_parts = line.split()
            break

        if first_parts is None:
            return Nsnpcounts, Ssnpcounts

        first_lower = [p.lower() for p in first_parts]
        has_named_header = all(col in first_lower for col in required_cols)
        if has_named_header:
            col_idx = {col: first_lower.index(col) for col in required_cols}
        else:
            # Backward-compatible fallback for headerless positional format.
            col_idx = {
                "codonpair": 0,
                "synnc": 1,
                "syncount": 2,
                "intronnc": 3,
                "introncount": 4,
            }
            add_row(first_parts, col_idx)

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            add_row(parts, col_idx)

    return Nsnpcounts, Ssnpcounts

def check_all_pairs(all_pairs):
    fulllist = ['AAA/AAG', 'AAC/AAT', 'AAG/AAA', 'AAT/AAC', 'ACA/ACC', 'ACA/ACG', 'ACA/ACT', 'ACC/ACA', 'ACC/ACG', 'ACC/ACT', 'ACG/ACA', 'ACG/ACC', 'ACG/ACT', 'ACT/ACA', 'ACT/ACC', 'ACT/ACG', 'AGA/AGG', 'AGA/CGA', 'AGC/AGT', 'AGG/AGA', 'AGG/CGG', 'AGT/AGC', 'ATA/ATC', 'ATA/ATT', 'ATC/ATA', 'ATC/ATT', 'ATT/ATA', 'ATT/ATC', 'CAA/CAG', 'CAC/CAT', 'CAG/CAA', 'CAT/CAC', 'CCA/CCC', 'CCA/CCG', 'CCA/CCT', 'CCC/CCA', 'CCC/CCG', 'CCC/CCT', 'CCG/CCA', 'CCG/CCC', 'CCG/CCT', 'CCT/CCA', 'CCT/CCC', 'CCT/CCG', 'CGA/AGA', 'CGA/CGC', 'CGA/CGG', 'CGA/CGT', 'CGC/CGA', 'CGC/CGG', 'CGC/CGT', 'CGG/AGG', 'CGG/CGA', 'CGG/CGC', 'CGG/CGT', 'CGT/CGA', 'CGT/CGC', 'CGT/CGG', 'CTA/CTC', 'CTA/CTG', 'CTA/CTT', 'CTA/TTA', 'CTC/CTA', 'CTC/CTG', 'CTC/CTT', 'CTG/CTA', 'CTG/CTC', 'CTG/CTT', 'CTG/TTG', 'CTT/CTA', 'CTT/CTC', 'CTT/CTG', 'GAA/GAG', 'GAC/GAT', 'GAG/GAA', 'GAT/GAC', 'GCA/GCC', 'GCA/GCG', 'GCA/GCT', 'GCC/GCA', 'GCC/GCG', 'GCC/GCT', 'GCG/GCA', 'GCG/GCC', 'GCG/GCT', 'GCT/GCA', 'GCT/GCC', 'GCT/GCG', 'GGA/GGC', 'GGA/GGG', 'GGA/GGT', 'GGC/GGA', 'GGC/GGG', 'GGC/GGT', 'GGG/GGA', 'GGG/GGC', 'GGG/GGT', 'GGT/GGA', 'GGT/GGC', 'GGT/GGG', 'GTA/GTC', 'GTA/GTG', 'GTA/GTT', 'GTC/GTA', 'GTC/GTG', 'GTC/GTT', 'GTG/GTA', 'GTG/GTC', 'GTG/GTT', 'GTT/GTA', 'GTT/GTC', 'GTT/GTG', 'TAC/TAT', 'TAT/TAC', 'TCA/TCC', 'TCA/TCG', 'TCA/TCT', 'TCC/TCA', 'TCC/TCG', 'TCC/TCT', 'TCG/TCA', 'TCG/TCC', 'TCG/TCT', 'TCT/TCA', 'TCT/TCC', 'TCT/TCG', 'TGC/TGT', 'TGT/TGC', 'TTA/CTA', 'TTA/TTG', 'TTC/TTT', 'TTG/CTG', 'TTG/TTA', 'TTT/TTC']

    diff = [item for item in fulllist if item not in all_pairs]
    return diff

def make_SFS_with_subsampling(args, snpcounts):
    """
    Build an SFS for each codonpair by subsampling to args.nc using a
    hypergeometric draw, aggregating counts across sites that share (nc, count).

    Args:
        snpcounts: dict[codonpair][snc][sa] -> multiplicity

    Returns:
        dict[codonpair] -> list[int] of length args.nc (bins 0..args.nc-1)
    """
    sfss = {}
    for codonpair in snpcounts:
        # Allocate bins 0..nc inclusive, then drop the fixed-derived bin at the end
        sfss[codonpair] = [0 for _ in range(args.nc + 1)]
        for snc in snpcounts[codonpair]:
            if snc >= args.nc:
                for sa, mult in snpcounts[codonpair][snc].items():
                    # successes=sa, failures=snc-sa, draws=args.nc
                    samples = np.random.hypergeometric(sa, snc - sa, args.nc, mult)
                    for c in samples:
                        sfss[codonpair][int(c)] += 1
    # Drop the last bin (fixed derived)
    for codonpair in list(sfss.keys()):
        sfss[codonpair] = sfss[codonpair][:-1]
    return sfss


def parsecommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="foldit", action="store_true", help="fold the resulting SFS")
    parser.add_argument("-i", dest="infilename", required=True, type=str, help="input file path/name")
    parser.add_argument("-o", dest="sfsfilename", required=True, type=str, help="output SFS file path")
    parser.add_argument("-n", dest="nc", required=True, type=int, help="# chromosomes to subsample")
    parser.add_argument("-e", dest="seed", type=int, default=1, help="random seed for subsampling")
    args = parser.parse_args(sys.argv[1:])
    args.commandstring = " ".join(sys.argv[1:])
    return args


def fold_sfs(sfs, nc):
    """Fold an unfolded SFS (length nc) into minor-allele bins.

    Preserves singletons at bin 0. Assumes sfs excludes the fixed-derived bin.
    Returns folded SFS of length 1+floor(nc/2).
    """
    nc_is_even = (nc % 2 == 0)
    if nc_is_even:
        fsfs = [0] + [sfs[j] + sfs[nc - j] for j in range(1, nc // 2)] + [sfs[nc // 2]]
    else:
        fsfs = [0] + [sfs[j] + sfs[nc - j] for j in range(1, 1 + nc // 2)]
    fsfs[0] = sfs[0]
    return fsfs


def run(args):
    random.seed(args.seed)
    np.random.seed(args.seed)

    Nsnpcounts, Ssnpcounts = readfile(args.infilename)

    Nsfss = make_SFS_with_subsampling(args, Nsnpcounts)
    Ssfss = make_SFS_with_subsampling(args, Ssnpcounts)

    # Use the union of codon pairs to be robust
    all_pairs = sorted(set(Nsfss.keys()) | set(Ssfss.keys()))
    missingpairs = check_all_pairs(all_pairs)
    if len(missingpairs) > 0:
        print ("Missing Codon Pairs: ",missingpairs)
    else:
        print("All Codon Pairs Found")
    with open(args.sfsfilename, 'w') as fout:
        for codonpair in all_pairs:
            syn_sfs = Ssfss.get(codonpair, [0] * args.nc)
            intron_sfs = Nsfss.get(codonpair, [0] * args.nc)

            if args.foldit:
                syn_sfs = fold_sfs(syn_sfs, args.nc)
                intron_sfs = fold_sfs(intron_sfs, args.nc)

            # Original ordering: Synonymous first, then Intron
            fout.write(f"Synonymous {codonpair}\n")
            fout.write(" ".join(map(str, syn_sfs)) + "\n")
            fout.write(f"Intron for {codonpair}\n")
            fout.write(" ".join(map(str, intron_sfs)) + "\n")


if __name__ == '__main__':
    args = parsecommandline()
    run(args)
