#!/usr/bin/env python3
"""
Compare rooting quality between whole-chromosome and gene-by-gene alignments.
"""

import pickle
import pandas as pd
from pathlib import Path

def compare_rooting_quality():
    base_path = Path("/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs")
    
    # Load SNP data
    snp_file = base_path / "ZI_SynSNPs_transcript_pos.txt"
    df = pd.read_csv(snp_file, sep='\t')
    
    # Load improved rooting lookup
    improved_rooting_file = base_path / "gene_by_gene_alignment/improved_rooting_lookup.pkl"
    
    if not improved_rooting_file.exists():
        print("Error: Run gene-by-gene alignment first!")
        return
    
    with open(improved_rooting_file, 'rb') as f:
        improved_rooting = pickle.load(f)
    
    # Load original alignment (from chromosome-wide approach)
    from Bio import SeqIO
    original_fasta = Path("/mnt/d/genemod/better_dNdS_models/drosophila/Dsim_resources/individual_chroms/dm6_aligned_to_drosim2_ASM75419v3.fa")
    
    # Read original alignment
    dmel_seqs = {}
    for record in SeqIO.parse(original_fasta, "fasta"):
        chrom = record.id.split('_')[0]
        dmel_seqs[chrom] = str(record.seq)
    
    # Compare rooting for overlapping SNPs
    original_stats = {'ref': 0, 'alt': 0, 'n': 0, 'other': 0}
    improved_stats = {'ref': 0, 'alt': 0, 'n': 0, 'other': 0}
    both_available = 0
    
    for idx, row in df.iterrows():
        chrom = row['Chrom'].replace('chr', '')
        pos = int(row['POS']) - 1  # 0-based
        ref = row['REF']
        alt = row['ALT']
        
        snp_key = f"{row['Chrom']}_{row['POS']}"
        
        # Get original rooting
        if chrom in dmel_seqs and pos < len(dmel_seqs[chrom]):
            original_base = dmel_seqs[chrom][pos].upper()
        else:
            original_base = 'N'
        
        # Get improved rooting
        if snp_key in improved_rooting:
            improved_base = improved_rooting[snp_key]['dmel_base']
            both_available += 1
            
            # Count original
            if original_base == 'N':
                original_stats['n'] += 1
            elif original_base == ref:
                original_stats['ref'] += 1
            elif original_base == alt:
                original_stats['alt'] += 1
            else:
                original_stats['other'] += 1
            
            # Count improved
            if improved_base == 'N':
                improved_stats['n'] += 1
            elif improved_base == ref:
                improved_stats['ref'] += 1
            elif improved_base == alt:
                improved_stats['alt'] += 1
            else:
                improved_stats['other'] += 1
    
    print("=== ROOTING QUALITY COMPARISON ===")
    print(f"SNPs with both alignments available: {both_available:,}")
    print()
    
    print("ORIGINAL (whole-chromosome) alignment:")
    total_orig = sum(original_stats.values())
    for category, count in original_stats.items():
        pct = count / total_orig * 100 if total_orig > 0 else 0
        print(f"  Dmel = {category.upper()}: {count:,} ({pct:.1f}%)")
    
    print()
    print("IMPROVED (gene-by-gene) alignment:")
    total_impr = sum(improved_stats.values())
    for category, count in improved_stats.items():
        pct = count / total_impr * 100 if total_impr > 0 else 0
        print(f"  Dmel = {category.upper()}: {count:,} ({pct:.1f}%)")
    
    print()
    print("=== IMPROVEMENT SUMMARY ===")
    orig_problematic = original_stats['alt'] / total_orig * 100 if total_orig > 0 else 0
    impr_problematic = improved_stats['alt'] / total_impr * 100 if total_impr > 0 else 0
    
    print(f"Original Dmel=ALT rate: {orig_problematic:.1f}%")
    print(f"Improved Dmel=ALT rate: {impr_problematic:.1f}%")
    
    if impr_problematic < orig_problematic:
        improvement = orig_problematic - impr_problematic
        print(f"✅ IMPROVEMENT: {improvement:.1f} percentage point reduction in problematic rooting!")
    else:
        print("⚠️  No improvement detected")
    
    # Count how many SNPs changed from ALT to REF
    changed_to_ref = 0
    stayed_problematic = 0
    
    for idx, row in df.iterrows():
        chrom = row['Chrom'].replace('chr', '')
        pos = int(row['POS']) - 1
        ref = row['REF']
        alt = row['ALT']
        snp_key = f"{row['Chrom']}_{row['POS']}"
        
        if snp_key in improved_rooting and chrom in dmel_seqs and pos < len(dmel_seqs[chrom]):
            original_base = dmel_seqs[chrom][pos].upper()
            improved_base = improved_rooting[snp_key]['dmel_base']
            
            if original_base == alt and improved_base == ref:
                changed_to_ref += 1
            elif original_base == alt and improved_base == alt:
                stayed_problematic += 1
    
    print(f"SNPs fixed (ALT → REF): {changed_to_ref:,}")
    print(f"SNPs still problematic (ALT → ALT): {stayed_problematic:,}")

if __name__ == "__main__":
    compare_rooting_quality()
