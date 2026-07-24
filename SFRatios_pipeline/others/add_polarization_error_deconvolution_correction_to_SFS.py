"""
Docstring for add_polarization_error_deconvolution_correction_to_SFS

gemini, codex and chatgpt all wrote function to do the deconvolution of polarization errors 
the chatgpt code did not work but the other two give very similar values

it took a while but got gemini and codex to generate functions that seem to do the appropriate deconvolution of errors, with swapping between reciprocal SFSs
	codex wrote this for the paper:
	
	"To approximate an error‑free site frequency spectrum (SFS) from data with possible mispolarization, we assumed a uniform flip rate r (< 0.5) and treated each mirror pair of derived‑allele count bins (i, N−i) as a 2×2 linear mixing between the forward and reverse spectra (fsfs, rsfs). Specifically, for i=1..N−1 with j=N−i, the observed counts satisfy F_obs[i] = (1−r)·F_true[i] + r·R_true[j] and R_obs[j] = (1−r)·R_true[j] + r·F_true[i] (and symmetrically for the pair R_obs[i], F_obs[j]). We inverted this system independently for each mirror pair to obtain deconvolved expectations; for example, F_hat[i] = ((1−r)·F_obs[i] − r·R_obs[j])/(1−2r), with the analogous expressions for R_hat[j], R_hat[i], and F_hat[j]. For even N, the central bin (i = N/2) couples only across fsfs and rsfs; bin 0 (invariant) is passed through unchanged and the fixed bin N is not represented in our arrays. We truncated small negative expectations to zero and converted to integer pseudo‑SFS via multinomial sampling so that each SFS preserves its original number of polymorphic sites (bins 1..N−1). This deconvolution yields a pseudo error‑free SFS under the assumed r, which we used to assess the impact of polarization errors.	"

usage: add_polarization_error_deconvolution_correction_to_SFS.py [-h] -i INFILENAME -o OUTFILENAME [-e ERRORRATE] [-s {g,c,a}]
                                                                 [-z SEED]

options:
  -h, --help      show this help message and exit
  -i INFILENAME
  -o OUTFILENAME
  -e ERRORRATE
  -s {g,c,a}      g for gemini function, c for codex function, a for both
  -z SEED         Random number seed for reproducibility

This only makes sense for unfolded SFSs
"""
import sys
import numpy as np
import argparse

# Comprehensive list of all codon pairs for matching
allpairs = ['AAA/AAG', 'AAC/AAT', 'AAG/AAA', 'AAT/AAC', 'ACA/ACC', 'ACA/ACG', 'ACA/ACT', 'ACC/ACA', 'ACC/ACG', 'ACC/ACT', 'ACG/ACA', 'ACG/ACC', 'ACG/ACT', 'ACT/ACA', 'ACT/ACC', 'ACT/ACG', 'AGA/AGG', 'AGA/CGA', 'AGC/AGT', 'AGG/AGA', 'AGG/CGG', 'AGT/AGC', 'ATA/ATC', 'ATA/ATT', 'ATC/ATA', 'ATC/ATT', 'ATT/ATA', 'ATT/ATC', 'CAA/CAG', 'CAC/CAT', 'CAG/CAA', 'CAT/CAC', 'CCA/CCC', 'CCA/CCG', 'CCA/CCT', 'CCC/CCA', 'CCC/CCG', 'CCC/CCT', 'CCG/CCA', 'CCG/CCC', 'CCG/CCT', 'CCT/CCA', 'CCT/CCC', 'CCT/CCG', 'CGA/AGA', 'CGA/CGC', 'CGA/CGG', 'CGA/CGT', 'CGC/CGA', 'CGC/CGG', 'CGC/CGT', 'CGG/AGG', 'CGG/CGA', 'CGG/CGC', 'CGG/CGT', 'CGT/CGA', 'CGT/CGC', 'CGT/CGG', 'CTA/CTC', 'CTA/CTG', 'CTA/CTT', 'CTA/TTA', 'CTC/CTA', 'CTC/CTG', 'CTC/CTT', 'CTG/CTA', 'CTG/CTC', 'CTG/CTT', 'CTG/TTG', 'CTT/CTA', 'CTT/CTC', 'CTT/CTG', 'GAA/GAG', 'GAC/GAT', 'GAG/GAA', 'GAT/GAC', 'GCA/GCC', 'GCA/GCG', 'GCA/GCT', 'GCC/GCA', 'GCC/GCG', 'GCC/GCT', 'GCG/GCA', 'GCG/GCC', 'GCG/GCT', 'GCT/GCA', 'GCT/GCC', 'GCT/GCG', 'GGA/GGC', 'GGA/GGG', 'GGA/GGT', 'GGC/GGA', 'GGC/GGG', 'GGC/GGT', 'GGG/GGA', 'GGG/GGC', 'GGG/GGT', 'GGT/GGA', 'GGT/GGC', 'GGT/GGG', 'GTA/GTC', 'GTA/GTG', 'GTA/GTT', 'GTC/GTA', 'GTC/GTG', 'GTC/GTT', 'GTG/GTA', 'GTG/GTC', 'GTG/GTT', 'GTT/GTA', 'GTT/GTC', 'GTT/GTG', 'TAC/TAT', 'TAT/TAC', 'TCA/TCC', 'TCA/TCG', 'TCA/TCT', 'TCC/TCA', 'TCC/TCG', 'TCC/TCT', 'TCG/TCA', 'TCG/TCC', 'TCG/TCT', 'TCT/TCA', 'TCT/TCC', 'TCT/TCG', 'TGC/TGT', 'TGT/TGC', 'TTA/CTA', 'TTA/TTG', 'TTC/TTT', 'TTG/CTG', 'TTG/TTA', 'TTT/TTC']

class codonpairdata():
    def __init__(self, flabel=None, fsfs=None, rlabel=None, rsfs=None):
        self.flabel = flabel
        self.fsfs = fsfs
        self.rlabel = rlabel
        self.rsfs = rsfs
        # Dictionaries to store results from multiple deconvolution methods
        self.simulated_f = {} 
        self.simulated_r = {}

    def __call__(self, flabel=None, fsfs=None, rlabel=None, rsfs=None):
        if flabel is not None: self.flabel = flabel
        if fsfs is not None:   self.fsfs = fsfs
        if rlabel is not None: self.rlabel = rlabel
        if rsfs is not None:   self.rsfs = rsfs
        return self

# --- DECONVOLUTION FUNCTIONS ---
def gemini_deconvolve_cross_sfs(fsfs, rsfs, error_rate=0.05):
    f_obs, r_obs = np.array(fsfs, dtype=float), np.array(rsfs, dtype=float)
    n, r = len(f_obs), error_rate
    denom = 1 - 2 * r
    f_hat, r_hat = np.zeros(n), np.zeros(n)
    for i in range(1, n):
        j = n - i
        f_hat[i] = ((1 - r) * f_obs[i] - r * r_obs[j]) / denom
        r_hat[j] = ((1 - r) * r_obs[j] - r * f_obs[i]) / denom
    f_hat, r_hat = np.maximum(f_hat, 0), np.maximum(r_hat, 0)
    def to_integers(original_total, inferred_counts):
        total_inferred = inferred_counts.sum()
        if total_inferred == 0: return np.zeros_like(inferred_counts, dtype=int)
        return np.random.multinomial(int(original_total), inferred_counts / total_inferred)
    return to_integers(f_obs[1:].sum(), f_hat).tolist(), to_integers(r_obs[1:].sum(), r_hat).tolist()

def codex_deconvolve_cross_sfs(fsfs, rsfs, errorrate=0.05):
    f_obs, r_obs = np.asarray(fsfs, dtype=float), np.asarray(rsfs, dtype=float)
    N, r = len(f_obs), float(errorrate)
    denom = 1.0 - 2.0 * r
    f_hat, r_hat = np.zeros_like(f_obs), np.zeros_like(r_obs)
    idx, mir = np.arange(1, N), N - np.arange(1, N)
    a, b = 1.0 - r, r
    Fi = (a * f_obs[idx] - b * r_obs[mir]) / denom
    Ri = (a * r_obs[idx] - b * f_obs[mir]) / denom
    f_hat[idx], r_hat[idx] = np.maximum(Fi, 0.0), np.maximum(Ri, 0.0)
    def to_pseudo(original, inferred):
        total = int(round(original[1:N].sum()))
        if inferred[1:N].sum() <= 0 or total <= 0: return np.zeros_like(original, dtype=int)
        out = np.zeros_like(original, dtype=int)
        out[1:N] = np.random.multinomial(total, inferred[1:N] / inferred[1:N].sum())
        return out
    return to_pseudo(f_obs, f_hat).tolist(), to_pseudo(r_obs, r_hat).tolist()

def main(args):
    # Set the random seed for NumPy if provided
    if args.seed is not None:
        np.random.seed(args.seed)

    with open(args.infilename) as f:
        ls = f.readlines()
    labels, sfss = ls[::2], [[int(x) for x in line.split()] for line in ls[1::2]]
    
    introncpairdic, syncpairdic = {}, {}
    for i, label in enumerate(labels):
        match = next((s for s in allpairs if s in label), None)
        clist = sorted([match[0:3], match[4:]])
        key = f"{clist[0]}/{clist[1]}"
        dic = introncpairdic if "Intron" in label else syncpairdic
        if key not in dic: dic[key] = codonpairdata(flabel=label, fsfs=sfss[i])
        else: dic[key](rlabel=label, rsfs=sfss[i])

    methods_map = {'g': ('Gemi', gemini_deconvolve_cross_sfs), 
                   'c': ('Code', codex_deconvolve_cross_sfs)}
    
    active_keys = ['g', 'c'] if args.aiscripter == 'a' else [args.aiscripter]

    # Process and print DEBUG table for AAA/AAG
    debug_key = 'AAA/AAG'
    if debug_key in syncpairdic:
        print(f"\n{'='*120}\nDEBUG: BIN-BY-BIN COMPARISON FOR PAIR: {debug_key}\n{'='*120}")
        for dtype, dic in [("SYNONYMOUS", syncpairdic), ("INTRON", introncpairdic)]:
            data = dic[debug_key]
            for mk in active_keys:
                data.simulated_f[mk], data.simulated_r[mk] = methods_map[mk][1](data.fsfs, data.rsfs, args.errorrate)
            
            header = f"Bin | Real_F | Real_R"
            for mk in active_keys: header += f" | {methods_map[mk][0]}_F | {methods_map[mk][0]}_R"
            print(f"\n--- {dtype} ---\n{header}\n{'-'*len(header)}")
            
            for i in range(min(160, len(data.fsfs))):
                row = f"{i:<3} | {data.fsfs[i]:<6} | {data.rsfs[i]:<6}"
                for mk in active_keys: row += f" | {data.simulated_f[mk][i]:<6} | {data.simulated_r[mk][i]:<6}"
                print(row)

    # File Writing (All results)
    with open(args.outfilename, 'w') as fo:
        for mk in active_keys:
            if len(active_keys) > 1:
                fo.write(f"# Method: {methods_map[mk][0]}\n")
            for k in sorted(syncpairdic.keys()):
                ds = syncpairdic[k]
                di = introncpairdic[k]
                if mk not in ds.simulated_f:
                    ds.simulated_f[mk], ds.simulated_r[mk] = methods_map[mk][1](ds.fsfs, ds.rsfs, args.errorrate)
                if mk not in di.simulated_f:
                    di.simulated_f[mk], di.simulated_r[mk] = methods_map[mk][1](di.fsfs, di.rsfs, args.errorrate)
                fo.write(f"{ds.flabel}{' '.join(map(str, ds.simulated_f[mk]))}\n")
                fo.write(f"{di.flabel}{' '.join(map(str, di.simulated_f[mk]))}\n")
                fo.write(f"{ds.rlabel}{' '.join(map(str, ds.simulated_r[mk]))}\n")
                fo.write(f"{di.rlabel}{' '.join(map(str, di.simulated_r[mk]))}\n")
                # for dic in [syncpairdic, introncpairdic]:
                #     d = dic[k]
                #     if mk not in d.simulated_f:
                #         d.simulated_f[mk], d.simulated_r[mk] = methods_map[mk][1](d.fsfs, d.rsfs, args.errorrate)
                #     fo.write(f"{d.flabel}{' '.join(map(str, d.simulated_f[mk]))}\n")
                #     fo.write(f"{d.rlabel}{' '.join(map(str, d.simulated_r[mk]))}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infilename", required=True)
    parser.add_argument("-o", dest="outfilename", required=True)
    parser.add_argument("-e", dest="errorrate", default=0.05, type=float)
    parser.add_argument("-s", dest="aiscripter", default='a', choices=['g', 'c', 'a'], help=" g for gemini function, c for codex function, a for both ")
    parser.add_argument("-z", dest="seed", type=int, default=None, help="Random number seed for reproducibility")
    main(parser.parse_args())

