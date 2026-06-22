import sys
import argparse
import os
import os.path as op
import subprocess 
import concurrent.futures as cf
from collections import Counter
from itertools import zip_longest
from statistics import fmean

codonpairs = [
                "AAA_AAG","AAC_AAT","AAG_AAA","AAT_AAC","ACA_ACC",
                "ACA_ACG","ACA_ACT","ACC_ACA","ACC_ACG","ACC_ACT",
                "ACG_ACA","ACG_ACC","ACG_ACT","ACT_ACA","ACT_ACC",
                "ACT_ACG","AGA_AGG","AGA_CGA","AGC_AGT","AGG_AGA",
                "AGG_CGG","AGT_AGC","ATA_ATC","ATA_ATT","ATC_ATA",
                "ATC_ATT","ATT_ATA","ATT_ATC","CAA_CAG","CAC_CAT",
                "CAG_CAA","CAT_CAC","CCA_CCC","CCA_CCG","CCA_CCT",
                "CCC_CCA","CCC_CCG","CCC_CCT","CCG_CCA","CCG_CCC",
                "CCG_CCT","CCT_CCA","CCT_CCC","CCT_CCG","CGA_AGA",
                "CGA_CGC","CGA_CGG","CGA_CGT","CGC_CGA","CGC_CGG",
                "CGC_CGT","CGG_AGG","CGG_CGA","CGG_CGC","CGG_CGT",
                "CGT_CGA","CGT_CGC","CGT_CGG","CTA_CTC","CTA_CTG",
                "CTA_CTT","CTA_TTA","CTC_CTA","CTC_CTG","CTC_CTT",
                "CTG_CTA","CTG_CTC","CTG_CTT","CTG_TTG","CTT_CTA",
                "CTT_CTC","CTT_CTG","GAA_GAG","GAC_GAT","GAG_GAA",
                "GAT_GAC","GCA_GCC","GCA_GCG","GCA_GCT","GCC_GCA",
                "GCC_GCG","GCC_GCT","GCG_GCA","GCG_GCC","GCG_GCT",
                "GCT_GCA","GCT_GCC","GCT_GCG","GGA_GGC","GGA_GGG",
                "GGA_GGT","GGC_GGA","GGC_GGG","GGC_GGT","GGG_GGA",
                "GGG_GGC","GGG_GGT","GGT_GGA","GGT_GGC","GGT_GGG",
                "GTA_GTC","GTA_GTG","GTA_GTT","GTC_GTA","GTC_GTG",
                "GTC_GTT","GTG_GTA","GTG_GTC","GTG_GTT","GTT_GTA",
                "GTT_GTC","GTT_GTG","TAC_TAT","TAT_TAC","TCA_TCC",
                "TCA_TCG","TCA_TCT","TCC_TCA","TCC_TCG","TCC_TCT",
                "TCG_TCA","TCG_TCC","TCG_TCT","TCT_TCA","TCT_TCC",
                "TCT_TCG","TGC_TGT","TGT_TGC","TTA_CTA","TTA_TTG",
                "TTC_TTT","TTG_CTG","TTG_TTA","TTT_TTC"
            ]

codons = [
		"AAG","AAA","AAT","AAC","ACG","ACC",
		"ACA","ACT","CGT","CGC","CGA","AGG",
		"AGA","CGG","TCG","AGT","TCA","AGC",
		"TCC","TCT","ATC","ATA","ATT","CAA",
		"CAG","CAT","CAC","CCG","CCA","CCC",
		"CCT","CTG","CTA","TTG","CTC","CTT",
		"TTA","GAG","GAA","GAT","GAC","GCC",
		"GCG","GCA","GCT","GGC","GGT","GGA",
		"GGG","GTC","GTG","GTT","GTA","TAT",
		"TAC","TGC","TGT","TTC","TTT"
		]

def get_unknown_prefix(folder_path, threshold=0.7):
    # 1. Get all filenames
    files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    if not files:
        return None, 0
    
    num_files = len(files)
    prefix_chars = []
    
    # 2. Zip filenames character-by-character
    # zip_longest handles different filename lengths by filling blanks with ''
    for chars in zip_longest(*files, fillvalue=''):
        # Count occurrences of each character at this specific position
        char_counts = Counter(chars)
        
        # Get the most common character and its frequency
        most_common_char, count = char_counts.most_common(1)[0]
        
        # 3. Check if the "majority" matches
        if (count / num_files) >= threshold:
            prefix_chars.append(most_common_char)
        else:
            # Stop as soon as the consensus breaks
            break
            
    prefix = "".join(prefix_chars)
    
    # 4. Count final files matching this prefix
    final_count = sum(1 for f in files if f.startswith(prefix))
    
    return prefix, final_count

def getprimary2Ns(summaryfilename,primary2Ns):
    with open(summaryfilename,'r') as fh:
        ls = fh.readlines()
    for l in ls[1:]:
        x = l.split()
        for cp in primary2Ns:
            if cp in x[0]:
                primary2Ns[cp].append(float(x[3]))
                break
    return primary2Ns
            
def getfitted2Ns(LSfilename,fitted2Ns,estg):
    with open(LSfilename,'r') as fh:
        ls = fh.readlines()
    for li,l in enumerate(ls):
        if "AA	Codon	2Ns" in l:
            nls = ls[li+1:]
            break
    tempestg = {}
    for l in nls:
        x = l.split()
        if len(x) > 1:
            estg[x[1]].append(float(x[2]))
            tempestg[x[1]] = float(x[2])
    for cp in codonpairs:
        fc = cp[0:3]
        tc = cp[4:]
        fitted2Ns[cp].append(tempestg[tc] - tempestg[fc])
    return fitted2Ns,estg

def parsecommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="rootfolder", type=str, help="folder containing subfolders, each with a set of SFRatios runs")
    parser.add_argument("-o", dest="outfilename", type=str, default = None, help="results file (not path),  goes in rootfolder")
    parser.add_argument("-j", dest="jobs", type=int, default=1, help="number of folders to process in parallel")
    parser.add_argument("-m", "--max-folders", dest="maxfolders", type=int, default=None,
                        help="maximum number of folders to process (default: all)")
    # parser.add_argument("-i", dest="infilename", required=True, type=str, help="input file path/name")
    # parser.add_argument("-o", dest="sfsfilename", required=True, type=str, help="output SFS file path")
    # parser.add_argument("-n", dest="nc", required=True, type=int, help="# chromosomes to subsample")
    # parser.add_argument("-e", dest="seed", type=int, default=1, help="random seed for subsampling")
    args = parser.parse_args(sys.argv[1:])
    if args.jobs < 1:
        parser.error("-j must be >= 1")
    if args.maxfolders is not None and args.maxfolders < 1:
        parser.error("-m/--max-folders must be >= 1")
    args.rootfolder = op.abspath(args.rootfolder)
    args.commandstring = " ".join(sys.argv[1:])
    return args

def process_one_folder(full_path, script_dir):
    name = op.basename(full_path)
    prefix, count = get_unknown_prefix(full_path)
    assert count>=134, "{}".format(count)

    prefixwc = "{}*.out".format(prefix)
    summaryfilename = prefix + "SFRatios_summary.txt"
    LSfilename = prefix + "LeastSquares_modelfitting.txt"

    summarize_script = op.join(script_dir, "summarize_multiple_SFRatios_runs.py")
    ls_script = op.join(script_dir, "Leastsquares_2Ns_estimates_with_masking_v2.py")

    command = [sys.executable, summarize_script,
               "-i","{}".format(prefixwc),
               "-o",summaryfilename
               ]
    subprocess.run(command, cwd=full_path, check=True, capture_output=True, text=True)

    command = [sys.executable, ls_script,
               "-f","summary","-n","1000","-s","1",
               "-a",summaryfilename,
               "-o",LSfilename
                ]
    subprocess.run(command, cwd=full_path, check=True, capture_output=True, text=True)

    return {
        "name": name,
        "summaryfilename": op.join(full_path, summaryfilename),
        "LSfilename": op.join(full_path, LSfilename),
    }

def main(args):

    primary2Ns = {cp:[] for cp in codonpairs}
    fitted2Ns  = {cp:[] for cp in codonpairs}
    estg = {c:[] for c in codons}
    folder_paths = []
    for root, dirs, files in os.walk(args.rootfolder):
        for name in dirs:
            folder_paths.append(op.join(root, name))
    folder_paths.sort()
    if args.maxfolders is not None:
        folder_paths = folder_paths[:args.maxfolders]

    numdatasets = len(folder_paths)
    if numdatasets == 0:
        raise RuntimeError("No subfolders found under {}".format(args.rootfolder))

    script_dir = op.dirname(op.abspath(__file__))
    done_results = []
    failures = []

    print("Processing {} folders with -j {}".format(numdatasets, args.jobs))
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for full_path in folder_paths:
            print("queue {}".format(op.basename(full_path)))
            futures[executor.submit(process_one_folder, full_path, script_dir)] = full_path
        for future in cf.as_completed(futures):
            full_path = futures[future]
            name = op.basename(full_path)
            try:
                result = future.result()
                done_results.append(result)
                print("done {}".format(name))
            except Exception as exc:
                failures.append((full_path, exc))
                print("failed {}".format(name))

    if failures:
        lines = ["{}: {}".format(path, exc) for path, exc in failures[:10]]
        raise RuntimeError("One or more folders failed:\n{}".format("\n".join(lines)))

    for result in done_results:
        primary2Ns = getprimary2Ns(result["summaryfilename"],primary2Ns)
        fitted2Ns, estg = getfitted2Ns(result["LSfilename"],fitted2Ns,estg)

    for values in primary2Ns.values():
        values.sort()
        numdatasets = len(values)
    for values in fitted2Ns.values():
        values.sort()
    for values in estg.values():
        values.sort()
    if args.outfilename == None:
        outfilename = op.join(args.rootfolder, "checkbootstrap_{}_datasets.txt".format(numdatasets))
    else:
        outfilename = args.outfilename
    with open(outfilename,'w') as tempout:
        tempout.write("Primary 2Ns\n")
        fmt_list = lambda x: [f"{val:.3g}" for val in x]
        for cp in primary2Ns:
            tempout.write("{} mean: {:.3g} low: {} high: {} full: {}\n".format(
                cp, 
                fmean(primary2Ns[cp]), 
                fmt_list(primary2Ns[cp][0:3]), 
                fmt_list(primary2Ns[cp][-3:]), 
                fmt_list(primary2Ns[cp])))
            
        # for cp in primary2Ns:
        #     tempout.write("{} mean: {} low: {} high: {} full: {}\n".format(cp,fmean(primary2Ns[cp]),primary2Ns[cp][0:3],primary2Ns[cp][-3:],primary2Ns[cp]))
        tempout.write("\nFitted 2Ns\n")
        
        for cp in fitted2Ns:
            tempout.write("{} mean: {:.3g} low: {} high: {} full: {}\n".format(
                cp, 
                fmean(fitted2Ns[cp]), 
                fmt_list(fitted2Ns[cp][0:3]), 
                fmt_list(fitted2Ns[cp][-3:]), 
                fmt_list(fitted2Ns[cp])))
        tempout.write("\nEstimated g\n")
        for c in estg:
            tempout.write("{} mean: {:.3g} low: {} high: {} full: {}\n".format(
                c, 
                fmean(estg[c]), 
                fmt_list(estg[c][0:3]), 
                fmt_list(estg[c][-3:]), 
                fmt_list(estg[c])))   
        lowcii = max(0,round(numdatasets*0.025)-1)
        hicii = min(numdatasets-1,round(numdatasets*0.975))
        tempout.write("\nPrimary 2Ns 95% confidence intervals\n")
        tempout.write("CodonPair\tMean\tLowCI\tHiCI\twidth\n")
        wsum = 0
        wi = 0
        for cp in primary2Ns:
            w = primary2Ns[cp][hicii] - primary2Ns[cp][lowcii]
            wi += 1
            wsum += w 
            tempout.write("{}\t{:.3g}\t{:.3g}\t{:.3g}\t{:.3g}\n".format(cp,fmean(primary2Ns[cp]),primary2Ns[cp][lowcii], primary2Ns[cp][hicii],w))
        tempout.write("Mean width: {:.3g}\n".format(wsum/wi))
        tempout.write("\nFitted 2Ns 95% confidence intervals\n")
        tempout.write("CodonPair\tMean\tLowCI\tHiCI\twidth\n")
        wsum = 0
        wi = 0
        for cp in fitted2Ns:
            w = fitted2Ns[cp][hicii] - fitted2Ns[cp][lowcii]
            wi += 1
            wsum += w 
            tempout.write("{}\t{:.3g}\t{:.3g}\t{:.3g}\t{:.3g}\n".format(cp,fmean(fitted2Ns[cp]),fitted2Ns[cp][lowcii], fitted2Ns[cp][hicii],w))
        tempout.write("Mean width: {:.3g}\n".format(wsum/wi))
        tempout.write("\nEstimated g 95% confidence intervals\n")
        tempout.write("Codon\tMean\tLowCI\tHiCI\twidth\n")
        wsum = 0
        wi = 0        
        for c in estg:
            w = estg[c][hicii] - estg[c][lowcii]
            wi += 1
            wsum += w
            tempout.write("{}\t{:.3g}\t{:.3g}\t{:.3g}\t{:.3g}\n".format(c,fmean(estg[c]),estg[c][lowcii], estg[c][hicii],w))
        tempout.write("Mean width: {:.3g}\n".format(wsum/wi))
 
if __name__ == '__main__':
    args = parsecommandline()
    main(args)
