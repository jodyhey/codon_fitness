"""
make 120 folders if not already present 
in each one put a bootstrap sample of  from a synonymous intron paired snp file 
run make_codon_pair_SFS_from_SNP_paired_allele_counts_fixed.py  on each one in their respective folder 

"""

import os
import random
import os.path as op
import subprocess 
import argparse 

def makebootstrapsample(source,dest):
    fi = open(source,"r")
    fils  = fi.readlines()
    fi.close()
    bootls = [fils[0]] + random.choices(fils[1:], k=len(fils[1:]))
    assert len(bootls) == len(fils)
    fo = open(dest,"w")
    for l in bootls:
        fo.write(l)
    fo.close()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract parameter values from SFRatios output files',
        formatter_class=argparse.RawDescriptionHelpFormatter )
    
    parser.add_argument('-i', dest='sourcefile', required=True,
                       help='Path to source paired count file')
    parser.add_argument('-b', dest='pathtobootstrapfolder', required=True,
                       help='path to folder to hold bootstrap folders')
    parser.add_argument('-n', dest='nbsamples',type=int, required=True,
                       help='number of bootstraps')
    
    
    return parser.parse_args()

def main():

    args = parse_arguments()
    curdir = os.getcwd()
    bootstrapsourcefilename  = args.sourcefile
    # for i in range(args.nbsamples):
    for i in [119]:
        npname = op.join(args.pathtobootstrapfolder,"b{}".format(i+1))
        if not op.exists(npname):
            os.makedirs(npname)
        os.chdir(npname)
        bootstrapcountfilename = "b{}_pairs.txt".format(i+1)
        makebootstrapsample(bootstrapsourcefilename,bootstrapcountfilename)
        assert op.exists(bootstrapcountfilename)
        bootstrapSFSfilename = "b{}_SFSs.txt".format(i+1)
        command = [
                "python", 
                "/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/scripts/make_codon_pair_SFS_from_SNP_paired_allele_counts_fixed.py", 
                "-n", "160", 
                "-e", "1", 
                "-i", bootstrapcountfilename, 
                "-o", bootstrapSFSfilename
                ]

        # check=True will raise an exception if the script fails (non-zero exit code)
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        # 2. Check if the specific string is in the captured stdout
        if "All Codon Pairs Found" in result.stdout:
            print("Success: Script confirmed all codon pairs were found for b{}".format(i+1))
        else:
            print("Warning: Missing codon pairs or script failed for b{}".format(i+1))
            # You can print the actual output to see what went wrong
            # print(f"Subprocess Output: {result.stdout}")
        os.chdir(curdir)




if __name__ == '__main__':
    main()
