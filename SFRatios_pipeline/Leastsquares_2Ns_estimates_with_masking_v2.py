"""
========================================================================================
Leastsquares_2Ns_estimates_with_masking_v2.py

Estimate codon fitness values from synonymous codon pair 2Ns estimates using least
squares optimization.

Based on: Leastsquares_2Ns_estimates_with_masking.py
Modified by: claude.ai
Date: 2025

========================================================================================
PURPOSE:
    This script fits a codon fitness model to empirical 2Ns (selection coefficient)
    estimates for synonymous codon pairs. For each amino acid, it estimates the
    relative fitness of each synonymous codon by minimizing the squared differences
    between observed and predicted 2Ns values.

========================================================================================
INPUT FORMATS:

    Two input formats are supported (specified with -f flag):

    1. 'original' format (default) - Output from Codonpair_2Ns_estimation.py:
        Header: Codon_pair	AA	maxbin	Lklhd	RMSE	invRMSE	ThetaRatio	2Ns
        Example: TTT>TTC	F	75	890.661	0.0210	5.9953	0.130	0.228
                 CTT>CTC	L	75	818.340	0.0219	5.7341	0.151	0.120

        - Codon pairs are explicitly provided in "C1>C2" format
        - 2Ns value is in the last column

    2. 'summary' format - Output from summarize_multiple_SFRatios_runs.py:
        Header: filename	NeutSNPcount	SelSNPcount	2Ns	AIC	likelihood	thetaratio
        Example: Dsimulans_imputed_Synonymous_Aga_Cga_Qratio_fixed2Ns_nc174_estimates.out
                 602	539	-4.8605	-77.179	40.59	0.7996

        - Codon pairs are extracted from filenames using pattern matching
        - Supports flexible filename formats (e.g., Aga_Cga, aaA_aaG)
        - Codons are automatically converted to uppercase
        - 2Ns value is in the '2Ns' column

========================================================================================
METHOD:

    The key assumption is that selection on synonymous codons can be modeled as
    differences in codon fitness values.

    1. DATA STRUCTURE:
       - For each amino acid, create a matrix where rows and columns represent codons
       - Matrix[i,j] contains the 2Ns estimate for transitioning from codon i to codon j
       - We assume 2Ns(j→i) = -2Ns(i→j) (directional symmetry)

    2. FITNESS MODEL:
       - Each codon has a fitness value p[i]
       - The 2Ns for codon i→j equals: 2Ns[i,j] = p[j] - p[i]
       - This implies: matrix[i, j] = p[j] - p[i]
                      matrix[j, i] = p[i] - p[j]

    3. LEAST SQUARES FITTING:
       - Minimize: Σ (observed_2Ns[i,j] - (p[j] - p[i]))²
       - First codon fitness is arbitrarily set to 0 (model is only identifiable up to
         a constant)
       - Uses scipy.optimize.minimize with Nelder-Mead method

    4. MASKING:
       - Optional filtering of extreme 2Ns values (|2Ns| > threshold) using -m flag
       - Masked values are excluded from the least squares calculation
       - Useful for removing outliers or unreliable estimates

    5. STATISTICAL ASSESSMENT:
       - Correlation between observed and predicted 2Ns values
       - Permutation test: shuffle 2Ns values and refit model (n simulations)
       - P-value: proportion of simulated correlations ≥ observed correlation
       - Fisher's combined test across all amino acids

========================================================================================
OUTPUT:

    The output file contains three sections:

    1. DETAILED RESULTS PER AMINO ACID:
       - List of codons sorted by fitness (low to high)
       - Mean fitness, fitness range, and range proportion
       - Least squares fit value
       - Correlation between observed and fitted 2Ns
       - Number of masked/excluded points

    2. SUMMARY TABLE:
       - One row per amino acid
       - Best and worst codon
       - Fitness range
       - Correlation coefficient
       - P-value from permutation test
       - Fisher's combined test p-value

    3. CODON FITNESS TABLE:
       - Tab-delimited table of amino acid, codon, and fitness value
       - Suitable for downstream analysis or plotting

    Optional: Scatter plots (with -p flag) showing observed vs. fitted 2Ns for each AA

========================================================================================
USAGE EXAMPLES:

    # Process summary format with default parameters
    python Leastsquares_2Ns_estimates_with_masking_v2.py \\
        -a Dsimulans_imputed_unfolded_SFRatios_runs.txt \\
        -f summary \\
        -o output.txt \\
        -n 1000

    # Process original format with masking and plots
    python Leastsquares_2Ns_estimates_with_masking_v2.py \\
        -a codon_pair_2Ns_estimates.txt \\
        -f original \\
        -o output.txt \\
        -m 50 \\
        -p \\
        -n 1000

    # Mask extreme values (|2Ns| > 20), 500 simulations, seed 42
    python Leastsquares_2Ns_estimates_with_masking_v2.py \\
        -a input.txt \\
        -f summary \\
        -o output.txt \\
        -m 20 \\
        -n 500 \\
        -s 42

========================================================================================
COMMAND LINE ARGUMENTS:

    Required:
        -a, --infilename    Path to input file with 2Ns estimates
        -o, --outfilename   Path to output results file

    Optional:
        -f, --format        Input format: 'original' or 'summary' (default: original)
        -m, --max2Ns        Mask 2Ns values with |value| > threshold (default: inf)
        -n, --numsim        Number of permutation simulations (default: 100)
        -s, --seed          Random number seed for reproducibility (default: 1)
        -p, --doplots       Generate scatter plots for each amino acid (default: False)

========================================================================================
NOTES:

    - Codon fitness values are relative (first codon = 0) and include a 2Ne scalar
    - The model assumes directional symmetry: 2Ns(A→B) = -2Ns(B→A)
    - For 2-fold degenerate amino acids, correlation is always 1.0 (not simulated)
    - Missing codon pairs are treated as 0 in the matrix
    - Warning messages indicate lines that could not be parsed or validated

========================================================================================
"""
import numpy as np
from scipy.optimize import minimize
import argparse
import sys
import matplotlib.pyplot as plt
from scipy import stats
import math
import statistics
import random
import bisect
import os.path as op
import re

bothdirections = True

bases = ["A","C","G","T"]

allbases = ["A","C","G","T","a","c","g","t"]

nonbases = ["N","-"]

complement = {"A":"T","a":"t","C":"G","c":"g","G":"C","g":"c","T":"A","t":"a","N":"N"}

AminoAcidList = ["CYS","ASP", "SER", "GLN","MET","ASN", "PRO", "LYS","THR", "PHE", "ALA", "GLY", "ILE", "LEU", "HIS", "ARG","TRP","VAL", "GLU", "TYR"]

AminoAcids3to1 = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "LYS": "K",
    "THR": "T",
    "PHE": "F",
    "ALA": "A",
    "GLY": "G",
    "ILE": "I",
    "LEU": "L",
    "HIS": "H",
    "ARG": "R",
    "TRP": "W",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "STOP": "*"
    }

AA3toSynonymousCodons = {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"],
    "STOP": ["TAG", "TGA", "TAA"]}

AminoAcids1to3 = {value.lower(): key for key, value in AminoAcids3to1.items()}
SynonymousCodonstoAA3 = {item: key for key, value_list in AA3toSynonymousCodons.items() for item in value_list}
SynonymousCodonstoAA1 = {key: AminoAcids3to1[value].lower() for key, value in SynonymousCodonstoAA3.items()}

# function to minimize
def LSQdistance(params, datam, mask):
    m = create_matrix_from_params(params)

    # Use the mask to ignore values with |2Ns| > max2Ns
    return np.sum(((m - datam) * mask)**2)

def create_matrix_from_params(p):
    # this is the heart of the method

    # Convert p to a list and insert 0 at the beginning
    p = np.array([0] + list(p)) # arbitrarily set the first codons 2Ns value to 0
    n = len(p)

    # Create an n x n matrix filled with zeros
    matrix = np.zeros((n, n))

    # Fill the matrix according to the specified rules
    for i in range(n):
        for j in range(n):
            if i != j:
                matrix[i, j] = p[j] - p[i]
                matrix[j, i] = p[i] - p[j]

    return matrix

def extract_codon_pair_from_filename(filename):
    """
    Extract codon pair from filename with flexible pattern matching.
    Handles formats like:
    - Dsimulans_imputed_Synonymous_Aga_Cga_Qratio_fixed2Ns_nc174_estimates.out
    - Dsimulans_imputed_Synonymous_aaA_aaG_Qratio_fixed2Ns_nc174_estimates.out

    Returns tuple of (codon1, codon2) in uppercase, or None if not found
    """
    # Pattern to match two codons separated by underscore
    # Codons can be 3 letters (e.g., Aga, Cga) or 2 letters + 1 letter (e.g., aa, A)
    # This pattern looks for: Synonymous_<codon1>_<codon2>_
    pattern = r'Synonymous_([a-zA-Z]{2,3})_([a-zA-Z]{2,3})_'

    match = re.search(pattern, filename)
    if match:
        c1 = match.group(1).upper()
        c2 = match.group(2).upper()
        return c1, c2

    return None, None

def process_file(filename, input_format='original'):
    """
    Builds dictionary to hold the results
    aa2Nsdict is a dictionary with amino acids as keys and a dictionary of codon_pairs as values
    aaCodondict is a dictionary with amino acids as keys and a list of codons as values (sorted alphabetically)

    input_format: 'original' for Codonpair_2Ns_estimation.py output
                  'summary' for summarize_multiple_SFRatios_runs.py output
    """
    global bothdirections

    aa2Nsdict = {}
    aaCodondict = {}
    skipped_lines = []

    # Read the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    if input_format == 'original':
        # Original format processing
        # Find the start of the table
        table_start = 0
        bothdirections = True
        for i, line in enumerate(lines):
            if line.upper().startswith("CODON_PAIR\t"):
                table_start = i + 1
                break

        # Process each line of the table
        for line in lines[table_start:]:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue

            codon_pair = fields[0]
            c1 = codon_pair[:3]
            c2 = codon_pair[4:]

            # Validate codons
            if c1 not in SynonymousCodonstoAA1 or c2 not in SynonymousCodonstoAA1:
                skipped_lines.append(f"Invalid codons: {codon_pair}")
                continue

            aa = SynonymousCodonstoAA1[c1].upper()
            try:
                ns_value = float(fields[-1])  # Convert 2Ns to float
            except ValueError:
                continue  # Skip lines where 2Ns is not a valid float

            # If AA is not in the dictionary, add it
            if aa not in aa2Nsdict:
                aa2Nsdict[aa] = {}
                aaCodondict[aa] = set()

            # Add the codon pair and 2Ns value (convert to format c1>c2)
            codon_pair_formatted = f"{c1}>{c2}"
            aa2Nsdict[aa][codon_pair_formatted] = ns_value
            aaCodondict[aa].add(c1)
            aaCodondict[aa].add(c2)

    elif input_format == 'summary':
        # New format processing from summarize_multiple_SFRatios_runs.py
        # Find header line
        header_idx = 0
        for i, line in enumerate(lines):
            if 'filename' in line.lower() and '2Ns' in line:
                header_idx = i
                break

        # Find column indices
        header = lines[header_idx].strip().split('\t')
        try:
            filename_col = header.index('filename')
            ns_col = header.index('2Ns')
        except ValueError:
            print("Error: Could not find required columns 'filename' and '2Ns' in header", file=sys.stderr)
            sys.exit(1)

        # Process data lines
        for line in lines[header_idx + 1:]:
            fields = line.strip().split('\t')
            if len(fields) <= max(filename_col, ns_col):
                continue

            # Extract codon pair from filename
            filename_val = fields[filename_col]
            c1, c2 = extract_codon_pair_from_filename(filename_val)

            if c1 is None or c2 is None:
                skipped_lines.append(f"Could not extract codon pair from: {filename_val}")
                print(filename_val,c1,c2)
                continue

            # Validate codons
            if c1 not in SynonymousCodonstoAA1 or c2 not in SynonymousCodonstoAA1:
                skipped_lines.append(f"Invalid codons extracted: {c1}, {c2} from {filename_val}")
                continue

            # Get amino acid
            aa = SynonymousCodonstoAA1[c1].upper()

            # Verify both codons code for same amino acid
            aa2 = SynonymousCodonstoAA1[c2].upper()
            if aa != aa2:
                skipped_lines.append(f"Codons {c1} and {c2} code for different amino acids: {aa} vs {aa2}")
                continue

            # Extract 2Ns value
            try:
                ns_value = float(fields[ns_col])
            except ValueError:
                skipped_lines.append(f"Invalid 2Ns value: {fields[ns_col]} for {filename_val}")
                continue

            # Add to dictionaries
            if aa not in aa2Nsdict:
                aa2Nsdict[aa] = {}
                aaCodondict[aa] = set()

            # Add the codon pair and 2Ns value
            codon_pair_formatted = f"{c1}>{c2}"
            aa2Nsdict[aa][codon_pair_formatted] = ns_value
            aaCodondict[aa].add(c1)
            aaCodondict[aa].add(c2)

    else:
        print(f"Error: Unknown input format '{input_format}'", file=sys.stderr)
        sys.exit(1)

    # Convert codon sets to sorted lists
    for aa in aaCodondict:
        aaCodondict[aa] = sorted(list(aaCodondict[aa]))

    # Report skipped lines if any
    if skipped_lines:
        print(f"Warning: Skipped {len(skipped_lines)} lines during processing:", file=sys.stderr)
        for skip_msg in skipped_lines[:10]:  # Show first 10
            print(f"  {skip_msg}", file=sys.stderr)
        if len(skipped_lines) > 10:
            print(f"  ... and {len(skipped_lines) - 10} more", file=sys.stderr)

    return aa2Nsdict, aaCodondict

def build_data_matrix(d2Ns, codons, max2Ns):
    """
    Builds data matrix and creates mask for values exceeding max2Ns threshold
    """
    n = len(codons)
    m = np.zeros((n, n))
    # Create a mask matrix, 1 for values to keep, 0 for values to ignore
    mask = np.ones((n, n))

    if bothdirections:
        for codon_pair in d2Ns:
            s = d2Ns[codon_pair]
            # Handle both ">" and underscore separators
            if '>' in codon_pair:
                c1 = codon_pair[:3]
                c2 = codon_pair[4:]
            else:
                parts = codon_pair.split('_')
                c1 = parts[0]
                c2 = parts[1]

            i = codons.index(c1)
            j = codons.index(c2)
            m[i][j] = s
            # Set mask to 0 for values with |2Ns| > max2Ns
            if abs(s) > max2Ns:
                mask[i][j] = 0
    else:
        for codon_pair in d2Ns:
            s = d2Ns[codon_pair]
            if '>' in codon_pair:
                c1 = codon_pair[:3]
                c2 = codon_pair[4:]
            else:
                parts = codon_pair.split('_')
                c1 = parts[0]
                c2 = parts[1]

            i = codons.index(c1)
            j = codons.index(c2)
            assert m[i][j] == 0 and m[j][i] == 0
            if i<j:
                m[i][j] = s
                m[j][i] = -s
                # Set mask to 0 for values with |2Ns| > max2Ns
                if abs(s) > max2Ns:
                    mask[i][j] = 0
                    mask[j][i] = 0
            else:
                m[i][j] = -s
                m[j][i] = s
                # Set mask to 0 for values with |2Ns| > max2Ns
                if abs(s) > max2Ns:
                    mask[i][j] = 0
                    mask[j][i] = 0

    return n, m, mask

def calc_correlation(datam, fitm, mask):
    # Extract non-diagonal elements and apply mask
    diag_mask = ~np.eye(datam.shape[0], dtype=bool)
    x = datam[diag_mask]
    y = fitm[diag_mask]
    m = mask[diag_mask]

    # Filter out points where mask is 0 (i.e., |2Ns| > max2Ns) or datam is zero
    valid_mask = (m != 0) & (x != 0)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]

    # Check if we have enough points for correlation
    if len(x_valid) < 2:
        return float('nan')  # Return NaN if not enough points

    # Calculate Pearson correlation for valid points
    try:
        correlation, _ = stats.pearsonr(x_valid, y_valid)
        return correlation
    except ValueError:
        # This will catch any other ValueError from pearsonr
        return float('nan')

def create_scatter_plot(datam, fitm, mask, AA,plotfilename):
    # Ensure the input arrays have the same shape
    if datam.shape != fitm.shape:
        raise ValueError("Input arrays must have the same shape")

    # Extract non-diagonal elements and apply mask
    diag_mask = ~np.eye(datam.shape[0], dtype=bool)
    x = datam[diag_mask]
    y = fitm[diag_mask]
    m = mask[diag_mask]

    # Filter out points where mask is 0 (i.e., |2Ns| > max2Ns) or datam is zero
    valid_mask = (m != 0) & (x != 0)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]

    # Get the excluded points to highlight differently
    excluded_mask = (m == 0) & (x != 0)
    x_excluded = x[excluded_mask]
    y_excluded = np.zeros_like(x_excluded)  # Dummy y values for excluded points

    # Create the scatter plot
    plt.figure(figsize=(10, 8))

    # Check if we have enough valid points
    not_enough_data = len(x_valid) < 2

    if not not_enough_data:
        plt.scatter(x_valid, y_valid, alpha=0.5, label='Included data')

        # Calculate Pearson correlation for valid points
        try:
            correlation, _ = stats.pearsonr(x_valid, y_valid)
            corr_text = f'Correlation: {correlation:.4f}'
        except ValueError:
            correlation = float('nan')
            corr_text = 'Correlation: N/A (insufficient data)'
    else:
        correlation = float('nan')
        corr_text = 'Correlation: N/A (insufficient data)'

    if len(x_excluded) > 0:
        plt.scatter(x_excluded, y_excluded, alpha=0.5, color='red', marker='x', label='Excluded (|2Ns| > threshold)')

    # Add labels and title
    plt.xlabel('Data Matrix Values')
    plt.ylabel('Fit Matrix Values')
    plt.title(f"Amino Acid: {AA}")

    # Add a diagonal line for reference
    if len(x_valid) > 0:
        min_val = min(x_valid.min(), y_valid.min())
        max_val = max(x_valid.max(), y_valid.max())
    else:
        min_val, max_val = -1, 1
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y = x')

    # Add correlation text to the plot
    plt.text(0.05, 0.95, corr_text,
             transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top')

    # Add count of excluded points
    plt.text(0.05, 0.90, f'Points excluded: {len(x_excluded)}',
             transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top')

    # Add info about valid points
    plt.text(0.05, 0.85, f'Valid points: {len(x_valid)}',
             transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top')

    # Add legend
    plt.legend()

    # Display the plot
    plt.grid(True)
    # plt.show()
    plt.savefig(plotfilename)

    return correlation

def writeAAresult(outfilename, aa, codons, result, rho,excluded_count, total_count):
    """
    write to outfilename a listing of the codons sorted by fitness
    also return a summary string for the aa
    """
    def sort_and_format_paired_data_with_mean(codons, params):
        # Pair the codons and params
        paired_data = list(zip(codons, params))

        # Sort the paired data based on the params
        sorted_data = sorted(paired_data, key=lambda x: x[1])

        # Calculate the mean of the params and normalize by the mean
        mean_value = statistics.mean(params)
        # temp = params - mean_value
        params -= mean_value
        # Subtract the mean from each value
        sorted_data = [(codon, value - mean_value) for codon, value in sorted_data]
        range_val = (max(params)-min(params))
        rangeprop = range_val/mean_value if mean_value != 0 else float('nan')

        # Prepare the output string and list of codon vals
        output_lines = []
        codon_sel_list = []
        mean_line_written = False
        for codon, param in sorted_data:
            if not mean_line_written and param > 0: # mean_value:
                # Add the mean separator line
                output_lines.append('-' * 30)  # Adjust the number of dashes as needed
                mean_line_written = True
            output_lines.append(f"{codon} {param}")
            codon_sel_list.append(f"{aa}\t{codon}\t{param}\n")

        # Join all lines into a single string
        output_string = '\n'.join(output_lines)


        return output_string, mean_value, rangeprop, range_val, sorted_data,codon_sel_list

    params = list(np.array([0.0] + list(result.x)))  # arbitrarily set the first codons 2Ns value to 0   same as is done in create_matrix_from_params()
    codontable, meanfit, rangeprop, range_val, sortedvals,aa_codon_sel_list = sort_and_format_paired_data_with_mean(codons, params)

    excluded_percent = (excluded_count / total_count) * 100 if total_count > 0 else 0

    # Handle NaN correlation (when not enough valid data points)
    if np.isnan(rho):
        corr_str = "NA (insufficient data)"
    else:
        corr_str = f"{rho:.4f}"

    outf = open(outfilename, "a")
    outf.write("Amino acid: {}  # codons: {}  Mean s: {:.3g} Range: {:.3g} RangeProp: {:.3f} LeastSquaresFit: {:.4f} Correlation: {} Points excluded: {} of {}\n{}\n\n".
            format(aa, len(codons), meanfit, range_val, rangeprop, result.fun, corr_str,
                    excluded_count, total_count, codontable))

    outf.close()

    # Format appropriate summary string
    if np.isnan(rho):
        summ_corr = "NA"
    else:
        summ_corr = f"{rho:.3f}"

    return aa_codon_sel_list,"{}\t{}\t{}\t{}\t{:.3g}\t{}\t({}/{})".format(
        aa, len(codons), sortedvals[-1][0], sortedvals[0][0], range_val, summ_corr,
        excluded_count, total_count)

def fitmodel(args, aa2Nsdict, aaCodondict, real, numsim=None):
    """
    fits the model to the 2Ns values for each aa
    runs on real data or simulated
    if real, runs just once and returns a string of results
    if simulated, does numsim trials and returns a dictionary with lists of simulated correlations
    """
    if real:
        realresultstrings = []
        codon_sel_list = ["AA\tCodon\t2Ns\n"]
    else:
        corrdict = {}
    ntrials = 1 if real else numsim
    for aa in aaCodondict:
        if real == False:
            corrdict[aa] = []
        for ni in range(ntrials):
            if real == False and len(aaCodondict[aa]) == 2:  # don't bother if only 2 fold aa
                continue
            if real:
                d2Ns = aa2Nsdict[aa]  # use real data
            else:  # shuffle the 2Ns values
                vals2Ns = list(aa2Nsdict[aa].values())
                # considered whether shuffling or using new random numbers drawn over the range.   Does not really seem to matter for outcomes.

                random.shuffle(vals2Ns)  # random shuffling
                # vals2Ns = [random.uniform(min(vals2Ns), max(vals2Ns)) for _ in range(len(vals2Ns))]  # new random numbers over the range

                d2Ns = dict(zip(aa2Nsdict[aa].keys(), vals2Ns))

            n, datam, mask = build_data_matrix(d2Ns, aaCodondict[aa], args.max2Ns)

            # Count excluded points (where mask is 0 and not on the diagonal)
            diag_mask = ~np.eye(datam.shape[0], dtype=bool)
            excluded_count = np.sum((mask == 0) & diag_mask)
            total_count = np.sum(diag_mask & (datam != 0))  # Only count non-zero elements

            initial_guess = [1] * (n - 1)
            result = minimize(LSQdistance, initial_guess, args=(datam, mask), method='Nelder-Mead')
            fitm = create_matrix_from_params(result.x)

            if real == False:
                rho = calc_correlation(datam, fitm, mask)
                corrdict[aa].append(rho)
            else:
                rho = calc_correlation(datam, fitm, mask)
                if args.doplots and n > 2:
                    plotfilename = "{}_AA_{}_plot.png".format(op.splitext(args.outfilename)[0],aa)
                    create_scatter_plot(datam, fitm, mask, aa,plotfilename)

            if real:
                aa_codon_sel_list,sumstring = writeAAresult(args.outfilename, aa, aaCodondict[aa], result, rho, excluded_count, total_count)
                realresultstrings.append(sumstring)
                codon_sel_list += aa_codon_sel_list

        if real == False:
            corrdict[aa].sort()

    if real:
        return realresultstrings,codon_sel_list
    else:
        return corrdict

def add_sim_corr_to_sumstrings(moresumstrings, simcorrelations):
    """
    kludgy code that for each amino acid:
        identifies the proportion of simulated correlations higher than the observed
        sticks this value on the end of the summary string for that amino acid
    returns updated summary strings
    """
    newmoresumstrings = []
    for i, s in enumerate(moresumstrings):
        aa = s.split()[0]
        corr_str = s.split()[5]  # Updated index to account for excluded count being added

        # Check if the correlation was successfully calculated
        if corr_str != "NA" and aa in simcorrelations and len(simcorrelations[aa]) > 0:
            try:
                realcorr = float(corr_str)
                if realcorr < 1.0:
                    start_index = bisect.bisect_left(simcorrelations[aa], realcorr)
                    end_index = bisect.bisect_right(simcorrelations[aa], realcorr)
                    # Calculate the middle index of the ties
                    middle_index_within_ties = (start_index + end_index) / 2
                    # Calculate the proportion of values greater than the middle of the ties
                    proportion_higher = (len(simcorrelations[aa]) - middle_index_within_ties) / len(simcorrelations[aa])

                    s += "\t{:.3f}".format(proportion_higher)
                else:
                    s += "\tna"
            except ValueError:
                s += "\tna"  # Handle case where corr_str can't be converted to float
        else:
            s += "\tna"  # Handle case with insufficient data for correlation

        newmoresumstrings.append(s)
    return newmoresumstrings

def doFishersCombinedtest(moresumstrings):
    clist = []
    for s in moresumstrings:
        parts = s.strip().split()
        if len(parts) > 0:  # Make sure we have elements
            v = parts[-1]  # Last element should be the p-value
            if v != "na":
                try:
                    clist.append(float(v))
                except (ValueError, IndexError):
                    # Skip values that can't be converted to float
                    pass

    if not clist:
        return " Fisher's Combined Test: No valid p-values"

    try:
        # Use SciPy to get Fisher's combined statistic and nominal p-value
        statistic, combined_p_value = stats.combine_pvalues(clist, method='fisher')

        # Degrees of freedom for Fisher's method is 2 * k
        df = 2 * len(clist)

        # If the direct p-value underflows to 0.0 or is not finite,
        # compute the log survival function to avoid underflow and format robustly.
        if (combined_p_value == 0.0) or (not np.isfinite(combined_p_value)):
            logp = stats.chi2.logsf(statistic, df)

            if np.isfinite(logp):
                # Convert natural log to base-10 components for a readable scientific string
                log10p = logp / math.log(10.0)
                exp10 = math.floor(log10p)
                # Guard against -inf -> cannot convert to int
                if not np.isfinite(exp10):
                    # Fall back to a conservative bound based on double precision limits
                    smallest = np.nextafter(0, 1)
                    return f" Fisher's Combined Test p<{smallest:.0e}"
                mantissa = 10 ** (log10p - exp10)
                return f" Fisher's Combined Test p≈{mantissa:.15g}e{int(exp10)}"

            # If logp is not finite, fall back to a conservative machine-limit bound
            smallest = np.nextafter(0, 1)
            return f" Fisher's Combined Test p<{smallest:.0e}"

        # Otherwise, report with full double precision in scientific/compact form
        return f" Fisher's Combined Test p={combined_p_value:.16g}"
    except Exception as e:
        # As a last resort, avoid crashing and emit a conservative upper bound
        try:
            smallest = np.nextafter(0, 1)
            return f" Fisher's Combined Test p<{smallest:.0e}"
        except Exception:
            return f" Fisher's Combined Test: Error - {str(e)}"

def run(args):
    np.random.seed(args.seed)
    aa2Nsdict, aaCodondict = process_file(args.infilename, input_format=args.format)

    outf = open(args.outfilename, "w")
    outf.write("Command line: " + args.commandstring + "\n")
    outf.write("Arguments:\n")
    for key, value in vars(args).items():
        outf.write("\t{}: {}\n".format(key, value))
    outf.write("\nSorted (low to high) population selection coefficients, 2Ns, where Fitness = 1-s\n")
    outf.write("\n")
    outf.close()

    sumstrings = ["\nAA\t#codons\tbest_codon\tworst_codon\ts_range\tCorrelation\tExcluded\tCorrSim_p_(prop._sim_higher)"]
    moresumstrings,codon_sel_list = fitmodel(args, aa2Nsdict, aaCodondict, True)  # fit to real data
    simcorrelations = fitmodel(args, aa2Nsdict, aaCodondict, False, numsim=args.numsim)  # run simulations if more than two codons
    moresumstrings = add_sim_corr_to_sumstrings(moresumstrings, simcorrelations)
    combinedresult = doFishersCombinedtest(moresumstrings)

    sumstrings += moresumstrings
    outf = open(args.outfilename, "a")
    outf.write("\n\nSummary Table " + combinedresult + "\n")
    outf.write("{}".format("\n".join(sumstrings)))
    outf.write("\n")
    outf.write("".join(codon_sel_list))
    outf.close()

def parsecommandline():
    parser = argparse.ArgumentParser(
        description="Fit codon fitness model to 2Ns estimates",
        epilog="""
Input formats:
  original: Output from Codonpair_2Ns_estimation.py (format: Codon_pair\\tAA\\t...\\t2Ns)
  summary:  Output from summarize_multiple_SFRatios_runs.py (format: filename\\t...\\t2Ns)
        """
    )
    parser.add_argument("-a", dest="infilename", required=True, type=str, help="Path for 2Ns estimate file")
    parser.add_argument("-f", dest="format", type=str, default="original", choices=["original", "summary"],
                        help="Input file format: 'original' (default) or 'summary'")
    parser.add_argument("-n", dest="numsim", type=int, default=100, help="# of simulations for assessing correlation when more than two codons per AA, default=100")
    parser.add_argument("-o", dest="outfilename", type=str, required=True, help="path to results file")
    parser.add_argument("-s", dest="seed", type=int, default=1, help="random number seed")
    parser.add_argument("-p", dest="doplots", action="store_true", default=False, help="do a plot for each amino acid, default is not to")
    parser.add_argument("-m", dest="max2Ns", type=float, default=float('inf'), help="maximum absolute value for 2Ns to include in analysis")

    args = parser.parse_args(sys.argv[1:])
    args.commandstring = " ".join(sys.argv[1:])
    print(args.commandstring)
    return args

if __name__ == '__main__':
    """
    Modified to handle:
    1. Cases where 2Ns values have an absolute value larger than max2Ns
    2. Two input formats: 'original' (Codonpair_2Ns_estimation.py) and 'summary' (summarize_multiple_SFRatios_runs.py)
    """
    args = parsecommandline()
    run(args)
