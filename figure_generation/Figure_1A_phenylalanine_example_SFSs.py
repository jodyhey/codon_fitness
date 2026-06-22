"""
    make SFS plot(s), cumulative and regular,  of SFSs built by downsample_from_vcf_consequence_pickle.py
    or the parallel version of that  program 
usage: plot_SFSs.py [-h] [-a] [-b] [-e NEUTRALLABEL] [-f] [-g] [-k KSTESTS] [-l] [-L LABELS [LABELS ...]] [-m] -o PLOTFILEPATH [-r] -s
                    SFSFILEPATH [-t] [-u XAXISUPPERLIMIT] [-w] [-x XAXISLOWERLIMIT] [-y YAXISLIMIT] [--usesubranges] [--minsubrangecount MINSUBRANGECOUNT]

options:
  -h, --help            show this help message and exit
  -a                    Use alternate plotting: pairs share colors, first gets pattern, second gets solid line
  -b                    Use the text in the file as the plot legend text
  -e NEUTRALLABEL       If using -k, the neutral SFS label
  -f                    Fold the SFSs
  -g                    Add gridlines
  -k KSTESTS            Do Kolmogorov-Smirnov test, with intergenic as neutral, -k 1 one sided -k 2 two sided, default = 0
  -l                    Plot log of SFS, 0's are skipped, does not work with -r
  -L LABELS [LABELS ...]
                        A series of labels, typically the same number as the number SFSs in the sfs file
  -m                    Plot the cumulative SFS, default is regular
  -o PLOTFILEPATH       Path and filename for plot figure
  -r                    Plot the SFS, whether reg or cumulative, proportional to the lowest bin, default is regular
  -s SFSFILEPATH        Path and filename for SFSs
  -t                    Plot ratios of pairs of SFSs (disables -k -l -m -r -x -y options)
  -u XAXISUPPERLIMIT    Highest x axis bin to include, default = None
  -w                    Show the plot on the screen
  -x XAXISLOWERLIMIT    Lowest x axis bin to include (can be 0 to include invariant sites), default = 1
  -y YAXISLIMIT         If '-m ' y axis lower limit, else upper limit
  --usesubranges        Use subrange binning for SFS data
  --minsubrangecount    Minimum count required for each subrange (default: 10)

"""
"""
    make SFS plot(s), cumulative and regular,  of SFSs built by downsample_from_vcf_consequence_pickle.py
    or the parallel version of that  program 
    
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse 
import os.path as op
import itertools
from  scipy.stats import ks_2samp
import math


def readSFS(fn, foldit):
    """
    Reads a file containing headers and SFS data in alternating lines.
    There may be additional nondata lines,  e.g. line 0
    If a nondata line begins with a digit, there will be a problem

    Headers are any non-numeric lines.
    SFS data are space-separated numbers (integer or float) on the next line after each header.
    All SFSs must be the same length
    
    Parameters:
        fn (str): Filename to read
        foldit (bool): Whether to fold the SFS
        
    Returns:
        tuple: (headers, SFSs) where:
            headers (list): List of header strings
            SFSs (list): List of SFS lists (numeric data)
    """
    with open(fn, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]  # Remove empty lines and whitespace
    
    headers = []
    SFSs = []
    sums = []
    ncall = None  # Initialize ncall
    # Check if line starts with a number, skipping first line 
    i = 1
    while i < len(lines):
        # Check if line starts with a number
        if lines[i][0].isdigit():
            break
        i += 1
    while i < len(lines):        
        # Found data
        sfs = list(map(float, lines[i].split()))
        if ncall is None:  # First SFS sets the expected length
            ncall = len(sfs)
        else:
            assert len(sfs) == ncall, f"SFS length mismatch: {len(sfs)} != {ncall}"
        
        if foldit:
            nc = len(sfs)
            if nc % 2 == 0:  # even length
                sfs = [sfs[0]] + [sfs[k] + sfs[nc-k] for k in range(1, nc//2)] + [sfs[nc//2]]
            else:  # odd length
                sfs = [sfs[0]] + [sfs[k] + sfs[nc-k] for k in range(1, 1 + nc//2)]
        sums.append(sum(sfs))        
        SFSs.append(sfs)        
        headers.append(lines[i-1].strip())
        i += 2
    assert len(headers) == len(SFSs), "Mismatch between number of headers and SFS counts"
    i = len(SFSs) -1
    while i >= 0:
        if sum(SFSs[i]) <= 0:
            print("No SNPs in SFS : ",headers[i])
            SFSs.pop(i)
            headers.pop(i)
        i -= 1
    
    return headers, SFSs, sums


def find_subranges1(sfs1, args):
    """
    Identifies subranges in a list where each subrange sums to at least a minimum value.
    
    Parameters:
    sfs1 (list): List of values (counts) where index corresponds directly to bin number (0-indexed)
    args: Object with attributes:
        minsubrangecount (float): Minimum sum required for each subrange
    
    Returns:
    tuple: (xi, subranges) where xi is the starting bin index and subranges is a list of 
           (start_bin, end_bin+1) tuples following Python convention of half-open intervals [start, end)
    """
    subranges = []
    
    if len(sfs1) <= 1:
        return None, []  # Not enough data for subranges
    
    # Find xi - first bin where count < minsubrangecount (starting from bin 1, skip bin 0)
    xi = None
    for i in range(1, len(sfs1)):
        if sfs1[i] < args.minsubrangecount:
            xi = i
            break
    
    if xi is None:
        return None, []  # No subranges needed - all bins meet threshold
    
    # Check if xi + 1 is beyond the array length
    if xi + 1 >= len(sfs1):
        return None, []  # No subranges possible
    
    # Start from xi + 1
    current_index = xi + 1
    
    while current_index < len(sfs1):
        sum_values = 0
        start_bin = current_index
        end_bin = start_bin
        
        # Continue adding bins until sum reaches or exceeds minsubrangecount
        while end_bin < len(sfs1):
            # Add the value at this bin
            sum_values += sfs1[end_bin]
            
            if sum_values >= args.minsubrangecount:
                # We've reached the minimum threshold
                subranges.append((start_bin, end_bin + 1))  # Python-style half-open interval
                # Move to the next bin for the next subrange
                current_index = end_bin + 1
                break
            
            # If we haven't reached minsubrangecount, move to the next bin
            end_bin += 1
        
        # If we didn't find a valid subrange starting from current_index
        # or if we've reached the end of the array, exit the loop
        if end_bin >= len(sfs1) or sum_values < args.minsubrangecount:
            break
    
    xi = subranges[0][0] if subranges else None
    return xi, subranges

def find_subranges2(sfs1, sfs2, args):
    """
    Identifies subranges in first and second SFS where each subrange meets criteria.
    
    Parameters:
    sfs1 (list): first SFS (counts) where index corresponds directly to bin number (0-indexed)
    sfs2 (list): second SFS (counts) where index corresponds directly to bin number (0-indexed)
    args: Object with attributes:
        minsubrangecount (float): Minimum sum required for each subrange in both SFSs
    
    Returns:
    tuple: (xi, subranges) where xi is the starting bin index and subranges is a list of 
           (start_bin, end_bin+1) tuples following Python convention of half-open intervals [start, end)
    """
    subranges = []
    
    # Validate inputs
    if len(sfs1) != len(sfs2):
        print(f"Error: SFS lengths don't match: {len(sfs1)} vs {len(sfs2)}")
        return None, []
    
    if len(sfs1) <= 1:
        return None, []  # Not enough data for subranges
    
    # Starting at bin 1, find first bin where sfs1 < args.minsubrangecount OR sfs2 == 0
    xi = None
    for i in range(1, len(sfs1)):  # Start at bin 1
        if sfs1[i] < args.minsubrangecount or sfs2[i] == 0:
            xi = i
            break
    
    # If no such bin found, no subranges needed
    if xi is None:
        return xi, subranges
        
    # Check if xi is beyond the array length
    if xi >= len(sfs1):
        return xi, subranges  # No subranges possible
    
    # Start from xi
    current_index = xi
    
    while current_index < len(sfs1):
        sfs1_sum = 0
        sfs2_sum = 0
        start_bin = current_index
        end_bin = start_bin
        
        # Continue adding bins until sum reaches criteria
        while end_bin < len(sfs1):
            # Add the values at this bin
            sfs1_sum += sfs1[end_bin]
            sfs2_sum += sfs2[end_bin]
            
            # Check if we meet both criteria
            if sfs1_sum >= args.minsubrangecount and sfs2_sum >= args.minsubrangecount:
                # We've reached the minimum requirements
                # Before adding this subrange, check if all remaining bins have insufficient sfs2
                # If so, extend this subrange to include them all
                extended_end = end_bin + 1
                
                # Check if all remaining bins have insufficient sfs2 to meet minsubrangecount requirement
                if end_bin + 1 < len(sfs2):
                    remaining_sfs2_sum = sum(sfs2[end_bin + 1:])
                    
                    if remaining_sfs2_sum < args.minsubrangecount:
                        # Extend the subrange to include all remaining bins
                        extended_end = len(sfs1)
                
                subranges.append((start_bin, extended_end))  # Python-style half-open interval
                
                # Move to the next bin for the next subrange
                current_index = extended_end
                break
            
            # If we haven't reached criteria, move to the next bin
            end_bin += 1
        
        # If we didn't find a valid subrange starting from current_index
        # or if we've reached the end of the array, exit the loop
        if end_bin >= len(sfs1) or (sfs1_sum < args.minsubrangecount or sfs2_sum < args.minsubrangecount):
            break
    
    return xi, subranges

def convert_sfs_to_subranges(sfs_data, xi, subranges):
    """
    Convert SFS data to subrange-binned data.
    
    Parameters:
        sfs_data (list): Original SFS data
        xi (int): Starting bin index  
        subranges (list): List of (start, end) tuples for subranges
        
    Returns:
        list: Subrange-binned data where each element is the sum within a subrange
    """
    if not subranges or xi is None:
        return list(sfs_data)  # Return original data as list
    
    subrange_data = []
    
    # Add bins before xi (individual bins)
    subrange_data.extend(sfs_data[:xi])
    
    # Add subrange sums
    for start, end in subranges:
        if start < len(sfs_data) and end <= len(sfs_data):
            subrange_sum = sum(sfs_data[start:end])
            subrange_data.append(subrange_sum)
        else:
            # Handle case where subrange extends beyond data
            available_end = min(end, len(sfs_data))
            if start < available_end:
                subrange_sum = sum(sfs_data[start:available_end])
                subrange_data.append(subrange_sum)
    
    return subrange_data

def process_subranges_for_data(data, args):
    """
    Process all SFS data for subranges if enabled.
    
    Parameters:
        data (list): List of SFS data
        args: Command line arguments
        
    Returns:
        tuple: (processed_data, subrange_info) where subrange_info contains xi and subranges
    """
    if not args.usesubranges:
        return data, None
    
    if len(data) == 0:
        return data, None
    
    print(f"Processing subranges with minsubrangecount={args.minsubrangecount}")
    
    # Process each SFS individually to find its subranges
    all_subrange_data = []
    max_subranges = 0
    subrange_infos = []
    
    for i, sfs in enumerate(data):
        xi, subranges = find_subranges1(sfs, args)
        subrange_sfs = convert_sfs_to_subranges(sfs, xi, subranges)
        
        all_subrange_data.append(subrange_sfs)
        subrange_infos.append({'xi': xi, 'subranges': subranges})
        max_subranges = max(max_subranges, len(subrange_sfs))
        
        print(f"SFS {i}: {len(subranges) if subranges else 0} subranges, total length: {len(subrange_sfs)}")
    
    # Pad all SFS data to the same length with NaN (won't be plotted)
    for i in range(len(all_subrange_data)):
        while len(all_subrange_data[i]) < max_subranges:
            all_subrange_data[i].append(float('nan'))
    
    print(f"All SFS data padded to {max_subranges} bins")
    
    # Return the processed data and info from the first SFS for reference
    subrange_info = {'max_length': max_subranges, 'individual_infos': subrange_infos}
    return all_subrange_data, subrange_info


def calculate_proportional_cumulative_sum(numbers):
    cumsum = np.cumsum(numbers)
    return cumsum / cumsum[-1]

def calculate_custom_sum(datai,numbers, xaxislowerlimit, xaxisupperlimit,cumulative, proportional, usesubranges=False):
    # Convert to numpy array to handle NaN values
    numbers = np.array(numbers)
    
    if not usesubranges:
        # Original logic for regular bins
        if xaxislowerlimit > 1:
            numbers = numbers[xaxislowerlimit-1:]    
        if xaxisupperlimit is not None:
            numbers = numbers[:xaxisupperlimit]
        if cumulative:
            # Use nansum for cumulative to handle NaN values
            result = np.nancumsum(numbers)
            if proportional:
                final_sum = result[-1] if len(result) > 0 and np.isfinite(result[-1]) else 0
                if final_sum <= 0.0:
                    print("problem plotting data set ",datai, "result[-1] ", final_sum)
                    print(numbers)
                    return np.zeros_like(result)
                else:
                    return result/final_sum
            else:
                return result
        else:
            # numbers = numbers[1:]
            numbers = numbers[xaxislowerlimit:]
            if proportional and len(numbers) > 0 and np.isfinite(numbers[0]) and numbers[0] > 0:
                return numbers / numbers[0]
            else:
                return numbers
    else:
        # Subrange logic - simpler since subranges are already pre-processed
        if xaxisupperlimit is not None:
            numbers = numbers[:xaxisupperlimit]
        
        if cumulative:
            # Use nansum for cumulative to handle NaN values
            result = np.nancumsum(numbers)
            if proportional:
                final_sum = result[-1] if len(result) > 0 and np.isfinite(result[-1]) else 0
                if final_sum <= 0.0:
                    print("problem plotting data set ",datai, "result[-1] ", final_sum)
                    print(numbers)
                    return np.zeros_like(result)
                else:
                    return result/final_sum
            else:
                return result
        else:
            # For subranges, skip the first bin (position 0) if it exists
            if len(numbers) > 1:
                numbers = numbers[1:]
            if proportional and len(numbers) > 0 and np.isfinite(numbers[0]) and numbers[0] > 0:
                return numbers / numbers[0]
            else:
                return numbers

def calculate_ratios(data, headers, args):
    """
    Calculate ratios of pairs of SFSs, with subrange support.
    
    Parameters:
        data (list): List of SFS data
        headers (list): List of headers corresponding to data
        args: Command line arguments
        
    Returns:
        tuple: (ratio_data, ratio_labels) where:
            ratio_data (list): List of ratio arrays
            ratio_labels (list): List of labels for ratios
    """
    if len(data) % 2 != 0:
        print(f"Error: For ratio plotting (-t), need an even number of SFSs.")
        print(f"Found {len(data)} SFSs, but need pairs for ratios.")
        sys.exit(1)
    
    ratio_data = []
    ratio_labels = []
    max_ratio_length = 0
    
    if args.usesubranges:
        print("Calculating ratios with subrange binning")
        
        # Process each pair with subranges
        for i in range(0, len(data), 2):
            numerator_sfs = data[i]
            denominator_sfs = data[i+1]
            
            # Find subranges for this pair
            xi, subranges = find_subranges2(numerator_sfs, denominator_sfs, args)
            
            if subranges:
                print(f"Pair {i//2 + 1}: Found {len(subranges)} subranges starting at bin {xi}")
                
                # Convert both SFSs to subrange format
                num_subrange = convert_sfs_to_subranges(numerator_sfs, xi, subranges)
                den_subrange = convert_sfs_to_subranges(denominator_sfs, xi, subranges)
                
                # Skip position 0 if it exists
                if len(num_subrange) > 1:
                    num_subrange = num_subrange[1:]
                    den_subrange = den_subrange[1:]
                
                # Calculate ratios
                numerator = np.array(num_subrange)
                denominator = np.array(den_subrange)
            else:
                print(f"Pair {i//2 + 1}: No subranges found, using individual bins")
                # Use original approach without subranges, skip position 0
                numerator = np.array(numerator_sfs[1:])
                denominator = np.array(denominator_sfs[1:])
            
            # Handle division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = numerator / denominator
                # Replace inf and nan with 0 for plotting
                ratio = np.where(np.isfinite(ratio), ratio, 0)
            
            ratio_data.append(ratio)
            max_ratio_length = max(max_ratio_length, len(ratio))
            
            # Create label
            num_label = headers[i] if i < len(headers) else f"SFS_{i}"
            den_label = headers[i+1] if i+1 < len(headers) else f"SFS_{i+1}"
            ratio_labels.append(f"Ratio {num_label}/{den_label}")
    
    else:
        # Original method without subranges
        for i in range(0, len(data), 2):
            numerator = np.array(data[i][1:])  # Skip position 0
            denominator = np.array(data[i+1][1:])  # Skip position 0
            
            # Handle division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = numerator / denominator
                # Replace inf and nan with 0 for plotting
                ratio = np.where(np.isfinite(ratio), ratio, 0)
            
            ratio_data.append(ratio)
            max_ratio_length = max(max_ratio_length, len(ratio))
            
            # Create label
            num_label = headers[i] if i < len(headers) else f"SFS_{i}"
            den_label = headers[i+1] if i+1 < len(headers) else f"SFS_{i+1}"
            ratio_labels.append(f"Ratio {num_label}/{den_label}")
    
    # Pad all ratio arrays to the same length with NaN
    for i in range(len(ratio_data)):
        current_length = len(ratio_data[i])
        if current_length < max_ratio_length:
            padded_array = np.full(max_ratio_length, float('nan'))
            padded_array[:current_length] = ratio_data[i]
            ratio_data[i] = padded_array
    
    if args.usesubranges:
        print(f"All ratio data padded to {max_ratio_length} subrange bins")
    
    return ratio_data, ratio_labels

def plot_ratios(ratio_data, ratio_labels, args, subrange_info=None):
    """
    Plot ratio data with optional subrange support.
    
    Parameters:
        ratio_data (list): List of ratio arrays
        ratio_labels (list): List of labels for ratios
        args: Command line arguments
        subrange_info (dict): Information about subranges if used
    """
    plt.rcParams.update({'font.size': 16*args.fontsizescalar})  # Set default font size
    
    if args.square:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Define colors and line styles
    # Use reordered tab20 to better separate similar hues when many lines are plotted
    cmap = plt.get_cmap('tab20')
    # Reordered to separate similar hues; 20 distinct before cycling
    colors = [cmap(i) for i in [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19]]

    # Favor readability: 6 solids first, then short/dense dash patterns
    line_styles = [
        '-', '-', '-', '-', '-', '-',
        (0, (1, 1)),          # very dense dots (shortest marks)
        (0, (2, 1)),          # very short dash
        (0, (3, 1)),          # short dash
        (0, (3, 1, 1, 1)),    # short dash-dot
        (0, (1, 1, 1, 1, 1, 1)),  # very dense dot-dot-dot
        (0, (2, 1, 1, 1)),    # short dash + short dot
    ]

    for i, (ratio, label) in enumerate(zip(ratio_data, ratio_labels)):
        # Convert to numpy array to handle NaN values
        ratio = np.array(ratio)
        
        # Apply upper limit if specified
        if args.xaxisupperlimit is not None:
            if args.usesubranges:
                # For subranges, limit is based on subrange count
                ratio = ratio[:args.xaxisupperlimit]
            else:
                ratio = ratio[:args.xaxisupperlimit-1]  # -1 because we skipped position 0
        
        # Create x-axis - matplotlib will handle NaN values automatically
        if args.usesubranges:
            # For subranges, x-axis represents subrange indices
            x = np.arange(1, len(ratio) + 1)
            xlabel = 'Subrange Index' if args.altXaxislabel is None else args.altXaxislabel
        else:
            # For regular bins, x-axis represents frequency bins (starting from 1)
            x = np.arange(1, len(ratio) + 1)
            xlabel = xlabel = 'Frequency Bin' if args.altXaxislabel is None else args.altXaxislabel            
        
        color = colors[i % len(colors)]

        line_style = line_styles[i % len(line_styles)]
        
        # Plot - matplotlib automatically handles NaN values by not plotting them
        ax.plot(x, ratio, label=label, color=color, linestyle=line_style, linewidth=2.5)

    # Set the labels and title
    ax.set_xlabel(xlabel, fontsize=16*args.fontsizescalar)
    ax.set_ylabel('Ratio', fontsize=16*args.fontsizescalar)
    
    if args.usesubranges:
        ax.set_title('SFS Ratios (Subrange Binned)', fontsize=15)
    else:
        ax.set_title('SFS Ratios', fontsize=16*args.fontsizescalar)

    # Set x-axis upper limit if specified
    if args.xaxisupperlimit is not None:
        ax.set_xlim(right=args.xaxisupperlimit)

    # Add legend
    if ratio_labels:
        ax.legend(loc='best', fontsize=14*args.fontsizescalar, frameon=True, handlelength=5)

    # Set tick marks
    ax.tick_params(axis='both', which='major', labelsize=16*args.fontsizescalar)
    
    # Set grid lines
    if args.gridlines:
        ax.grid(True, which='major', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)

    # Save the figure
    plt.savefig(args.plotfilepath, dpi=300, bbox_inches='tight')
    
    # Show the plot if necessary
    if args.plot_to_screen:
        plt.show()

def plot_data(data, counts, labels, args,ksresults=None):

    # if len(labels) < len(data):
    #     print("problem,  len(labels) != len(data)")
    #     print(" was -b invoked ?")
    #     print("headers :",labels)
    #     exit()
        # for j in range(len(labels),len(data)):
        #     labels.append("dataset_{}".format(j))
    if ksresults:
        labels = ["{} ({}){}".format(l,round(counts[i]),ksresults[i]) for i,l in enumerate(labels)]
    plt.rcParams.update({'font.size': 16*args.fontsizescalar})  # Set default font size
    if args.square:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Define colors and line styles
    # Use reordered tab20 to better separate similar hues when many lines are plotted
    cmap = plt.get_cmap('tab20')
    # Reordered to separate similar hues; 20 distinct before cycling
    colors = [cmap(i) for i in [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19]]

    if args.alternate_plotline:
        # For alternate plotting: each pair gets same color, different styles
        # Emphasize discriminability: shorter dash patterns first; second in pair stays solid
        pair_line_styles = [
            (0, (1, 1)),          # very dense dots (shortest marks)
            (0, (2, 1)),          # very short dash
            (0, (3, 1)),          # short dash
            (0, (3, 1, 1, 1)),    # short dash-dot
            (0, (1, 1, 1, 1, 1, 1)),  # very dense dot-dot-dot
            (0, (2, 1, 1, 1)),    # short dash + short dot
            (0, (4, 1)),          # short-ish dash
            (0, (5, 1)),          # medium-short dash
        ]
    else:
        # Favor readability: 6 solids first, then short/dense dash patterns
        line_styles = [
            '-', '-', '-', '-', '-', '-',
            (0, (1, 1)),          # very dense dots (shortest marks)
            (0, (2, 1)),          # very short dash
            (0, (3, 1)),          # short dash
            (0, (3, 1, 1, 1)),    # short dash-dot
            (0, (1, 1, 1, 1, 1, 1)),  # very dense dot-dot-dot
            (0, (2, 1, 1, 1)),    # short dash + short dot
        ]
    colors = ["green","green","orange","orange","black"]
    # lines_styles = ["-",(0, (3, 1)),"-",(0, (3, 1)),"-"]
    # Keep original order from input file
    sorted_labels = labels
    sorted_data = data

    for i, (numbers, label) in enumerate(zip(sorted_data, sorted_labels)):
        # Convert to numpy array to handle NaN values properly
        numbers = np.array(numbers)
        
        plotvals = calculate_custom_sum(i,numbers, args.xaxislowerlimit,args.xaxisupperlimit, args.plotcumulative, args.plotproportional, args.usesubranges)
        plotvals = np.array(plotvals)
        
        if args.usesubranges:
            # For subranges, x-axis represents subrange indices
            x = np.arange(1, len(plotvals) + 1)  # Start from 1 for subrange indices
            xlabel = 'Subrange Index' if args.altXaxislabel is None else args.altXaxislabel
        else:
            # For regular bins, use original logic
            x = np.arange(args.xaxislowerlimit, len(plotvals) + args.xaxislowerlimit)
            xlabel = 'Index' if args.altXaxislabel is None else args.altXaxislabel
        
        if args.alternate_plotline:
            # Use the current index for pairing logic (since we're not sorting anymore)
            pair_index = i // 2
            is_second_in_pair = i % 2 == 1
            
            color = colors[pair_index % len(colors)]
            if is_second_in_pair:
                line_style = '-'  # solid line for second in pair
            else:
                line_style = pair_line_styles[pair_index % len(pair_line_styles)]
        else:
            # Original behavior: cycle through colors and line styles in original order
            color = colors[i % len(colors)]
            line_style = line_styles[i % len(line_styles)]
        if i % 2 == 0 :
            line_style = "-"
        else:
            line_style =  (0, (3, 1))
        if args.plotlogsfs:
            # Filter out values where plotvals[i] <= 0 or is NaN
            valid_mask = (plotvals > 0) & np.isfinite(plotvals)
            x_log = x[valid_mask]
            y_log = np.log(plotvals[valid_mask])
            if len(x_log) > 0:  # Only plot if we have valid data
                ax.plot(x_log, y_log, label=label, color=color, linestyle=line_style, linewidth=2.5)
        else:
            # matplotlib automatically handles NaN values by not plotting them
            ax.plot(x, plotvals, label=label, color=color, linestyle=line_style, linewidth=2.5)

    # Set the labels and title
    ax.set_xlabel(xlabel, fontsize=16*args.fontsizescalar)
    if args.plotcumulative:
        if args.plotproportional:
            if args.yaxislimit is not None:
                ax.set_ylim(args.yaxislimit, 1.001)
            ax.set_ylabel('Proportional Cumulative Sum', fontsize=16*args.fontsizescalar)
        else:
            if args.plotlogsfs:
                ax.set_ylabel('Log Cumulative Sum', fontsize=16*args.fontsizescalar)
            else:
                ax.set_ylabel('Cumulative Sum', fontsize=16*args.fontsizescalar)
    else:
        if args.plotproportional:
            if args.yaxislimit is not None:
                ax.set_ylim(0.0,args.yaxislimit)
            ax.set_ylabel("Proportional to lowest frequency bin", fontsize=16*args.fontsizescalar)
        else:
            if args.plotlogsfs:
                ax.set_ylabel('Log Count', fontsize=16*args.fontsizescalar)
            else:
                ax.set_ylabel('Count', fontsize=16*args.fontsizescalar)
            if args.yaxislimit is not None:
                ax.set_ylim(0.0,args.yaxislimit)
            
            if args.usesubranges:
                ax.set_title('SNP Count (Subrange Binned)', fontsize=16*args.fontsizescalar)
            else:
                ax.set_title('SNP Count', fontsize=16*args.fontsizescalar)

    # Add legend with smaller font size and longer line samples
    if sorted_labels:
        ax.legend(loc='best', fontsize=14*args.fontsizescalar, frameon=True, handlelength=5)

    # Set tick marks
    ax.tick_params(axis='both', which='major', labelsize=16*args.fontsizescalar)
    
    # Set grid lines
    if args.gridlines:
        ax.grid(True, which='major', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)

    # Save the figure
    plt.savefig(args.plotfilepath, dpi=300, bbox_inches='tight')
    
    # Show the plot if necessary
    if args.plot_to_screen:
        plt.show()

def kstest(counts1,counts2,alternative='greater'):
    d1 = []
    for ci,count in enumerate(counts1):
        val = ci + 1
        for i in range(round(count)):
            d1.append(val)
    d2 = []
    for ci,count in enumerate(counts2):
        val = ci + 1
        for i in range(round(count)):
            d2.append(val)            
    res = ks_2samp(d1,d2,alternative='greater' if args.KStests==1 else "two-sided")
    
    # Determine direction by comparing cumulative distributions at midpoint
    cum1 = np.cumsum(counts1)
    cum2 = np.cumsum(counts2)
    # Normalize to proportions
    cum1_norm = cum1 / cum1[-1] if cum1[-1] > 0 else cum1
    cum2_norm = cum2 / cum2[-1] if cum2[-1] > 0 else cum2
    
    # Compare at midpoint
    midpoint = len(cum1_norm) // 2
    if midpoint < len(cum1_norm) and midpoint < len(cum2_norm):
        if cum1_norm[midpoint] > cum2_norm[midpoint]:
            direction = "sel>neut@low"  # Selected has more low-frequency variants
        elif cum1_norm[midpoint] < cum2_norm[midpoint]:
            direction = "sel<neut@low"  # Selected has fewer low-frequency variants  
        else:
            direction = "sel≈neut@low"  # Similar at midpoint
    else:
        direction = "unclear"
    
    return res.pvalue, res.statistic, direction

def parsecommandline():
    parser = argparse.ArgumentParser()
    
    # Alphabetized argparse options
    parser.add_argument("-a", dest="alternate_plotline", action="store_true", default=False,
                       help="Use alternate plotting: pairs share colors, first gets pattern, second gets solid line")
    parser.add_argument("-b", dest="useheaderlabels", action="store_true", default=False,
                       help="Use the text in the file as the plot legend text")
    parser.add_argument("-e", dest="neutrallabel", default=None,
                       help="If using -k,and a specific SFS is the neutral control, give the neutral SFS label, if None, assume alternating selected, neutral")
    parser.add_argument("-f", dest="foldit", action="store_true", default=False,
                       help="Fold the SFSs")
    parser.add_argument("-g", dest="gridlines", action="store_true", default=False,
                       help="Add gridlines")
    parser.add_argument("-i", dest="altXaxislabel", type=str,default = None,help = "Alternative X axis label")
    parser.add_argument("-k", dest="KStests", type=int, default=0,
                       help="Do Kolmogorov-Smirnov test, with intergenic as neutral, -k 1 one sided -k 2 two sided, default = 0, does not work with -r")
    parser.add_argument("-l", dest="plotlogsfs", action="store_true", default=False,
                       help="Plot log of SFS, 0's are skipped, does not work with -r")
    parser.add_argument("-L", dest="labels", nargs="+", default=[],
                       help="A series of labels, typically the same number as the number of SFSs in the sfs file")
    parser.add_argument("-m", dest="plotcumulative", action="store_true", default=False,
                       help="Plot the cumulative SFS, default is regular")
    parser.add_argument("-o", dest="plotfilepath", type=str, required=True,
                       help="Path and filename for plot figure")
    parser.add_argument("-q", dest="square", action="store_true", default=False,
                       help="Make a square plot, default is rectangular")
    parser.add_argument("-r", dest="plotproportional", action="store_true", default=False,
                       help="Plot the SFS, whether reg or cumulative, proportional to the lowest bin, default is regular")
    parser.add_argument("-s", dest="sfsfilepath", type=str, required=True,
                       help="Path and filename for SFSs")
    parser.add_argument("-t", dest="plotratios", action="store_true", default=False,
                       help="Plot ratios of pairs of SFSs (disables -k -l -m -r -x -y options)")
    parser.add_argument("-u", dest="xaxisupperlimit", type=int, default=None,
                       help="Highest x axis bin to include, default = None")
    parser.add_argument("-w", dest="plot_to_screen", action="store_true", default=False,
                       help="Show the plot on the screen")
    parser.add_argument("-x", dest="xaxislowerlimit", type=int, default=1,
                       help="Lowest x axis bin to include (can be 0 to include invariant sites), default = 1")
    parser.add_argument("-y", dest="yaxislimit", type=float, default=None,
                       help="If '-m ' y axis lower limit, else upper limit")
    parser.add_argument("-z", dest="fontsizescalar", type=float, default=1.0,
                       help="Multiplier of font sizes, default = 1.0")
    parser.add_argument("--usesubranges", dest="usesubranges", action="store_true", default=False,
                       help="Use subrange binning for SFS data")
    parser.add_argument("--minsubrangecount", dest="minsubrangecount", type=float, default=10,
                       help="Minimum count required for each subrange (default: 10)")

    args = parser.parse_args(sys.argv[1:])   
    args.commandstring = " ".join(sys.argv[1:])
    
    # Validate incompatible options with -t
    if args.plotratios:
        incompatible_options = []
        if args.KStests != 0:
            incompatible_options.append("-k")
        if args.plotlogsfs:
            incompatible_options.append("-l") 
        if args.plotcumulative:
            incompatible_options.append("-m")
        if args.plotproportional:
            incompatible_options.append("-r")
        if args.xaxislowerlimit != 1:
            incompatible_options.append("-x")
        if args.yaxislimit is not None:
            incompatible_options.append("-y")
        
        if incompatible_options:
            print(f"Error: Option -t (plot ratios) is incompatible with: {', '.join(incompatible_options)}")
            sys.exit(1)
    
    return args

if __name__ == '__main__':
    """

    """
    args = parsecommandline()
    
    # Validate subrange parameters
    if args.usesubranges:
        if args.minsubrangecount <= 0:
            print("Error: --minsubrangecount must be positive")
            sys.exit(1)
        if args.minsubrangecount < 1:
            print("Warning: --minsubrangecount < 1 may result in very small subranges")
    
    headers, data, counts = readSFS(args.sfsfilepath,args.foldit)
    
    # Process subranges if enabled (note: folding already done in readSFS)
    if args.usesubranges and not args.plotratios:
        # For regular plotting with subranges
        print("Note: Subranges calculated on folded SFS data" if args.foldit else "Note: Subranges calculated on unfolded SFS data")
        data, subrange_info = process_subranges_for_data(data, args)
        if subrange_info:
            max_length = subrange_info['max_length']
            print(f"All SFS data converted to {max_length} bins (subrange + individual bins)")
        else:
            print("No subranges found, using original data")
    
    if args.plotratios:
        # Calculate and plot ratios (handles subranges internally)
        ratio_data, ratio_labels = calculate_ratios(data, headers, args)
        
        # Use custom labels if provided and -b not used
        if args.labels and not args.useheaderlabels:
            if len(args.labels) >= len(data):
                # Create ratio labels from custom labels
                custom_ratio_labels = []
                for i in range(0, len(data), 2):
                    num_label = args.labels[i] if i < len(args.labels) else f"SFS_{i}"
                    den_label = args.labels[i+1] if i+1 < len(args.labels) else f"SFS_{i+1}"
                    custom_ratio_labels.append(f"Ratio {num_label}/{den_label}")
                ratio_labels = custom_ratio_labels
        
        plot_ratios(ratio_data, ratio_labels, args)
        
    else:
        # Original plotting functionality
        if args.useheaderlabels:
            args.labels = headers
        if args.KStests:
            if args.neutrallabel is None:
                # Alternating pairs mode: second, fourth, sixth... rows are neutral controls
                # First check that we have an even number of datasets
                if len(data) % 2 != 0:
                    print(f"Error: When using -k without -e, expecting alternating pairs.")
                    print(f"Found {len(data)} datasets, but need an even number for pairs.")
                    print(f"Available labels: {args.labels}")
                    sys.exit(1)
                
                print(f"Using alternating pairs mode: {len(data)//2} selected-neutral pairs")
                
                ksresults = []
                for di in range(len(data)):
                    if di % 2 == 0:
                        # This is a selected dataset (even indices: 0, 2, 4, ...)
                        # Its neutral control is at di + 1
                        neutral_index = di + 1
                        if neutral_index < len(data):
                            p, stat, direction = kstest(data[di][1:], data[neutral_index][1:], 
                                           alternative='greater' if args.KStests==1 else "two-sided")
                            ksresults.append(", p={:.3g}, {}".format(p, direction))
                        else:
                            print(f"Error: Selected dataset at index {di} has no corresponding neutral control")
                            sys.exit(1)
                    else:
                        # This is a neutral dataset (odd indices: 1, 3, 5, ...)
                        ksresults.append("")  # No p-value for neutral datasets
            else:
                # Single neutral mode: find the specified neutral label
                if args.useheaderlabels:
                    # If using header labels, find neutral in headers which matches data order
                    try:
                        ni = headers.index(args.neutrallabel)
                    except ValueError:
                        print(f"Error: Neutral label '{args.neutrallabel}' not found in data headers")
                        print(f"Available labels: {headers}")
                        sys.exit(1)
                else:
                    # If using custom labels, find neutral in args.labels
                    try:
                        ni = args.labels.index(args.neutrallabel)
                    except ValueError:
                        print(f"Error: Neutral label '{args.neutrallabel}' not found in provided labels")
                        print(f"Available labels: {args.labels}")
                        sys.exit(1)
                
                ksresults = []
                for di,d in enumerate(data):
                    if di == ni:
                        ksresults.append("")
                    else:
                        p, stat, direction = kstest(d[1:],data[ni][1:], alternative='greater' if args.KStests==1 else "two-sided")
                        ksresults.append(", p={:.3g}, {}".format(p, direction))

            plot_data(data, counts, args.labels,args,ksresults=ksresults)    
        else:
            plot_data(data,counts, args.labels,args,ksresults=None)
