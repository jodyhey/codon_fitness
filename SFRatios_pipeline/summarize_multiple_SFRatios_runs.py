#!/usr/bin/env python3
"""
Parameter Scraper Script
Extracts parameter values from SFRatios output files and creates a tab-delimited table.

intended for running after Run_multiple_SFRatios_jobs.py
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract parameter values from SFRatios output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scraper.py -i "/path/to/files/*.out" -o results.txt
  python scraper.py -i "results_*.txt" -o summary.tsv
        """
    )
    
    parser.add_argument('-i', dest='infilestring', required=True,
                       help='Path pattern to input files (can contain wildcards * and ?)')
    parser.add_argument('-o', dest='outfilename', required=True,
                       help='Output filename for tab-delimited table')
    
    return parser.parse_args()

def extract_parameters_from_file(filepath):
    """
    Extract parameters from the AIC block and SNP counts from data table in a single file.
    
    Args:
        filepath (str): Path to the file to process
        
    Returns:
        dict: Dictionary of parameter names to values, including SNP counts, or None if error
    """
    parameters = {}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Find the line starting with "AIC"
        aic_start_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith('AIC'):
                aic_start_idx = i
                break
        
        if aic_start_idx is None:
            print(f"Warning: No AIC block found in {filepath}", file=sys.stderr)
            return parameters
        
        # Process lines from AIC until first blank line to get parameters
        aic_end_idx = None
        for i in range(aic_start_idx, len(lines)):
            line = lines[i].strip()
            
            # Stop at first blank line
            if not line or line[0].isspace() or line[0] == '\t' or line[0] == '*':
                aic_end_idx = i
                break
            
            # Split on tabs or whitespace and look for parameter name and value
            parts = re.split(r'\s+', line, maxsplit=1)  # Split on first whitespace
            if len(parts) >= 2:
                param_name = parts[0].strip()
                value_part = parts[1].strip()
                
                # Extract the first floating point number from the value part
                # Handle formats like: "2212.549", "-4419.098", "0.0138	(0.0137 - 0.0139)"
                float_match = re.search(r'(-?\d+\.?\d*(?:[eE][+-]?\d+)?)', value_part)
                if float_match:
                    try:
                        value = float(float_match.group(1))
                        parameters[param_name] = value
                    except ValueError:
                        print(f"Warning: Could not parse value '{float_match.group(1)}' for parameter '{param_name}' in {filepath}", file=sys.stderr)
                else:
                    print(f"Warning: No numeric value found for parameter '{param_name}' in line: '{line}' in {filepath}", file=sys.stderr)
                # Capture confidence intervals in parentheses and normalize as "low,high"
                ci_match = re.search(r'\(([^()]*)\)', line)
                if ci_match:
                    ci_raw = ci_match.group(1).strip()
                    # Split only on a hyphen that has whitespace on both sides to avoid
                    # treating minus signs in negative numbers as separators.
                    parts = [p.strip() for p in re.split(r'\s+-\s+', ci_raw) if p.strip()]
                    if len(parts) == 2:
                        ci_text = ','.join(parts)
                    else:
                        # Fallback: remove all whitespace
                        ci_text = re.sub(r'\s+', '', ci_raw)
                    parameters[param_name + '_CI'] = ci_text
            else:
                print(f"Warning: Could not parse line: '{line}' in {filepath}", file=sys.stderr)
        
        # Now look for SNP counts in the data table after the dashes
        neut_snp_count = 0
        sel_snp_count = 0
        found_snp_counts = False
        
        # Start searching after the AIC block
        search_start = aic_end_idx if aic_end_idx is not None else aic_start_idx + 10
        
        for i in range(search_start, len(lines)):
            line = lines[i].strip()
            
            # Look for line of dashes
            if line.count('-') > 10:  # Line of many dashes
                # Check the next few lines for the data table header
                for j in range(1, 5):  # Check next 4 lines for header
                    if i + j < len(lines):
                        header_line = lines[i + j].strip()
                        if 'DataN' in header_line and 'DataS' in header_line:
                            # Found the header line, determine column indices
                            header_parts = header_line.split('\t')
                            # Remove empty parts and clean up
                            header_parts = [part.strip() for part in header_parts if part.strip()]
                            
                            try:
                                data_n_col_index = header_parts.index('DataN')
                                data_s_col_index = header_parts.index('DataS')
                                
                                # Now read the data rows starting from the line after the header
                                for k in range(i + j + 1, len(lines)):
                                    data_line = lines[k].strip()
                                    
                                    if not data_line:  # Empty line, end of table
                                        break
                                    
                                    # Skip lines that don't look like data (e.g., contain letters except in first column)
                                    if any(char.isalpha() for char in data_line[2:]):  # Allow letters in first 2 chars for "i" column
                                        continue
                                    
                                    data_parts = data_line.split('\t')
                                    # Remove empty parts and clean up
                                    data_parts = [part.strip() for part in data_parts if part.strip()]
                                    
                                    if len(data_parts) > max(data_n_col_index, data_s_col_index):
                                        try:
                                            # Extract DataN and DataS values
                                            data_n_value = float(data_parts[data_n_col_index])
                                            data_s_value = float(data_parts[data_s_col_index])
                                            
                                            neut_snp_count += data_n_value
                                            sel_snp_count += data_s_value
                                            
                                        except (ValueError, IndexError) as e:
                                            # If we can't parse this line, we might have reached the end of the table
                                            break
                                    else:
                                        # Line doesn't have enough columns, probably end of table
                                        break
                                
                                found_snp_counts = True
                                break
                                
                            except ValueError:
                                # DataN or DataS not found in header
                                continue
                
                if found_snp_counts:
                    break
        
        # Add SNP counts to parameters
        if found_snp_counts:
            parameters['NeutSNPcount'] = int(neut_snp_count) if neut_snp_count > 0 else 0
            parameters['SelSNPcount'] = int(sel_snp_count) if sel_snp_count > 0 else 0
        else:
            print(f"Warning: Could not find SNP count data table in {filepath}", file=sys.stderr)
            parameters['NeutSNPcount'] = 'na'
            parameters['SelSNPcount'] = 'na'
    
    except Exception as e:
        print(f"Error reading file {filepath}: {e}", file=sys.stderr)
        return {}
    
    return parameters

def process_all_files(file_pattern):
    """
    Process all files matching the pattern and extract parameters.
    
    Args:
        file_pattern (str): Glob pattern for input files
        
    Returns:
        tuple: (file_data, all_parameters) where:
            file_data (dict): filename -> parameter dict
            all_parameters (set): set of all unique parameter names
    """
    # Find all matching files
    matching_files = glob.glob(file_pattern)
    
    if not matching_files:
        print(f"Error: No files found matching pattern: {file_pattern}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(matching_files)} files matching pattern")
    
    file_data = {}
    all_parameters = set()
    
    for filepath in matching_files:
        filename = os.path.basename(filepath)
        print(f"Processing: {filename}")
        
        parameters = extract_parameters_from_file(filepath)
        
        if parameters:
            file_data[filename] = parameters
            all_parameters.update(parameters.keys())
            print(f"  Found {len(parameters)} parameters")
            
            # Show SNP counts if found
            if 'NeutSNPcount' in parameters and 'SelSNPcount' in parameters:
                print(f"    NeutSNPcount: {parameters['NeutSNPcount']}, SelSNPcount: {parameters['SelSNPcount']}")
        else:
            print(f"  No parameters found")
            file_data[filename] = {}
    
    print(f"\nTotal unique parameters found: {len(all_parameters)}")
    print(f"Parameters: {sorted(all_parameters)}")
    
    return file_data, all_parameters

def write_output_table(file_data, all_parameters, output_filename):
    """
    Write the results table to output file with SNP count columns after filename.
    
    Args:
        file_data (dict): filename -> parameter dict
        all_parameters (set): set of all unique parameter names
        output_filename (str): path to output file
    """
    # Create ordered list of parameters with SNP counts first, then others
    snp_count_params = ['NeutSNPcount', 'SelSNPcount']
    other_params = sorted([p for p in all_parameters if p not in snp_count_params])
    
    # Column order: filename, SNP counts, then other parameters
    ordered_parameters = snp_count_params + other_params
    
    try:
        with open(output_filename, 'w', encoding='utf-8') as f:
            # Write header
            header = ['filename'] + ordered_parameters
            f.write('\t'.join(header) + '\n')
            
            # Write data rows
            for filename in sorted(file_data.keys()):
                row = [filename]
                
                for param in ordered_parameters:
                    if param in file_data[filename]:
                        # Format floating point numbers consistently
                        value = file_data[filename][param]
                        if isinstance(value, float):
                            # Use scientific notation for very small/large numbers, regular for others
                            if abs(value) < 0.001 or abs(value) > 10000:
                                row.append(f"{value:.5e}")
                            else:
                                row.append(f"{value:.6g}")
                        else:
                            row.append(str(value))
                    else:
                        row.append('na')
                
                f.write('\t'.join(row) + '\n')
        
        print(f"\nResults written to: {output_filename}")
        
    except Exception as e:
        print(f"Error writing output file {output_filename}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main execution function"""
    args = parse_arguments()
    
    print(f"Input pattern: {args.infilestring}")
    print(f"Output file: {args.outfilename}")
    
    # Process all matching files
    file_data, all_parameters = process_all_files(args.infilestring)
    
    if not file_data:
        print("Error: No data extracted from any files", file=sys.stderr)
        sys.exit(1)
    
    # Write output table
    write_output_table(file_data, all_parameters, args.outfilename)
    
    print(f"\nProcessing complete!")
    print(f"Files processed: {len(file_data)}")
    print(f"Parameters extracted: {len(all_parameters)}")

if __name__ == '__main__':
    main()
