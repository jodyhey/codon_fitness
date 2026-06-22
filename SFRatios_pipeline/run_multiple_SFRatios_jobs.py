#!/usr/bin/env python3
"""
Multi-SFS Ratios Analysis Script
Runs multiple SFRatios.py jobs in parallel on pairs of SFS data
Command line interface matches SFRatios.py exactly, except:
1. Input file contains multiple SFS pairs (4 lines each)
2. Additional -j argument for parallel jobs
"""

import argparse
import os
import sys
import subprocess
import tempfile
import concurrent.futures
import re
import glob
from pathlib import Path
import time

# Add parent directory to Python path to find SFRatios.py and its modules
# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)
# sys.path.insert(0, parent_dir)

# Get the absolute directory of the script currently running
current_script_dir = os.path.dirname(os.path.abspath(__file__))

# Build the path to SFRatios.py
# Get the parent directory (Go up one level from the script)
# parent_dir = os.path.dirname(current_script_dir)
# SFRATIOS_PATH = os.path.join(parent_dir, 'SFRatios.py')
SFRatios_dir = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/SFRatios"
SFRATIOS_PATH = os.path.join(SFRatios_dir, 'SFRatios.py')
# print(SFRATIOS_PATH)
# exit()


def check_sfratios_availability():
    """Check if SFRatios.py is available in the expected directory"""
    if not os.path.exists(SFRATIOS_PATH):
        print(f"Error: SFRatios.py not found at {SFRATIOS_PATH}", file=sys.stderr)
        print(f"Expected location: {SFRatios_dir}", file=sys.stderr)
        sys.exit(1)
    try:
        result = subprocess.run(
            [sys.executable, SFRATIOS_PATH, '--help'],
            capture_output=True,
            text=True,
            cwd=SFRatios_dir
        )
        if result.returncode != 0:
            print(f"Warning: SFRatios.py found but may have issues: {result.stderr}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not test SFRatios.py: {e}", file=sys.stderr)

def parse_arguments():
    """Parse command line arguments - matches SFRatios.py exactly plus -j for jobs"""
    parser = argparse.ArgumentParser(
        description='Run multiple SFRatios.py jobs on pairs of SFS data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input file format:
Each SFS pair takes 4 lines:
  Line 1: Title for selected SFS
  Line 2: Selected SFS data
  Line 3: Title for neutral SFS  
  Line 4: Neutral SFS data
No blank lines between pairs.
        """
    )
    
    # Arguments matching SFRatios.py exactly
    parser.add_argument('-a', dest='inputfile', required=True,
                       help='Input file with multiple SFS pairs')
    parser.add_argument('-c', dest='fix_theta_ratio', type=float,
                       help='set the fixed value of thetaS/thetaN (i.e. mutation rate ratio)')
    parser.add_argument('-d', dest='densityof2Ns', 
                       help='gamma, lognormal, normal, fixed2Ns')
    parser.add_argument('-g', dest='globaloptimization', action='store_true',
                       help='turn on optimization using basinhopping and dualannealing (very slow, often finds better optimum)')
    parser.add_argument('-f', dest='foldstatus', required=True,
                       help="usage regarding folded or unfolded SFS distribution, 'isfolded', 'foldit' or 'unfolded'")
    parser.add_argument('-m', dest='setmax2Ns', type=float,
                       help='optional setting for 2Ns maximum, default = 0, use with -d lognormal or -d gamma')
    parser.add_argument('-M', dest='maxi', type=int,
                       help='the maximum bin index to include in the calculations, default=None')
    parser.add_argument('-p', dest='poplabel', required=True,
                       help='a population name and/or data type or other label for the start of the output filename')
    parser.add_argument('-t', dest='estimate_max2Ns', action='store_true',
                       help='if -d lognormal or -d gamma, estimate the maximum 2Ns value')
    parser.add_argument('-r', dest='resultsdir',
                       help='results directory')
    parser.add_argument('-x', dest='exitifexists', action='store_true',
                       help='if true and output file already exists, the run is stopped, else a new numbered output file is made')
    parser.add_argument('-z', dest='estimate_pointmass0', action='store_true',
                       help='include a proportion of the mass at zero in the density model')
    parser.add_argument('-Q', dest='thetaratio_range', nargs='+', type=float,
                       help='optional range for thetaratio (i.e. mutation rate ratio), low end followed by high end')
    parser.add_argument('-v', dest='mindenom', type=int,
                       help='optional setting the minimum count in a denominator, e.g. 10. If used, args.maxi is set to this bin number')
    parser.add_argument('-w', dest='sumbins', action='store_true',
                       help='requires -v, if true sums neutral SFS bins till mindenom is reached')
    parser.add_argument('-s', dest='sumintronSFSs', action='store_true',
                       help='use the sum of all intron SFSs for the intron SFS in all runs')
    
    # Additional argument for parallel processing
    parser.add_argument('-j', dest='parallel_jobs', type=int, default=1,
                       help='Number of parallel jobs to run')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode with detailed output and tracing')
    # Allow overriding the SFRatios directory (where SFRatios.py is located)
    parser.add_argument('-S', dest='sfratios_dir', default=None,
                       help=f'Path to SFRatios directory (default: {SFRatios_dir})')
    # New: indicate what the FIRST SFS represents in each 4-line block
    parser.add_argument('--first-sfs', dest='first_sfs', choices=['neutral', 'selected'], default='selected',
                        help="Meaning of the FIRST SFS in each 4-line block: 'selected' (default) or 'neutral'.")
    
    args = parser.parse_args()
    # normalize
    if hasattr(args, 'first_sfs') and isinstance(args.first_sfs, str):
        args.first_sfs = args.first_sfs.lower()
    return args

def parse_sfs_file(filename, sumintronSFS=False, first_sfs_default: str = 'selected'):
    """Parse the multi-SFS input file (robust to order by using labels).

    Accepts 4-line blocks: label1, sfs1, label2, sfs2. Detects which label
    corresponds to intron/neutral vs synonymous/selected by keyword match.
    """
    def is_neutral_label(lbl: str) -> bool:
        L = lbl.lower()
        return ('intron' in L) or ('neutral' in L)

    def is_selected_label(lbl: str) -> bool:
        L = lbl.lower()
        return ('synonymous' in L) or ('selected' in L)

    sfs_pairs = []
    try:
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        if len(lines) % 4 != 0:
            raise ValueError(f"Input file must have lines divisible by 4, got {len(lines)} lines")

        # If summing intron SFS is requested, accumulate over neutral slots found per block
        total_neutral = None
        if sumintronSFS:
            for i in range(0, len(lines), 4):
                l1, s1 = lines[i], lines[i + 1]
                l2, s2 = lines[i + 2], lines[i + 3]
                neutral_sfs_line = s1 if is_neutral_label(l1) else (s2 if is_neutral_label(l2) else None)
                if neutral_sfs_line is None:
                    # Fallback: assume second pair is neutral
                    neutral_sfs_line = s2
                arr = list(map(int, neutral_sfs_line.split()))
                if total_neutral is None:
                    total_neutral = arr[:]
                else:
                    total_neutral = [a + b for a, b in zip(total_neutral, arr)]
            total_neutral_str = ' '.join(map(str, total_neutral)) if total_neutral is not None else None

        for i in range(0, len(lines), 4):
            l1, s1 = lines[i], lines[i + 1]
            l2, s2 = lines[i + 2], lines[i + 3]
            # Determine roles by label
            if is_neutral_label(l1) and is_selected_label(l2):
                neutral_title, neutral_sfs = l1, s1
                selected_title, selected_sfs = l2, s2
            elif is_selected_label(l1) and is_neutral_label(l2):
                selected_title, selected_sfs = l1, s1
                neutral_title, neutral_sfs = l2, s2
            else:
                # Fallback: use user-declared ordering for FIRST SFS
                if (first_sfs_default or 'selected').lower() == 'neutral':
                    neutral_title, neutral_sfs = l1, s1
                    selected_title, selected_sfs = l2, s2
                else:
                    selected_title, selected_sfs = l1, s1
                    neutral_title, neutral_sfs = l2, s2

            if sumintronSFS and total_neutral_str is not None:
                neutral_title = "Introns"
                neutral_sfs = total_neutral_str

            sfs_pairs.append({
                'selected_title': selected_title,
                'selected_sfs': selected_sfs,
                'neutral_title': neutral_title,
                'neutral_sfs': neutral_sfs,
            })

    except Exception as e:
        print(f"Error parsing SFS file {filename}: {e}", file=sys.stderr)
        sys.exit(1)

    return sfs_pairs

def sanitize_label(label):
    """Convert label to safe filename component"""
    # Replace spaces with underscores and remove problematic characters
    sanitized = re.sub(r'[^\w\-_.]', '_', label.replace(' ', '_'))
    # Remove multiple consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    return sanitized.strip('_')

def create_temp_sfs_file(neutral_sfs, selected_sfs, neutral_title, selected_title):
    """Create temporary SFS input file for SFRatios.py.

    SFRatios.py expects line 2 = NEUTRAL SFS and line 4 = SELECTED SFS.
    Line 1 is a header, line 3 is ignored.
    """
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.sfs', delete=False)
    
    try:
        # Header: include both roles to avoid confusion in temp file
        header = f"Selected: {selected_title} | Neutral: {neutral_title}"
        temp_file.write(f"{header}\n")
        # NEUTRAL SFS on line 2
        temp_file.write(f"{neutral_sfs}\n")
        # Line 3 ignored by SFRatios.py; write a placeholder or blank
        temp_file.write("\n")
        # SELECTED SFS on line 4
        temp_file.write(f"{selected_sfs}\n")
        temp_file.flush()
        return temp_file.name
    finally:
        temp_file.close()

def run_sfratios_job(args, sfs_pair, job_index):
    """Run a single SFRatios.py job with identical arguments"""
    selected_label = sanitize_label(sfs_pair['selected_title'])
    job_poplabel = f"{args.poplabel}_{selected_label}"
    
    # Convert results directory to absolute path relative to where the batch script was called
    # This ensures output files go to the correct directory regardless of SFRatios.py's working directory
    if args.resultsdir:
        if os.path.isabs(args.resultsdir):
            absolute_resultsdir = args.resultsdir
        else:
            # Convert relative path to absolute path based on current working directory
            # (where the batch script was called from)
            absolute_resultsdir = os.path.abspath(args.resultsdir)
    else:
        absolute_resultsdir = None
    
    # Create temporary input file (ensure NEUTRAL first for SFRatios.py)
    if hasattr(args, 'debug') and args.debug:
        print(f"Mapping for SFRatios input: NEUTRAL first -> '{sfs_pair['neutral_title']}', SELECTED second -> '{sfs_pair['selected_title']}'")
    temp_input = create_temp_sfs_file(
        sfs_pair['neutral_sfs'],
        sfs_pair['selected_sfs'],
        sfs_pair['neutral_title'],
        sfs_pair['selected_title'],
    )
    
    try:
        # Build SFRatios.py command with correct path
        cmd = [
            sys.executable, SFRATIOS_PATH,  # Use full path to SFRatios.py
            '-a', temp_input,
            '-p', job_poplabel,
            '-f', args.foldstatus
        ]
        
        # Add debug/trace options if requested
        if hasattr(args, 'debug') and args.debug:
            # Insert trace option before the script
            cmd = [sys.executable, '-u', '-v', SFRATIOS_PATH] + cmd[2:]
        
        # Add all optional arguments exactly as they would be passed to SFRatios.py
        if args.fix_theta_ratio is not None:
            cmd.extend(['-c', str(args.fix_theta_ratio)])
        if args.densityof2Ns:
            cmd.extend(['-d', args.densityof2Ns])
        if args.globaloptimization:
            cmd.append('-g')
        if args.setmax2Ns is not None:
            cmd.extend(['-m', str(args.setmax2Ns)])
        if args.maxi is not None:
            cmd.extend(['-M', str(args.maxi)])
        if args.estimate_max2Ns:
            cmd.append('-t')
        if absolute_resultsdir:
            cmd.extend(['-r', absolute_resultsdir])
        if args.exitifexists:
            cmd.append('-x')
        if args.estimate_pointmass0:
            cmd.append('-z')
        if args.thetaratio_range:
            cmd.extend(['-Q'] + [str(x) for x in args.thetaratio_range])
        if args.mindenom is not None:
            cmd.extend(['-v', str(args.mindenom)])
        if args.sumbins:
            cmd.append('-w')
        
        # Run the command from the parent directory (where SFRatios.py and its modules are)
        print(f"Starting job {job_index + 1}: {selected_label}")
        if hasattr(args, 'debug') and args.debug:
            print(f"Working directory: {SFRatios_dir}")
            print(f"Temp input file: {temp_input}")
            print(f"Expected output pattern: {job_poplabel}*.out")
            print(f"Output directory (absolute): {absolute_resultsdir or 'current directory'}")
        print(f"Command: {' '.join(cmd)}")  # Debug: show the actual command
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            # timeout=3600,  # 1 hour timeout
            cwd=SFRatios_dir  # Run from parent directory where modules are available
        )
        
        # Check if job completed successfully
        success, error_msg, actual_output_file = check_job_success(result, job_poplabel, absolute_resultsdir)
        
        if hasattr(args, 'debug') and args.debug:
            # List files in output directory for debugging
            search_dir = absolute_resultsdir if absolute_resultsdir else '.'
            try:
                files_in_dir = os.listdir(search_dir)
                matching_files = [f for f in files_in_dir if job_poplabel in f and f.endswith('.out')]
                print(f"  Files in output dir matching pattern: {matching_files}")
            except Exception as e:
                print(f"  Could not list output directory: {e}")
        
        if not success:
            print(f"Job {job_index + 1} failed: {error_msg}", file=sys.stderr)
            print(f"  Return code: {result.returncode}", file=sys.stderr)
            if result.stderr:
                print(f"  stderr: {result.stderr[:300]}...", file=sys.stderr)
            if result.stdout:
                print(f"  stdout preview: {result.stdout[:300]}...", file=sys.stderr)
            # Return error information for inclusion in summary
            return {
                'job_index': job_index,
                'selected_label': selected_label,
                'poplabel': job_poplabel,
                'success': False,
                'error_msg': error_msg,
                'stdout_preview': result.stdout[:500] if result.stdout else "",
                'actual_output_file': actual_output_file,
                'temp_input': temp_input
            }
        
        print(f"Completed job {job_index + 1}: {selected_label}")
        print(f"  Output file: {actual_output_file}")
        
        # Return job results for successful job
        return {
            'job_index': job_index,
            'selected_label': selected_label,
            'poplabel': job_poplabel,
            'success': True,
            'error_msg': None,
            'stdout_preview': result.stdout[:500] if result.stdout else "",
            'actual_output_file': actual_output_file,
            'temp_input': temp_input
        }
    
    except subprocess.TimeoutExpired:
        error_msg = f"Job timed out after 1 hour"
        print(f"Job {job_index + 1} timed out after 1 hour", file=sys.stderr)
        return {
            'job_index': job_index,
            'selected_label': selected_label,
            'poplabel': job_poplabel,
            'success': False,
            'error_msg': error_msg,
            'stdout_preview': "",
            'actual_output_file': None,
            'temp_input': temp_input
        }
    except Exception as e:
        error_msg = f"Job exception: {str(e)}"
        print(f"Job {job_index + 1} error: {e}", file=sys.stderr)
        return {
            'job_index': job_index,
            'selected_label': selected_label,
            'poplabel': job_poplabel,
            'success': False,
            'error_msg': error_msg,
            'stdout_preview': "",
            'actual_output_file': None,
            'temp_input': temp_input
        }
    finally:
        # Clean up temporary file
        if os.path.exists(temp_input):
            os.unlink(temp_input)

def check_job_success(result, job_poplabel, resultsdir=None):
    """
    Check if SFRatios.py job completed successfully and capture any error information.
    
    Args:
        result: subprocess result object
        job_poplabel: the job's population label
        resultsdir: results directory (if specified)
        
    Returns:
        tuple: (success: bool, error_message: str or None, actual_output_file: str or None)
    """
    
    # Check for non-zero return code
    if result.returncode != 0:
        error_msg = f"Exit code {result.returncode}"
        if result.stderr:
            error_msg += f": {result.stderr.strip()[:200]}"
        return False, error_msg, None
    
    # Check for critical errors in stderr
    if result.stderr:
        stderr_text = result.stderr.strip()
        stderr_lower = stderr_text.lower()
        critical_errors = ['error', 'exception', 'traceback', 'failed', 'abort']
        for error in critical_errors:
            if error in stderr_lower:
                return False, f"Error in stderr: {stderr_text[:200]}", None
    
    # Check for error patterns in stdout
    if result.stdout:
        stdout_text = result.stdout
        stdout_lower = stdout_text.lower()
        error_patterns = [
            'nan', 'inf', 'convergence failed', 'optimization failed', 
            'no solution found', 'matrix is singular', 'error:', 'exception:'
        ]
        for pattern in error_patterns:
            if pattern in stdout_lower:
                # Find the line containing the error for better context
                for line in stdout_text.split('\n'):
                    if pattern in line.lower():
                        return False, f"Error detected: {line.strip()[:200]}", None
                return False, f"Error pattern detected: {pattern}", None
    
    # Try to find the actual output file created by SFRatios.py
    # SFRatios.py might create files with more complex naming schemes
    search_dir = resultsdir if resultsdir else '.'
    actual_output_file = None
    
    try:
        # Look for files that start with the job_poplabel and end with .out
        pattern = os.path.join(search_dir, f"{job_poplabel}*.out")
        matching_files = glob.glob(pattern)
        
        if matching_files:
            # Use the most recently created file
            actual_output_file = max(matching_files, key=os.path.getctime)
        else:
            # Try to parse the output filename from stdout
            if result.stdout and 'done' in result.stdout:
                lines = result.stdout.split('\n')
                for line in lines:
                    if 'done' in line and '.out' in line:
                        # Extract filename from line like "done  ./filename.out"
                        parts = line.strip().split()
                        for part in parts:
                            if part.endswith('.out') and job_poplabel in part:
                                actual_output_file = part
                                break
                        if actual_output_file:
                            break
        
        if not actual_output_file:
            return False, f"No output file found matching pattern: {pattern}", None
            
    except Exception as e:
        return False, f"Error searching for output file: {str(e)}", None
    
    # Check if the found output file has reasonable content
    try:
        with open(actual_output_file, 'r') as f:
            content = f.read()
        
        if len(content.strip()) < 50:  # Very short file suggests incomplete run
            return False, f"Output file too short ({len(content)} chars)", actual_output_file
        
        # Look for error indicators in the output file
        content_lower = content.lower()
        file_error_patterns = ['error', 'failed', 'exception', 'nan', 'inf']
        for pattern in file_error_patterns:
            if pattern in content_lower:
                # Find the specific line with the error
                for line in content.split('\n'):
                    if pattern in line.lower():
                        return False, f"Error in output file: {line.strip()[:200]}", actual_output_file
                return False, f"Error pattern in output file: {pattern}", actual_output_file
        
    except Exception as e:
        return False, f"Could not read output file {actual_output_file}: {str(e)}", actual_output_file
    
    # If we get here, the job appears to have completed successfully
    return True, None, actual_output_file

def parse_sfratios_output(output_file_path=None, output_text=None, debug_mode=False):
    """Parse the output from SFRatios.py to extract key results and SNP counts"""
    results = {}
    
    try:
        if output_file_path and os.path.exists(output_file_path):
            # Read from the actual output file
            with open(output_file_path, 'r') as f:
                content = f.read()
        elif output_text:
            # Use provided text (fallback)
            content = output_text
        else:
            return results
        
        lines = content.strip().split('\n')
        
        # First pass: extract basic results
        for line in lines:
            # Look for key result patterns - adapt based on actual SFRatios.py output format
            if 'Likelihood:' in line or 'likelihood' in line.lower():
                try:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        results['Likelihood'] = float(parts[1].strip())
                except ValueError:
                    pass
            elif 'AIC:' in line or 'aic' in line.lower():
                try:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        results['AIC'] = float(parts[1].strip())
                except ValueError:
                    pass
            elif 'Theta ratio:' in line or 'Thetaratio:' in line or 'theta ratio' in line.lower():
                try:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        results['Thetaratio'] = float(parts[1].strip())
                except ValueError:
                    pass
            elif line.startswith('Param'):
                # Handle parameter lines like "Param1: 1.234"
                try:
                    parts = line.split(':')
                    if len(parts) == 2:
                        param_name = parts[0].strip()
                        param_value = float(parts[1].strip())
                        results[param_name] = param_value
                except ValueError:
                    pass
        
        # Second pass: extract SNP counts from data table
        # Look for the data table after AIC block and dashes
        in_data_table = False
        data_n_col_index = None
        data_s_col_index = None
        neut_snp_count = 0
        sel_snp_count = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Look for line of dashes after AIC section
            if not in_data_table and line.count('-') > 10:  # Line of many dashes
                # Check if the next few lines contain a table header
                for j in range(1, 5):  # Check next 4 lines for header
                    if i + j < len(lines):
                        next_line = lines[i + j].strip()
                        if 'DataN' in next_line and 'DataS' in next_line:
                            # Found the header line
                            header_parts = next_line.split('\t')
                            try:
                                data_n_col_index = header_parts.index('DataN')
                                data_s_col_index = header_parts.index('DataS')
                                in_data_table = True
                                break
                            except ValueError:
                                # DataN or DataS not found in header
                                pass
            
            elif in_data_table and line:
                # We're in the data table, try to parse the line
                if line.startswith('i\t') or 'DataN' in line:
                    # Skip header lines
                    continue
                
                parts = line.split('\t')
                if len(parts) > max(data_n_col_index or 0, data_s_col_index or 0):
                    try:
                        # Extract DataN and DataS values
                        if data_n_col_index is not None and data_n_col_index < len(parts):
                            data_n_value = float(parts[data_n_col_index])
                            neut_snp_count += data_n_value
                        
                        if data_s_col_index is not None and data_s_col_index < len(parts):
                            data_s_value = float(parts[data_s_col_index])
                            sel_snp_count += data_s_value
                            
                    except (ValueError, IndexError):
                        # If we can't parse this line, we might have reached the end of the table
                        break
                else:
                    # Line doesn't have enough columns, probably end of table
                    break
        
        # Add SNP counts to results
        results['NeutSNPcount'] = int(neut_snp_count) if neut_snp_count > 0 else 'NA'
        results['SelSNPcount'] = int(sel_snp_count) if sel_snp_count > 0 else 'NA'
        
        # Debug output for SNP count extraction
        if debug_mode:
            print(f"Debug: Found DataN column at index {data_n_col_index}, DataS at index {data_s_col_index}")
            print(f"Debug: NeutSNPcount = {results['NeutSNPcount']}, SelSNPcount = {results['SelSNPcount']}")
    
    except Exception as e:
        print(f"Warning: Could not parse output: {e}", file=sys.stderr)
        # Set default values if parsing fails
        results['NeutSNPcount'] = 'NA'
        results['SelSNPcount'] = 'NA'
    
    return results

def main():
    """Main execution function"""
    args = parse_arguments()
    # Apply user override for SFRatios directory before checking availability
    if getattr(args, 'sfratios_dir', None):
        global SFRatios_dir, SFRATIOS_PATH
        SFRatios_dir = args.sfratios_dir
        SFRATIOS_PATH = os.path.join(SFRatios_dir, 'SFRatios.py')
    
    # Check if SFRatios.py is available
    print(f"Looking for SFRatios.py in parent directory: {SFRatios_dir}")
    check_sfratios_availability()
    print(f"Found SFRatios.py at: {SFRATIOS_PATH}")
    
    # Validate inputs
    if not os.path.exists(args.inputfile):
        print(f"Error: Input file {args.inputfile} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if args.parallel_jobs < 1:
        print("Error: Number of jobs must be at least 1", file=sys.stderr)
        sys.exit(1)
    
    # Parse SFS pairs
    print(f"Parsing SFS pairs from {args.inputfile}")
    sfs_pairs = parse_sfs_file(
        args.inputfile,
        sumintronSFS=getattr(args, 'sumintronSFS', False),
        first_sfs_default=getattr(args, 'first_sfs', 'selected'),
    )
    print(f"Found {len(sfs_pairs)} SFS pairs")
    
    if len(sfs_pairs) == 0:
        print("Error: No SFS pairs found in input file", file=sys.stderr)
        sys.exit(1)
    
    # Run jobs in parallel
    print(f"Running {len(sfs_pairs)} jobs with {args.parallel_jobs} parallel workers")
    start_time = time.time()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel_jobs) as executor:
        # Submit all jobs
        futures = []
        for i, sfs_pair in enumerate(sfs_pairs):
            future = executor.submit(run_sfratios_job, args, sfs_pair, i)
            futures.append(future)
        
        # Collect results
        job_results = []
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            job_results.append(result)
    
    end_time = time.time()
    print(f"All jobs completed in {end_time - start_time:.2f} seconds")
    
    # Sort results by job index
    job_results.sort(key=lambda x: x['job_index'] if x else float('inf'))
    
    # Count successful jobs (only those explicitly marked success=True)
    successful_jobs = sum(1 for result in job_results if result and result.get('success'))
    print(f"Successful jobs: {successful_jobs}/{len(sfs_pairs)}")
    
    if successful_jobs == 0:
        print("Error: No jobs completed successfully", file=sys.stderr)
        sys.exit(1)
    
    print("Multi-SFS analysis complete")

if __name__ == '__main__':
    main()
