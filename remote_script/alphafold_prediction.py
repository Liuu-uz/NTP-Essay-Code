#!/usr/bin/env python3
"""
AlphaFold batch processing script - Sequence length filter: 200-1000 AA
Added MSA rate limiting handling and intelligent rest strategy
"""

import os
import json
import subprocess
import time
import random
import re
from pathlib import Path
import argparse

def check_system_environment():
    """Check system environment"""
    print("\nSystem Environment Check:")
    
    # Check Python version
    import sys
    print(f"   Python version: {sys.version.split()[0]}")
    
    # Check CPU information
    try:
        import multiprocessing
        cpu_count = multiprocessing.cpu_count()
        print(f"   CPU cores: {cpu_count}")
    except:
        print("   CPU info: Unable to retrieve")
    
    # Check memory information
    try:
        with open('/proc/meminfo', 'r') as f:
            meminfo = f.read()
            for line in meminfo.split('\n'):
                if 'MemTotal' in line:
                    mem_total = int(line.split()[1]) // 1024  # Convert to MB
                    print(f"   System memory: {mem_total//1024:.1f} GB")
                    break
                elif 'MemAvailable' in line:
                    mem_available = int(line.split()[1]) // 1024  # Convert to MB
                    print(f"   Available memory: {mem_available//1024:.1f} GB")
                    break
    except:
        print("   Memory info: Unable to retrieve")

def convert_cif_to_pdb(output_dir, sequence_name, global_pdb_dir="./all_pdb_files"):
    """Automatically convert CIF files to PDB files"""
    output_path = Path(output_dir)
    
    # Find CIF files
    sample_0_cif = None
    all_cif_files = list(output_path.glob("**/*.cif"))
    
    for cif_file in all_cif_files:
        if "sample_0" in cif_file.name or len(all_cif_files) == 1:
            sample_0_cif = cif_file
            break
    
    if not sample_0_cif:
        print(f"Warning: No CIF files found")
        return 0
    
    print(f"Converting CIF file found: {sample_0_cif.name}")
    
    # Find conversion script
    script_dir = Path(__file__).parent
    converter_script = script_dir / "convert_cif_to_pdb.py"
    
    if not converter_script.exists():
        print(f"Warning: Conversion script not found: {converter_script}")
        return 0
    
    # Create global PDB directory
    global_pdb_path = Path(global_pdb_dir)
    global_pdb_path.mkdir(exist_ok=True)
    
    try:
        # Generate PDB filename
        pdb_filename = f"{sequence_name}.pdb"
        pdb_file = global_pdb_path / pdb_filename
        
        # Call conversion script
        cmd = ['python3', str(converter_script), str(sample_0_cif), str(pdb_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print(f"   Success: {sample_0_cif.name} -> all_pdb_files/{pdb_filename}")
            
            # Verify PDB file content
            if pdb_file.exists() and pdb_file.stat().st_size > 0:
                with open(pdb_file, 'r') as f:
                    atom_lines = [line for line in f if line.startswith('ATOM')]
                    print(f"   Generated PDB contains {len(atom_lines)} atom records")
                return 1
            else:
                print(f"   Error: PDB file is empty")
                return 0
        else:
            print(f"   Conversion failed: {result.stderr.strip()}")
            return 0
            
    except Exception as e:
        print(f"   Conversion exception: {e}")
        return 0

def cleanup_intermediate_files(output_dir):
    """Clean up intermediate files"""
    output_path = Path(output_dir)
    cleaned_files = 0
    saved_space = 0
    
    keep_extensions = {'.pdb', '.json', '.cif'}
    keep_patterns = ['confidence', 'summary', 'result']
    protected_dirs = {'pdb_files'}
    
    try:
        for file_path in output_path.rglob('*'):
            if file_path.is_file():
                is_protected = any(protected_dir in file_path.parts for protected_dir in protected_dirs)
                
                if is_protected:
                    continue
                
                should_keep = (
                    file_path.suffix.lower() in keep_extensions or
                    any(pattern in file_path.name.lower() for pattern in keep_patterns)
                )
                
                if not should_keep:
                    try:
                        file_size = file_path.stat().st_size
                        file_path.unlink()
                        cleaned_files += 1
                        saved_space += file_size
                    except Exception as e:
                        print(f"Warning: Cannot delete {file_path}: {e}")
        
        if cleaned_files > 0:
            print(f"Cleanup complete: Deleted {cleaned_files} intermediate files, saved {saved_space/1024/1024:.1f} MB")
    
    except Exception as e:
        print(f"Warning: Error during cleanup: {e}")

def get_sequence_length(json_file):
    """Get sequence length"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list) and len(data) > 0:
            sequence = data[0]['sequences'][0]['proteinChain']['sequence']
        else:
            sequence = data['sequences'][0]['proteinChain']['sequence']
        
        return len(sequence)
    except:
        return 0

def parse_msa_status_from_output(output_text):
    """Parse MSA status from output"""
    msa_info = {
        'pending_detected': False,
        'sleep_time': 0,
        'msa_success': False,
        'error_detected': False
    }
    
    lines = output_text.split('\n')
    
    for line in lines:
        # Check PENDING status
        if 'Sleeping for' in line and 'PENDING' in line:
            match = re.search(r'Sleeping for (\d+)s', line)
            if match:
                msa_info['pending_detected'] = True
                msa_info['sleep_time'] = int(match.group(1))
        
        # Check MSA success
        elif 'update msa result success' in line:
            msa_info['msa_success'] = True
        
        # Check other errors
        elif any(error_keyword in line.lower() for error_keyword in ['error', 'failed', 'timeout']):
            msa_info['error_detected'] = True
    
    return msa_info

def calculate_smart_delay(consecutive_pendings, seq_length, is_first_request=False):
    """Calculate intelligent delay time"""
    
    if is_first_request:
        # No delay for first request
        return 0
    
    # Adjust base delay and multiplier based on consecutive PENDING count
    if consecutive_pendings == 0:
        # Normal case, base delay 60 seconds
        base_delay = 60
        pending_factor = 1.0
    elif consecutive_pendings == 1:
        # First PENDING, rest about 3 minutes
        base_delay = 180
        pending_factor = 1.0
    elif consecutive_pendings == 2:
        # Consecutive PENDING, rest about 5 minutes
        base_delay = 300
        pending_factor = 1.0
    else:
        # Severe rate limiting, rest about 8 minutes
        base_delay = 480
        pending_factor = 1.0
    
    # Fine-tune based on sequence length (longer sequences require more MSA time)
    length_factor = min(seq_length / 300, 1.5)
    
    # Randomize to avoid synchronized requests
    random_factor = random.uniform(0.8, 1.2)
    
    total_delay = base_delay * pending_factor * length_factor * random_factor
    
    # Limit delay range: minimum 15 seconds, maximum 600 seconds (10 minutes)
    return max(15, min(total_delay, 600))

def run_alphafold_prediction_optimized(json_file, output_base_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files", previous_msa_info=None):
    """Optimized AlphaFold prediction with MSA status handling"""
    sequence_name = json_file.stem
    output_dir = Path(output_base_dir) / f"{sequence_name}_output"
    
    # Get sequence length
    seq_length = get_sequence_length(json_file)
    
    print(f"\nPredicting: {json_file.name}")
    print(f"Sequence length: {seq_length} AA")
    
    # Build command - use MSA server
    cmd = [
        'protenix', 'predict',
        '--input', str(json_file),
        '--out_dir', str(output_dir),
        '--seeds', '101',
        '--use_msa_server'
    ]
    
    print(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        # Run prediction with real-time output monitoring
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        output_lines = []
        msa_start_time = None
        pending_detected = False
        
        # Read output in real-time
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                output_lines.append(output.strip())
                
                # Real-time detection of key status
                if 'starting to update msa result' in output:
                    msa_start_time = time.time()
                    print(f"   Starting MSA search...")
                
                elif 'Sleeping for' in output and 'PENDING' in output:
                    match = re.search(r'Sleeping for (\d+)s', output)
                    if match:
                        sleep_time = int(match.group(1))
                        print(f"   MSA server PENDING, waiting {sleep_time} seconds...")
                        pending_detected = True
                
                elif 'Files downloaded and extracted successfully' in output:
                    if msa_start_time:
                        msa_duration = time.time() - msa_start_time
                        print(f"   MSA download complete, took {msa_duration:.1f} seconds")
                
                elif 'update msa result success' in output:
                    if msa_start_time:
                        total_msa_time = time.time() - msa_start_time
                        print(f"   MSA processing complete, total time {total_msa_time:.1f} seconds")
        
        # Wait for process completion
        return_code = process.wait()
        duration = time.time() - start_time
        
        # Parse output
        full_output = '\n'.join(output_lines)
        msa_info = parse_msa_status_from_output(full_output)
        
        if return_code == 0:
            # Check output files
            pdb_files = list(output_dir.glob("**/*.pdb"))
            cif_files = list(output_dir.glob("**/*.cif"))
            structure_files = list(output_dir.glob("**/*structure*"))
            model_files = list(output_dir.glob("**/*model*"))
            all_files = list(output_dir.glob("**/*"))
            
            print(f"Output file check:")
            print(f"   PDB files: {len(pdb_files)}")
            print(f"   CIF files: {len(cif_files)}") 
            print(f"   Structure files: {len(structure_files)}")
            print(f"   Model files: {len(model_files)}")
            print(f"   Total files: {len([f for f in all_files if f.is_file()])}")
            
            # Convert CIF files
            converted_count = 0
            if cif_files:
                converted_count = convert_cif_to_pdb(output_dir, sequence_name, global_pdb_dir)
            else:
                print("   No CIF files found, skipping conversion")
            
            total_structure_files = len(pdb_files) + len(cif_files) + len(structure_files) + len(model_files)
            
            # Clean up intermediate files
            if cleanup_intermediate and total_structure_files > 0:
                print(f"Starting cleanup of intermediate files...")
                cleanup_intermediate_files(output_dir)
            
            print(f"Success! Duration: {duration/60:.1f} minutes, Structure files: {total_structure_files}")
            
            return {
                'status': 'success',
                'duration': duration,
                'pdb_count': len(pdb_files),
                'cif_count': len(cif_files),
                'structure_count': len(structure_files),
                'model_count': len(model_files),
                'converted_count': converted_count,
                'total_files': len([f for f in all_files if f.is_file()]),
                'sequence_length': seq_length,
                'cleaned_up': cleanup_intermediate and total_structure_files > 0,
                'msa_info': msa_info,
                'pending_detected': pending_detected
            }
        else:
            print(f"Failed! Duration: {duration/60:.1f} minutes")
            print(f"Return code: {return_code}")
            
            # Show partial output for debugging
            if full_output:
                print(f"Output summary:")
                print(full_output[-300:])  # Show last 300 characters
            
            return {
                'status': 'failed',
                'duration': duration,
                'error': full_output[-500:] if full_output else '',
                'sequence_length': seq_length,
                'msa_info': msa_info,
                'pending_detected': pending_detected
            }
    
    except Exception as e:
        print(f"Exception: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'sequence_length': seq_length
        }

def process_batch_optimized(input_dir, start_from=0, limit=None, output_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files"):
    """Optimized batch processing with intelligent rest strategy"""
    input_path = Path(input_dir)
    json_files = sorted(list(input_path.glob("*.json")))
    
    if not json_files:
        print(f"Error: No JSON files found in {input_dir}")
        return
    
    print(f"Found {len(json_files)} JSON files")
    
    # Filter and sort by sequence length (258-1000 AA)
    print("Sorting and filtering by sequence length (processing only 258-1000 AA)...")
    files_with_length = []
    too_short_count = 0
    too_long_count = 0
    
    for json_file in json_files:
        seq_length = get_sequence_length(json_file)
        if seq_length > 0:
            if 0 <= seq_length <= 1000:
                files_with_length.append((json_file, seq_length))
            elif seq_length < 0:
                too_short_count += 1
                if too_short_count <= 3:
                    print(f"Skipping short sequence: {json_file.name} ({seq_length} AA)")
            else:  # seq_length > 1000
                too_long_count += 1
                if too_long_count <= 3:
                    print(f"Skipping long sequence: {json_file.name} ({seq_length} AA)")
        else:
            print(f"Warning: Cannot get sequence length: {json_file.name}")
    
    if too_short_count > 3:
        print(f"... {too_short_count - 3} more short sequences skipped")
    if too_long_count > 3:
        print(f"... {too_long_count - 3} more long sequences skipped")
    
    print(f"Filter results:")
    print(f"  - Skipped short sequences: {too_short_count} (<258 AA)")
    print(f"  - Skipped long sequences: {too_long_count} (>1000 AA)")
    print(f"  - Processing queue: {len(files_with_length)} sequences (258-1000 AA)")
    
    if not files_with_length:
        print("Error: No sequences found meeting criteria")
        return
    
    # Sort by length (short to long)
    files_with_length.sort(key=lambda x: x[1])
    sorted_files = [item[0] for item in files_with_length]
    
    # Apply start position and limit
    if start_from > 0:
        sorted_files = sorted_files[start_from:]
        print(f"Starting from file {start_from + 1}")
    
    if limit:
        sorted_files = sorted_files[:limit]
        print(f"Limited to processing {limit} files")
    
    total_files = len(sorted_files)
    if total_files == 0:
        print("Error: No files to process")
        return
    
    # Check current time and give suggestions
    current_hour = time.localtime().tm_hour
    if current_hour in [9, 10, 14, 15, 20, 21, 22]:
        print(f"\nCurrent time slot ({current_hour}:00) may be peak hours for MSA server usage")
        print("   Recommend running during late night (2-6 AM) or midday (11 AM-1 PM) for best performance")
        proceed = input("   Continue? (y/n): ")
        if proceed.lower() != 'y':
            print("Batch processing cancelled")
            return
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"\nStarting MSA-optimized batch processing of {total_files} sequences...")
    print("=" * 80)
    
    results = []
    success_count = 0
    total_duration = 0
    start_time = time.time()
    consecutive_pendings = 0  # Consecutive PENDING counter
    previous_msa_info = None
    
    for i, json_file in enumerate(sorted_files, 1):
        current_length = get_sequence_length(json_file)
        
        print(f"\nProgress: [{i}/{total_files}] - {json_file.name} ({current_length} AA)")
        
        # Intelligent delay strategy
        if i > 1:  # No delay for first sequence
            delay_time = calculate_smart_delay(consecutive_pendings, current_length, is_first_request=False)
            
            if delay_time > 60:  # Show detailed info for delays over 1 minute
                print(f"   Intelligent rest: {delay_time/60:.1f} minutes (consecutive PENDING: {consecutive_pendings} times)")
            else:
                print(f"   Preventive delay: {delay_time:.0f} seconds")
            
            time.sleep(delay_time)
        
        # Execute prediction
        result = run_alphafold_prediction_optimized(
            json_file, output_dir, cleanup_intermediate, global_pdb_dir, previous_msa_info
        )
        result['file'] = json_file.name
        results.append(result)
        
        # Update statistics
        if result['status'] == 'success':
            success_count += 1
            total_duration += result['duration']
            avg_duration = total_duration / success_count
            
            remaining_files = total_files - i
            estimated_remaining_time = (remaining_files * avg_duration) / 3600
            
            print(f"Statistics: Success {success_count}/{i}, Average {avg_duration/60:.1f} min/file")
            print(f"Estimated remaining: {estimated_remaining_time:.1f} hours")
        
        # Update PENDING status tracking
        if result.get('pending_detected'):
            consecutive_pendings += 1
            print(f"   PENDING detected (consecutive #{consecutive_pendings})")
        else:
            consecutive_pendings = 0  # Reset counter
        
        # Force long rest when too many consecutive PENDING
        if consecutive_pendings >= 3:
            long_rest = 360 + random.uniform(60, 180)  # 6-9 minutes
            print(f"   Too many consecutive PENDING, forcing long rest {long_rest/60:.1f} minutes...")
            time.sleep(long_rest)
            consecutive_pendings = 0  # Reset counter
        
        # Save MSA info for next use
        if 'msa_info' in result:
            previous_msa_info = result['msa_info']
        
        print("-" * 50)
    
    # Final statistics
    batch_duration = time.time() - start_time
    failed_count = total_files - success_count
    
    print(f"\nMSA-optimized batch processing complete!")
    print("=" * 80)
    print(f"Summary:")
    print(f"  - Total files: {total_files}")
    print(f"  - Success: {success_count} ({success_count/total_files*100:.1f}%)")
    print(f"  - Failed: {failed_count} ({failed_count/total_files*100:.1f}%)")
    print(f"  - Total time: {batch_duration/3600:.1f} hours")
    if success_count > 0:
        print(f"  - Average prediction time: {total_duration/success_count/60:.1f} minutes")
        print(f"  - Actual average (including rest): {batch_duration/total_files/60:.1f} minutes/file")
    
    # Save results
    save_results_optimized(results, output_dir)

def save_results_optimized(results, output_dir):
    """Save optimized results report"""
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    report_file = Path(output_dir) / f"batch_results_optimized_{timestamp}.json"
    
    success_results = [r for r in results if r['status'] == 'success']
    pending_count = sum(1 for r in results if r.get('pending_detected', False))
    
    report = {
        'metadata': {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_files': len(results),
            'success_count': len(success_results),
            'failed_count': len(results) - len(success_results),
            'pending_encounters': pending_count,
            'length_filter': '258-1000 AA only',
            'optimization': 'MSA smart delay enabled'
        },
        'statistics': {
            'avg_duration': sum(r.get('duration', 0) for r in success_results) / len(success_results) if success_results else 0,
            'total_prediction_time': sum(r.get('duration', 0) for r in success_results),
            'pending_rate': pending_count / len(results) if results else 0,
            'sequence_length_stats': {
                'min': min(r.get('sequence_length', 0) for r in results if r.get('sequence_length', 0) > 0) if results else 0,
                'max': max(r.get('sequence_length', 0) for r in results),
                'avg': sum(r.get('sequence_length', 0) for r in results) / len(results) if results else 0
            }
        },
        'results': results
    }
    
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    
    print(f"Detailed report saved: {report_file}")
    print(f"PENDING encounter rate: {pending_count}/{len(results)} ({pending_count/len(results)*100:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description='AlphaFold batch processing script (MSA optimized)')
    parser.add_argument('--input_dir', '-i', default='/home/webserver/student/students_webserver/zhijing/input_jsons', help='JSON file input directory')
    parser.add_argument('--output_dir', '-o', default='./alphafold_predictions', help='Output directory')
    parser.add_argument('--start_from', '-f', type=int, default=0, help='Start from which file number')
    parser.add_argument('--limit', '-l', type=int, help='Limit number of files to process')
    parser.add_argument('--no_cleanup', action='store_true', help='Do not clean up intermediate files')
    
    args = parser.parse_args()
    
    print("AlphaFold Batch Processor (MSA Optimized)")
    print("Sequence length filter: 258-1000 amino acids")
    print("Intelligent features: MSA rate limiting detection + adaptive rest strategy")
    print("=" * 60)
    
    # Check system environment
    check_system_environment()
    
    print("Preparing to start batch processing...")
    
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory does not exist: {args.input_dir}")
        return
    
    print(f"\nInput directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Processing range: Starting from #{args.start_from + 1}" + (f", limited to {args.limit} files" if args.limit else ", no limit"))
    print(f"Auto cleanup: {'Disabled' if args.no_cleanup else 'Enabled'}")
    print("=" * 60)
    
    # Start optimized batch processing
    process_batch_optimized(
        input_dir=args.input_dir,
        start_from=args.start_from,
        limit=args.limit,
        output_dir=args.output_dir,
        cleanup_intermediate=not args.no_cleanup,
        global_pdb_dir="./all_pdb_files"
    )

if __name__ == "__main__":
    main()