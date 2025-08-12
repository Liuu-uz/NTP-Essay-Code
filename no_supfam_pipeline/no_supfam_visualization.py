import os
import csv
import pymol
from pymol import cmd
import glob

def parse_extracted_csv(csv_file_path):
    """Parse extracted_first_rows.csv file"""
    results = []
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                # Extract query ID from Source_File (e.g.: Q67ZM7_zscores.csv -> Q67ZM7)
                query_id = row['Source_File'].replace('_zscores.csv', '')
                
                # Extract template ID from the last 5 characters of filename, take first 4
                filename = row['filename']
                # Remove .txt suffix
                filename_no_ext = filename.replace('.txt', '')
                
                if len(filename_no_ext) >= 5:
                    # Take last 5 characters, then take first 4
                    last_5_chars = filename_no_ext[-5:]
                    template_id = last_5_chars[:4]
                    
                    print(f"  Filename: {filename} -> Last 5 chars: {last_5_chars} -> Template ID: {template_id}")
                    
                    results.append({
                        'query_id': query_id,
                        'template_id': template_id,
                        'z_score': float(row['Z-score']),
                        'filename': filename
                    })
                else:
                    print(f"  Warning: Filename {filename} has less than 5 characters")
                    
        return results
    except Exception as e:
        print(f"Error parsing CSV file: {e}")
        return []

def find_query_pdb(predicted_dir, query_id):
    """Find query structure PDB file in predicted_structures directory"""
    # Try multiple possible filename formats
    possible_patterns = [
        f"{query_id.lower()}*.pdb",
        f"{query_id[:4].lower()}*.pdb",
        f"{query_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        query_pattern = os.path.join(predicted_dir, pattern)
        query_files = glob.glob(query_pattern)
        if query_files:
            return query_files[0]
    
    return None

def find_template_pdb(template_dir, template_id):
    """Find template PDB file in template directory"""
    if not os.path.exists(template_dir):
        print(f"    Template directory does not exist: {template_dir}")
        return None
    
    # Try multiple filename formats
    possible_patterns = [
        f"{template_id.lower()}*.pdb",
        f"{template_id.upper()}*.pdb", 
        f"{template_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        template_pattern = os.path.join(template_dir, pattern)
        template_files = glob.glob(template_pattern)
        if template_files:
            print(f"    Found template file: {template_files[0]} (using pattern: {pattern})")
            return template_files[0]
    
    print(f"    Template file not found, tried patterns: {possible_patterns}")
    return None

def setup_pymol_session():
    """Initialize PyMOL session"""
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()

def visualize_structures_with_atp(query_pdb, template_pdb, output_prefix, query_id, template_id, z_score):
    """Visualize structures using PyMOL, including ATP and active sites, without text labels"""
    
    # Clear current session
    cmd.delete('all')
    
    # Load structure files
    query_name = 'query_structure'
    template_name = 'template_structure'
    
    try:
        cmd.load(query_pdb, query_name)
        cmd.load(template_pdb, template_name)
    except Exception as e:
        print(f"    Failed to load PDB files: {e}")
        return None
    
    # Structure alignment - use CA atoms for alignment
    try:
        alignment_result = cmd.align(query_name, template_name)
        print(f"    Structure alignment result: RMSD = {alignment_result[0]:.3f} Å, Aligned atoms = {alignment_result[1]}")
    except Exception as e:
        print(f"    Structure alignment failed: {e}")
        return None
    
    # Set structure display style
    cmd.hide('everything')
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # Find and display ATP molecules
    atp_selection = 'resn ATP or resn ADP or resn AMP'
    nucleotide_found = False
    
    if cmd.count_atoms(atp_selection) > 0:
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        print("    Found ATP/ADP/AMP molecules")
        nucleotide_found = True
        
        # Show active site residues around ATP (within 5Å range)
        active_site = f'byres ({atp_selection} around 5)'
        cmd.show('sticks', active_site)
        cmd.color('yellow', f'{active_site} and not {atp_selection}')
        print("    Showing active site residues")
    else:
        # If no ATP, try to display other nucleotides
        other_nucleotides = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
        if cmd.count_atoms(other_nucleotides) > 0:
            cmd.show('sticks', other_nucleotides)
            cmd.color('red', other_nucleotides)
            active_site = f'byres ({other_nucleotides} around 5)'
            cmd.show('sticks', active_site)
            cmd.color('yellow', f'{active_site} and not {other_nucleotides}')
            print("    Found other nucleotide molecules")
            nucleotide_found = True
        else:
            print("    No nucleotide molecules found")
    
    # Set view
    cmd.zoom('all')
    cmd.orient()
    
    # Set high quality rendering parameters
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)  # Enable shadows for 3D effect
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    
    # Save white background version
    cmd.bg_color('white')
    png_file_white = f'{output_prefix}_clean_white.png'
    cmd.png(png_file_white, width=1800, height=1400, dpi=300, ray=1)
    
    # Save aligned structures
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    
    # Add title information to terminal output (not displayed in image)
    title_text = f"Query: {query_id} vs Template: {template_id} (Z-score: {z_score:.2f})"
    print(f"    {title_text}")
    
    return {
        'rmsd': alignment_result[0],
        'aligned_atoms': alignment_result[1],
        'nucleotide_found': nucleotide_found,
        'png_file_white': png_file_white,
        'pdb_file': pdb_file
    }

def process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """Process structure visualization based on CSV data"""
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize PyMOL
    setup_pymol_session()
    
    # Parse CSV file
    csv_data = parse_extracted_csv(csv_file_path)
    if not csv_data:
        print("Cannot parse CSV file or file is empty")
        return []
    
    print(f"Parsed {len(csv_data)} records from CSV file")
    
    # If maximum processing number is set, only process first N
    if max_structures:
        csv_data = csv_data[:max_structures]
        print(f"Limited to processing first {max_structures} structures")
    
    results = []
    successful_count = 0
    
    for i, data in enumerate(csv_data, 1):
        query_id = data['query_id']
        template_id = data['template_id']
        z_score = data['z_score']
        filename = data['filename']
        
        print(f"\n[{i}/{len(csv_data)}] Processing: {query_id} vs {template_id} (Z-score: {z_score:.2f})")
        print(f"  Original filename: {filename}")
        
        # Find query structure PDB file
        query_pdb = find_query_pdb(predicted_dir, query_id)
        if not query_pdb:
            print(f"    Query structure not found: {query_id}")
            continue
        
        # Find template structure PDB file
        template_pdb = find_template_pdb(template_dir, template_id)
        if not template_pdb:
            print(f"    Template structure not found: {template_id}")
            continue
        
        print(f"    Query structure: {query_pdb}")
        print(f"    Template structure: {template_pdb}")
        
        # Create output file prefix
        output_prefix = os.path.join(output_dir, f"{query_id}_vs_{template_id}_zscore{z_score:.1f}")
        
        try:
            # Perform structure visualization
            vis_result = visualize_structures_with_atp(
                query_pdb, template_pdb, output_prefix, 
                query_id, template_id, z_score
            )
            
            if vis_result:
                # Record results
                result = {
                    'query_id': query_id,
                    'template_id': template_id,
                    'z_score': z_score,
                    'filename': filename,
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb,
                    'rmsd': vis_result['rmsd'],
                    'aligned_atoms': vis_result['aligned_atoms'],
                    'nucleotide_found': vis_result['nucleotide_found'],
                    'png_file_white': vis_result['png_file_white'],
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    Success ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    Processing failed: {e}")
            continue
    
    # Save processing results summary
    summary_file = os.path.join(output_dir, 'visualization_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Query_ID', 'Template_ID', 'Z_Score', 'Original_Filename', 'Query_PDB', 'Template_PDB', 
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'PNG_Clean_White', 'PDB_File'
        ])
        
        for result in results:
            writer.writerow([
                result['query_id'],
                result['template_id'],
                f"{result['z_score']:.2f}",
                result['filename'],
                result['query_pdb'],
                result['template_pdb'],
                f"{result['rmsd']:.3f}",
                result['aligned_atoms'],
                result['nucleotide_found'],
                os.path.basename(result['png_file_white']),
                os.path.basename(result['pdb_file'])
            ])
    
    print(f"\nProcessing complete!")
    print(f"Total records: {len(csv_data)}")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed count: {len(csv_data) - successful_count}")
    print(f"Results summary saved to: {summary_file}")
    
    return results

def main():
    """Main function"""
    # Set paths
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    pipeline_dir = os.path.join(base_dir, "no_supfam_pipeline")
    
    csv_file_path = os.path.join(pipeline_dir, "extracted_first_rows.csv")
    predicted_dir = os.path.join(pipeline_dir, "predicted_structures")
    template_dir = os.path.join(pipeline_dir, "template")  # Explicitly specify template directory
    output_dir = os.path.join(pipeline_dir, "pymol_visualization_output")
    
    # Check if files and directories exist
    if not os.path.exists(csv_file_path):
        print(f"Error: CSV file does not exist: {csv_file_path}")
        return
    
    if not os.path.exists(predicted_dir):
        print(f"Error: Predicted structures directory does not exist: {predicted_dir}")
        return
    
    if not os.path.exists(template_dir):
        print(f"Error: Template directory does not exist: {template_dir}")
        return
    
    print("PyMOL structure visualization based on CSV file (using first 4 of last 5 characters from filename as template ID)")
    print(f"CSV file: {csv_file_path}")
    print(f"Predicted structures directory: {predicted_dir}")
    print(f"Template directory: {template_dir}")
    print(f"Output directory: {output_dir}")
    
    # Ask whether to limit processing count (for testing)
    response = input("\nLimit processing count? (Enter number to limit, press Enter to process all): ").strip()
    max_structures = None
    if response.isdigit():
        max_structures = int(response)
        print(f"Will only process first {max_structures} structures")
    
    # Process files
    results = process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures)
    
    # Display statistics
    if results:
        z_scores = [r['z_score'] for r in results]
        rmsds = [r['rmsd'] for r in results]
        nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
        
        print(f"\nStatistics:")
        print(f"Z-score range: {min(z_scores):.2f} - {max(z_scores):.2f}")
        print(f"Average Z-score: {sum(z_scores)/len(z_scores):.2f}")
        print(f"RMSD range: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"Average RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"Structures containing nucleotides: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
    
    # Quit PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()