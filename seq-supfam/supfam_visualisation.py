import os
import csv
import pymol
from pymol import cmd
import glob
import pandas as pd

def parse_xlsx_file(xlsx_file_path):
    """Parse xlsx file to extract EC number, sequenceID and PDBID+chain information"""
    results = []
    try:
        # Read xlsx file
        df = pd.read_excel(xlsx_file_path)
        
        print(f"Read {len(df)} rows of data")
        print(f"Column names: {list(df.columns)}")
        
        for index, row in df.iterrows():
            try:
                # Column 1: EC number
                ec_number = str(row.iloc[0]).strip()
                
                # Column 3: sequenceID  
                sequence_id = str(row.iloc[2]).strip()
                
                # Column 7: Matched PDBID+chain
                pdb_chain = str(row.iloc[6]).strip()
                
                # Validate data
                if ec_number and sequence_id and pdb_chain and pdb_chain != 'nan':
                    # Extract PDB ID (first 4 characters) and chain ID
                    if len(pdb_chain) >= 4:
                        pdb_id = pdb_chain[:4].lower()
                        chain_id = pdb_chain[4:] if len(pdb_chain) > 4 else 'A'
                        
                        results.append({
                            'ec_number': ec_number,
                            'sequence_id': sequence_id,
                            'pdb_chain': pdb_chain,
                            'pdb_id': pdb_id,
                            'chain_id': chain_id
                        })
                        
                        print(f"  Parsed: EC={ec_number}, SeqID={sequence_id}, PDB={pdb_id}, Chain={chain_id}")
                
            except Exception as e:
                print(f"  Skipping row {index+1}, parsing error: {e}")
                continue
                
        return results
        
    except Exception as e:
        print(f"Error parsing xlsx file: {e}")
        return []

def find_prediction_pdb(predicted_dir, sequence_id, ec_number):
    """Find predicted structure PDB file in supfam_pdb_files directory"""
    # File naming convention: sequenceID_ECnumber.pdb
    expected_filename = f"{sequence_id}_{ec_number}.pdb"
    expected_path = os.path.join(predicted_dir, expected_filename)
    
    if os.path.exists(expected_path):
        print(f"    Found prediction file: {expected_filename}")
        return expected_path
    
    # Try other possible formats
    possible_patterns = [
        f"{sequence_id}_{ec_number}*.pdb",
        f"{sequence_id}_*.pdb",
        f"*{sequence_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        search_pattern = os.path.join(predicted_dir, pattern)
        files = glob.glob(search_pattern)
        if files:
            print(f"    Found prediction file: {os.path.basename(files[0])} (using pattern: {pattern})")
            return files[0]
    
    print(f"    Prediction file not found, expected: {expected_filename}")
    return None

def find_template_pdb(template_dir, pdb_id, chain_id):
    """Find template PDB file in pdb_files directory"""
    if not os.path.exists(template_dir):
        print(f"    Template directory does not exist: {template_dir}")
        return None
    
    # Try multiple filename formats
    possible_patterns = [
        f"{pdb_id}.pdb",
        f"{pdb_id.upper()}.pdb",
        f"pdb{pdb_id}.ent",
        f"{pdb_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        template_pattern = os.path.join(template_dir, pattern)
        template_files = glob.glob(template_pattern)
        if template_files:
            print(f"    Found template file: {template_files[0]} (using pattern: {pattern})")
            return template_files[0]
    
    print(f"    Template file not found, tried PDB ID: {pdb_id}")
    return None

def setup_pymol_session():
    """Initialize PyMOL session"""
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()

def extract_chain_from_pdb(pdb_file, chain_id, output_path):
    """Extract specified chain from PDB file and save"""
    try:
        # Load PDB file
        temp_name = 'temp_structure'
        cmd.load(pdb_file, temp_name)
        
        # Check if specified chain exists
        chains = cmd.get_chains(temp_name)
        print(f"    Chains in template file: {chains}")
        
        if chain_id.upper() in [c.upper() for c in chains]:
            # Select specified chain
            chain_selection = f'{temp_name} and chain {chain_id.upper()}'
            cmd.create('extracted_chain', chain_selection)
            
            # Save extracted chain
            cmd.save(output_path, 'extracted_chain')
            
            # Cleanup
            cmd.delete(temp_name)
            cmd.delete('extracted_chain')
            
            print(f"    Successfully extracted chain {chain_id.upper()} to: {output_path}")
            return True
        else:
            print(f"    Warning: Chain {chain_id.upper()} does not exist, using complete structure")
            cmd.save(output_path, temp_name)
            cmd.delete(temp_name)
            return True
            
    except Exception as e:
        print(f"    Chain extraction failed: {e}")
        return False

def visualize_structures_with_atp(query_pdb, template_pdb, output_prefix, sequence_id, pdb_id, chain_id, ec_number):
    """Visualize structures using PyMOL, including ATP and active sites"""
    
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
    
    # Disable cartoon arrows/directional display
    cmd.set('cartoon_fancy_helices', 0)  # Disable fancy helix display
    cmd.set('cartoon_fancy_sheets', 0)   # Disable fancy sheet display
    cmd.set('cartoon_flat_sheets', 1)    # Use flat sheets
    cmd.set('cartoon_smooth_loops', 1)   # Smooth loop regions
    cmd.set('cartoon_highlight_color', 'grey50')  # Set highlight color
    cmd.set('cartoon_putty_transform', 0)  # Disable putty transformation
    cmd.set('cartoon_ring_mode', 0)      # Disable ring mode
    cmd.set('cartoon_tube_radius', 0.4)  # Set tube radius
    cmd.set('cartoon_oval_length', 1.2)  # Set oval length
    cmd.set('cartoon_oval_width', 0.3)   # Set oval width
    
    # Find and display ATP molecules
    atp_selection = 'resn ATP or resn ADP or resn AMP'
    nucleotide_found = False
    
    if cmd.count_atoms(atp_selection) > 0:
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        cmd.show('spheres', f'{atp_selection} and name P*')
        cmd.set('sphere_scale', 0.3, f'{atp_selection} and name P*')
        print("    Found ATP/ADP/AMP molecules")
        nucleotide_found = True
        
        # Show active site residues around ATP (within 4Å range)
        active_site = f'byres ({atp_selection} around 4)'
        cmd.show('sticks', f'{active_site} and sidechain')
        cmd.color('yellow', f'{active_site} and sidechain and not {atp_selection}')
        print("    Showing active site residues")
    else:
        # If no ATP, try to display other nucleotides
        other_nucleotides = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
        if cmd.count_atoms(other_nucleotides) > 0:
            cmd.show('sticks', other_nucleotides)
            cmd.color('red', other_nucleotides)
            cmd.show('spheres', f'{other_nucleotides} and name P*')
            cmd.set('sphere_scale', 0.3, f'{other_nucleotides} and name P*')
            active_site = f'byres ({other_nucleotides} around 4)'
            cmd.show('sticks', f'{active_site} and sidechain')
            cmd.color('yellow', f'{active_site} and sidechain and not {other_nucleotides}')
            print("    Found other nucleotide molecules")
            nucleotide_found = True
        else:
            print("    No nucleotide molecules found")
    
    # Set view
    cmd.zoom('all')
    cmd.orient()
    
    # Disable all settings that might produce arrows
    cmd.set('cgo_line_width', 1.0)
    cmd.set('dash_gap', 0.0)
    cmd.set('dash_length', 0.25)
    cmd.set('dash_round_ends', 1)
    cmd.hide('labels')  # Hide all labels
    cmd.hide('nonbonded')  # Hide non-bonded atoms
    
    # Set high quality rendering parameters, avoid arrow display
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    cmd.set('ray_opaque_background', 1)  # Ensure opaque background
    
    # Save image
    cmd.bg_color('white')
    png_file = f'{output_prefix}_visualization.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    
    # Save aligned structures
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    
    # Output information
    title_text = f"Query: {sequence_id} (EC: {ec_number}) vs Template: {pdb_id}_{chain_id}"
    print(f"    {title_text}")
    
    return {
        'rmsd': alignment_result[0],
        'aligned_atoms': alignment_result[1],
        'nucleotide_found': nucleotide_found,
        'png_file': png_file,
        'pdb_file': pdb_file
    }

def process_xlsx_data(xlsx_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """Process structure visualization based on xlsx data"""
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize PyMOL
    setup_pymol_session()
    
    # Parse xlsx file
    xlsx_data = parse_xlsx_file(xlsx_file_path)
    if not xlsx_data:
        print("Cannot parse xlsx file or file is empty")
        return []
    
    print(f"Parsed {len(xlsx_data)} records from xlsx file")
    
    # If maximum processing number is set, only process first N
    if max_structures:
        xlsx_data = xlsx_data[:max_structures]
        print(f"Limited to processing first {max_structures} structures")
    
    results = []
    successful_count = 0
    
    for i, data in enumerate(xlsx_data, 1):
        sequence_id = data['sequence_id']
        ec_number = data['ec_number']
        pdb_id = data['pdb_id']
        chain_id = data['chain_id']
        pdb_chain = data['pdb_chain']
        
        print(f"\n[{i}/{len(xlsx_data)}] Processing: {sequence_id} (EC: {ec_number}) vs {pdb_id}_{chain_id}")
        
        # Find predicted structure PDB file
        query_pdb = find_prediction_pdb(predicted_dir, sequence_id, ec_number)
        if not query_pdb:
            print(f"    Predicted structure not found: {sequence_id}_{ec_number}")
            continue
        
        # Find template structure PDB file
        template_pdb_raw = find_template_pdb(template_dir, pdb_id, chain_id)
        if not template_pdb_raw:
            print(f"    Template structure not found: {pdb_id}")
            continue
        
        # Extract specified chain from template PDB
        template_pdb = os.path.join(output_dir, f"temp_{pdb_id}_{chain_id}.pdb")
        if not extract_chain_from_pdb(template_pdb_raw, chain_id, template_pdb):
            print(f"    Template chain extraction failed: {pdb_id}_{chain_id}")
            continue
        
        print(f"    Predicted structure: {query_pdb}")
        print(f"    Template structure: {template_pdb}")
        
        # Create output file prefix
        output_prefix = os.path.join(output_dir, f"{sequence_id}_{ec_number}_vs_{pdb_id}_{chain_id}")
        
        try:
            # Perform structure visualization
            vis_result = visualize_structures_with_atp(
                query_pdb, template_pdb, output_prefix, 
                sequence_id, pdb_id, chain_id, ec_number
            )
            
            if vis_result:
                # Record results
                result = {
                    'sequence_id': sequence_id,
                    'ec_number': ec_number,
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'pdb_chain': pdb_chain,
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb,
                    'rmsd': vis_result['rmsd'],
                    'aligned_atoms': vis_result['aligned_atoms'],
                    'nucleotide_found': vis_result['nucleotide_found'],
                    'png_file': vis_result['png_file'],
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    Success ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    Processing failed: {e}")
            continue
        
        finally:
            # Clean up temporary files
            if os.path.exists(template_pdb):
                try:
                    os.remove(template_pdb)
                except:
                    pass
    
    # Save processing results summary
    summary_file = os.path.join(output_dir, 'visualization_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Sequence_ID', 'EC_Number', 'PDB_ID', 'Chain_ID', 'PDB_Chain', 'Query_PDB', 'Template_PDB',
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'PNG_File', 'PDB_File'
        ])
        
        for result in results:
            writer.writerow([
                result['sequence_id'],
                result['ec_number'],
                result['pdb_id'],
                result['chain_id'],
                result['pdb_chain'],
                result['query_pdb'],
                result['template_pdb'],
                f"{result['rmsd']:.3f}",
                result['aligned_atoms'],
                result['nucleotide_found'],
                os.path.basename(result['png_file']),
                os.path.basename(result['pdb_file'])
            ])
    
    print(f"\nProcessing complete!")
    print(f"Total records: {len(xlsx_data)}")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed count: {len(xlsx_data) - successful_count}")
    print(f"Results summary saved to: {summary_file}")
    
    return results

def main():
    """Main function"""
    # Set paths
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    
    xlsx_file_path = os.path.join(base_dir, "seq-supfam", "NTP_Analysis_Report.xlsx")
    predicted_dir = os.path.join(base_dir, "supfam_pdb_files")
    template_dir = os.path.join(base_dir, "pdb_files")
    output_dir = os.path.join(base_dir, "structure_visualization_output")
    
    # Check if files and directories exist
    if not os.path.exists(xlsx_file_path):
        print(f"Error: xlsx file does not exist: {xlsx_file_path}")
        return
    
    if not os.path.exists(predicted_dir):
        print(f"Error: Predicted structures directory does not exist: {predicted_dir}")
        return
    
    if not os.path.exists(template_dir):
        print(f"Error: Template directory does not exist: {template_dir}")
        return
    
    print("PyMOL structure visualization based on xlsx file")
    print(f"xlsx file: {xlsx_file_path}")
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
    results = process_xlsx_data(xlsx_file_path, predicted_dir, template_dir, output_dir, max_structures)
    
    # Display statistics
    if results:
        rmsds = [r['rmsd'] for r in results]
        nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
        
        print(f"\nStatistics:")
        print(f"RMSD range: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"Average RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"Structures containing nucleotides: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
    
    # Quit PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()