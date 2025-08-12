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

def set_high_quality_render():
    """Set high quality rendering parameters"""
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    cmd.set('antialias', 2)
    cmd.set('ray_opaque_background', 1)

def create_view1_protein_chains(query_name, template_name, output_prefix):
    """View 1: Comparison of two protein chains"""
    print("    Generating View 1: Protein chain comparison")
    
    # Hide all displays
    cmd.hide('everything')
    
    # Show protein backbone
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # Set view
    cmd.zoom('all')
    cmd.orient()
    
    # Set rendering quality
    set_high_quality_render()
    cmd.bg_color('white')
    
    # Save image
    png_file = f'{output_prefix}_view1_chains.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    Saved: {png_file}")
    
    return png_file

def create_view2_with_atp(query_name, template_name, output_prefix):
    """View 2: Show only overlapping regions of two chains + ATP molecules in overlapping region"""
    print("    Generating View 2: Overlapping region protein chains + ATP")
    
    # Hide all displays
    cmd.hide('everything')
    
    # Find overlapping region of two chains (using 8Å distance to define overlap)
    overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
    overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
    overlap_region = f'({overlap_query}) or ({overlap_template})'
    
    print(f"    Query chain overlapping residues: {cmd.count_atoms(f'{overlap_query} and name CA')}")
    print(f"    Template chain overlapping residues: {cmd.count_atoms(f'{overlap_template} and name CA')}")
    
    # Show overlapping region protein chains
    if cmd.count_atoms(overlap_region) > 0:
        cmd.show('cartoon', overlap_region)
        cmd.color('cyan', overlap_query)
        cmd.color('orange', overlap_template)
        print("    Showing overlapping region protein chains")
    else:
        print("    Warning: No overlapping region found, showing complete protein chains")
        # If no overlapping region, show complete chains
        cmd.show('cartoon', 'all')
        cmd.color('cyan', query_name)
        cmd.color('orange', template_name)
        overlap_region = 'all'
    
    # Find all ATP molecules
    all_atp_selection = 'resn ATP or resn ADP or resn AMP'
    if cmd.count_atoms(all_atp_selection) == 0:
        all_atp_selection = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
    
    nucleotide_found = False
    
    if cmd.count_atoms(all_atp_selection) > 0:
        # Only show ATP in overlapping region (ATP within 6Å of overlapping region)
        overlap_atp = f'{all_atp_selection} and (({overlap_region}) around 6.0)'
        
        if cmd.count_atoms(overlap_atp) > 0:
            # Further filter: select 1-2 ATP molecules closest to overlapping region
            cmd.select('temp_overlap_atp', overlap_atp)
            atp_list = []
            cmd.iterate('temp_overlap_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
            
            if atp_list:
                # Limit ATP display, keep only first 2 (if multiple exist)
                selected_atp_residues = atp_list[:2]  # Maximum 2 ATP molecules
                
                # Build precise ATP selection
                atp_parts = []
                for chain, resi in selected_atp_residues:
                    atp_parts.append(f'(chain {chain} and resi {resi})')
                
                if atp_parts:
                    final_atp = f'({" or ".join(atp_parts)}) and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    
                    cmd.show('sticks', final_atp)
                    cmd.color('red', final_atp)
                    cmd.show('spheres', f'{final_atp} and name P*')
                    cmd.set('sphere_scale', 0.3, f'{final_atp} and name P*')
                    nucleotide_found = True
                    atp_selection = final_atp
                    
                    print(f"    Showing filtered ATP: {len(selected_atp_residues)} molecules, {cmd.count_atoms(final_atp)} atoms")
                    
                    # Statistics
                    total_atp_molecules = len(atp_list)
                    displayed_molecules = len(selected_atp_residues)
                    hidden_molecules = total_atp_molecules - displayed_molecules
                    print(f"    ATP molecule statistics: Total in overlap region {total_atp_molecules}, Displayed {displayed_molecules}, Hidden {hidden_molecules}")
            
            cmd.delete('temp_overlap_atp')
        else:
            print("    No ATP molecules found within 6Å of overlapping region")
            # If no ATP within 6Å, try larger range but limit quantity
            extended_atp = f'{all_atp_selection} and (({overlap_region}) around 12.0)'
            if cmd.count_atoms(extended_atp) > 0:
                cmd.select('temp_extended_atp', extended_atp)
                atp_list = []
                cmd.iterate('temp_extended_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
                
                if atp_list:
                    # Show only closest 1 ATP
                    selected_atp = atp_list[:1]
                    chain, resi = selected_atp[0]
                    final_atp = f'chain {chain} and resi {resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    
                    cmd.show('sticks', final_atp)
                    cmd.color('red', final_atp)
                    cmd.show('spheres', f'{final_atp} and name P*')
                    cmd.set('sphere_scale', 0.3, f'{final_atp} and name P*')
                    nucleotide_found = True
                    atp_selection = final_atp
                    print(f"    Showing extended range closest ATP: 1 molecule")
                
                cmd.delete('temp_extended_atp')
            else:
                atp_selection = ""
    else:
        print("    Warning: No nucleotide molecules found")
        atp_selection = ""
    
    # Set view - focus on overlapping region
    if nucleotide_found:
        focus_selection = f'{overlap_region} or {atp_selection}'
    else:
        focus_selection = overlap_region
    
    if cmd.count_atoms(focus_selection) > 0:
        cmd.zoom(focus_selection)
        cmd.orient(focus_selection)
    else:
        cmd.zoom('all')
        cmd.orient()
    
    # Set rendering quality
    set_high_quality_render()
    cmd.bg_color('white')
    
    # Save image
    png_file = f'{output_prefix}_view2_overlap_with_atp.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    Saved: {png_file}")
    
    return png_file, nucleotide_found

def create_view3_atp_closeup(query_name, template_name, output_prefix):
    """View 3: Direct magnification of ATP region based on View 2 - maintain same colors and style"""
    print("    Generating View 3: ATP region magnification based on View 2")
    
    # Don't clear display! Keep all content and colors from View 2
    # cmd.hide('everything') - commented out
    
    # Find ATP (using same logic as View 2)
    all_atp_selection = 'resn ATP or resn ADP or resn AMP'
    if cmd.count_atoms(all_atp_selection) == 0:
        all_atp_selection = 'resn GTP or resn GDP or resn GMP+CTP+CDP+CMP+UTP+UDP+UMP'
    
    # Find overlapping region (exactly same as View 2)
    overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
    overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
    overlap_region = f'({overlap_query}) or ({overlap_template})'
    
    # Find ATP in overlapping region, but select only one
    overlap_atp = f'{all_atp_selection} and (({overlap_region}) around 6.0)'
    target_atp = ""
    
    if cmd.count_atoms(overlap_atp) > 0:
        cmd.select('temp_overlap_atp', overlap_atp)
        atp_list = []
        cmd.iterate('temp_overlap_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
        
        if atp_list:
            # Select only first ATP molecule
            chain, resi = atp_list[0]
            target_atp = f'chain {chain} and resi {resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
            print(f"    Selected single ATP molecule: chain {chain} residue {resi}")
            
            # If multiple ATPs exist, hide others
            if len(atp_list) > 1:
                for other_chain, other_resi in atp_list[1:]:
                    other_atp = f'chain {other_chain} and resi {other_resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    cmd.hide('everything', other_atp)
                print(f"    Hidden other {len(atp_list)-1} ATP molecules")
        
        cmd.delete('temp_overlap_atp')
    
    if not target_atp:
        print("    Cannot find target ATP, unable to generate View 3")
        return None
    
    # Define amino acids within 4Å of ATP (add sidechain display, but don't change main colors)
    active_site_4a = f'byres ({target_atp} around 4.0) and polymer.protein'
    
    # Add sidechains of amino acids within 4Å range based on existing display
    if cmd.count_atoms(active_site_4a) > 0:
        # Show sidechains within 4Å range
        cmd.show('sticks', f'{active_site_4a} and sidechain')
        
        # Sidechain colors slightly deeper but coordinated with backbone
        query_sidechains = f'{active_site_4a} and sidechain and {query_name}'
        template_sidechains = f'{active_site_4a} and sidechain and {template_name}'
        
        # Use colors coordinated with backbone
        cmd.color('forest', query_sidechains)      # Deep cyan (coordinated with cyan backbone)
        cmd.color('gray50', template_sidechains)   # Deep orange (coordinated with orange backbone)
        
        # Set sidechain stick size
        cmd.set('stick_radius', 0.15, f'{active_site_4a} and sidechain')
        
        print(f"    Added 4Å amino acid sidechains: {cmd.count_atoms(f'{active_site_4a} and name CA')} residues")
    
    # Add hydrogen bond display (maintain continuity with View 2)
    try:
        cmd.distance('atp_closeup_hbonds', f'{target_atp}', f'{active_site_4a}', mode=2, cutoff=3.2)
        cmd.hide('labels', 'atp_closeup_hbonds')
        cmd.color('yellow', 'atp_closeup_hbonds')
        cmd.set('dash_width', 2, 'atp_closeup_hbonds')
        print("    Added ATP interaction hydrogen bonds")
    except:
        print("    Hydrogen bond display failed")
    
    # Key: Super magnification, let ATP+amino acids occupy 1/2 of the screen
    focus_selection = f'{target_atp} or {active_site_4a}'
    
    try:
        # Super extreme focus - ATP+amino acids occupy 1/2 of screen
        cmd.zoom(target_atp, buffer=0.1)  # Extreme magnification of ATP itself
        
        # Then include amino acids, but still very large overall (occupy 1/2 of screen)
        cmd.zoom(focus_selection, buffer=0.3)  # Very small buffer, ATP+amino acids occupy 1/2 of screen
        
        # Adjust viewing angle to ensure ATP and amino acids are clearly visible
        cmd.turn('x', 10)   # Tilt down to see ATP structure
        cmd.turn('y', 20)   # Side angle to show interactions
        cmd.turn('z', -10)  # Fine-tune rotation for optimal orientation
        
        print("    Super magnification: ATP+amino acids occupy 1/2 of screen size")
        
    except Exception as view_error:
        print(f"    View focusing failed: {view_error}")
        # Most extreme backup plan
        cmd.zoom(target_atp, buffer=0.05)  # Maximum magnification
        cmd.zoom(focus_selection, buffer=0.2)
    
    # Further optimize display, highlight ATP and key amino acids
    try:
        # Make ATP more prominent
        cmd.set('stick_radius', 0.3, target_atp)  # Thicker ATP sticks
        cmd.set('sphere_scale', 0.6, f'{target_atp} and name P*')  # Larger phosphate groups
        
        # Amino acid sidechains should also be clear
        cmd.set('stick_radius', 0.25, f'{active_site_4a} and sidechain')  # Thicker sidechain sticks
        
        # Hide protein chains that might block ATP view
        # Define region in front of ATP (judged by distance and position)
        foreground_region = f'byres ({target_atp} around 8.0) and polymer.protein'
        
        # Selectively hide parts of chains that might block ATP
        # Method 1: Hide chain segments above and in front of ATP
        blocking_chains = f'{foreground_region} and not ({active_site_4a})'
        
        if cmd.count_atoms(blocking_chains) > 0:
            # First try making these chain segments transparent
            cmd.set('cartoon_transparency', 0.8, blocking_chains)
            print(f"    Set foreground chain segments to high transparency")
            
            # If still blocking, can completely hide some chain segments
            # Select parts far from ATP but still in view
            distant_blocking = f'byres ({target_atp} around 10.0) and polymer.protein and not (byres ({target_atp} around 6.0))'
            if cmd.count_atoms(distant_blocking) > 0:
                cmd.hide('cartoon', distant_blocking)
                print(f"    Hidden some chain segments blocking ATP")
        
        print("    Enhanced ATP and amino acid visualization, reduced blocking")
    except Exception as e:
        print(f"    Error during display optimization: {e}")
        pass
    
    # Set rendering quality
    set_high_quality_render()
    cmd.bg_color('white')
    
    # Save giant magnification image - ATP dominates the screen
    png_file = f'{output_prefix}_view3_atp_closeup.png'
    try:
        # Ultra-high resolution display of giant ATP details
        cmd.png(png_file, width=4000, height=3000, dpi=300, ray=1)  # 4K resolution
        print(f"    Saved ATP giant close-up image (occupying 1/2 of screen): {png_file}")
    except Exception as save_error:
        print(f"    Image save failed: {save_error}")
        return None
    
    # Statistics
    try:
        active_site_count = cmd.count_atoms(active_site_4a)
        atp_count = cmd.count_atoms(target_atp)
        query_active_residues = cmd.count_atoms(f'{active_site_4a} and {query_name} and name CA')
        template_active_residues = cmd.count_atoms(f'{active_site_4a} and {template_name} and name CA')
        
        print(f"    ATP magnification statistics:")
        print(f"      ATP atoms: {atp_count}")
        print(f"      Query chain active site residues (4Å): {query_active_residues}")
        print(f"      Template chain active site residues (4Å): {template_active_residues}")
        
        return {
            'png_file': png_file,
            'active_site_atoms': active_site_count,
            'atp_atoms': atp_count,
            'query_active_residues': query_active_residues,
            'template_active_residues': template_active_residues,
            'overlap_atp_atoms': atp_count
        }
    except Exception as stat_error:
        print(f"    Statistics calculation failed: {stat_error}")
        return {
            'png_file': png_file,
            'active_site_atoms': 0,
            'atp_atoms': 0,
            'query_active_residues': 0,
            'template_active_residues': 0,
            'overlap_atp_atoms': 0
        }

def visualize_structures_three_views(query_pdb, template_pdb, output_prefix, query_id, template_id, z_score):
    """Generate three-view structure visualization using PyMOL"""
    
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
    
    # Generate three views
    results = {}
    
    # View 1: Protein chain comparison
    view1_file = create_view1_protein_chains(query_name, template_name, output_prefix)
    results['view1_chains'] = view1_file
    
    # View 2: Protein chains + ATP
    view2_file, nucleotide_found = create_view2_with_atp(query_name, template_name, output_prefix)
    results['view2_with_atp'] = view2_file
    results['nucleotide_found'] = nucleotide_found
    
    # View 3: ATP active site magnification
    if nucleotide_found:
        view3_result = create_view3_atp_closeup(query_name, template_name, output_prefix)
        if view3_result:
            results['view3_atp_closeup'] = view3_result['png_file']
            results['active_site_atoms'] = view3_result['active_site_atoms']
            results['atp_atoms'] = view3_result['atp_atoms']
            results['query_active_residues'] = view3_result['query_active_residues']
            results['template_active_residues'] = view3_result['template_active_residues']
            results['overlap_atp_atoms'] = view3_result['overlap_atp_atoms']
        else:
            results['view3_atp_closeup'] = None
    else:
        results['view3_atp_closeup'] = None
        print("    Skipping View 3: No nucleotide molecules found")
    
    # Save aligned structures
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    results['pdb_file'] = pdb_file
    
    # Add alignment information
    results['rmsd'] = alignment_result[0]
    results['aligned_atoms'] = alignment_result[1]
    
    # Output summary information
    title_text = f"Query: {query_id} vs Template: {template_id} (Z-score: {z_score:.2f})"
    print(f"    {title_text}")
    print(f"    Generated files:")
    print(f"      View 1 (protein chains): {os.path.basename(results['view1_chains'])}")
    print(f"      View 2 (with ATP): {os.path.basename(results['view2_with_atp'])}")
    if results['view3_atp_closeup']:
        print(f"      View 3 (ATP magnification): {os.path.basename(results['view3_atp_closeup'])}")
        print(f"      Active site atoms: {results.get('active_site_atoms', 'N/A')}")
    
    return results

def process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """Process structure visualization based on CSV data - three views version"""
    
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
            # Perform three-view structure visualization
            vis_result = visualize_structures_three_views(
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
                    'view1_chains': vis_result['view1_chains'],
                    'view2_with_atp': vis_result['view2_with_atp'],
                    'view3_atp_closeup': vis_result['view3_atp_closeup'],
                    'active_site_atoms': vis_result.get('active_site_atoms', None),
                    'atp_atoms': vis_result.get('atp_atoms', None),
                    'query_active_residues': vis_result.get('query_active_residues', None),
                    'template_active_residues': vis_result.get('template_active_residues', None),
                    'overlap_atp_atoms': vis_result.get('overlap_atp_atoms', None),
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    Success ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    Processing failed: {e}")
            continue
    
    # Save processing results summary
    summary_file = os.path.join(output_dir, 'visualization_summary_three_views.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Query_ID', 'Template_ID', 'Z_Score', 'Original_Filename', 'Query_PDB', 'Template_PDB', 
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'View1_Chains', 'View2_With_ATP', 
            'View3_ATP_Closeup', 'Active_Site_Atoms', 'ATP_Atoms', 'Query_Active_Residues', 
            'Template_Active_Residues', 'Overlap_ATP_Atoms', 'PDB_File'
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
                os.path.basename(result['view1_chains']) if result['view1_chains'] else '',
                os.path.basename(result['view2_with_atp']) if result['view2_with_atp'] else '',
                os.path.basename(result['view3_atp_closeup']) if result['view3_atp_closeup'] else '',
                result['active_site_atoms'] if result['active_site_atoms'] else '',
                result['atp_atoms'] if result['atp_atoms'] else '',
                result['query_active_residues'] if result['query_active_residues'] else '',
                result['template_active_residues'] if result['template_active_residues'] else '',
                result['overlap_atp_atoms'] if result['overlap_atp_atoms'] else '',
                os.path.basename(result['pdb_file']) if result['pdb_file'] else ''
            ])
    
    print(f"\nThree-view visualization processing complete!")
    print(f"Total records: {len(csv_data)}")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed count: {len(csv_data) - successful_count}")
    print(f"Results summary saved to: {summary_file}")
    
    # Statistics on view generation
    view1_count = sum(1 for r in results if r['view1_chains'])
    view2_count = sum(1 for r in results if r['view2_with_atp'])
    view3_count = sum(1 for r in results if r['view3_atp_closeup'])
    nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
    
    print(f"\nView generation statistics:")
    print(f"View 1 (protein chains): {view1_count}/{successful_count}")
    print(f"View 2 (with ATP): {view2_count}/{successful_count}")
    print(f"View 3 (ATP magnification): {view3_count}/{successful_count}")
    print(f"Structures containing nucleotides: {nucleotide_count}/{successful_count}")
    
    return results

def main():
    """Main function"""
    # Set paths
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    pipeline_dir = os.path.join(base_dir, "NO_highscore_result")
    
    csv_file_path = os.path.join(pipeline_dir, "zscore.csv")
    predicted_dir = os.path.join(pipeline_dir, "pdb_results")
    template_dir = os.path.join(pipeline_dir, "template")
    output_dir = os.path.join(pipeline_dir, "pymol_three_views_output")
    
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
    
    print("PyMOL three-view structure visualization based on CSV file")
    print("Generate three views:")
    print("  1. Protein chain comparison")
    print("  2. Protein chains + ATP")
    print("  3. ATP active site magnification (4Å range)")
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
        
        print(f"\nFinal statistics:")
        print(f"Z-score range: {min(z_scores):.2f} - {max(z_scores):.2f}")
        print(f"Average Z-score: {sum(z_scores)/len(z_scores):.2f}")
        print(f"RMSD range: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"Average RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"Structures containing nucleotides: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
        
        # Active site statistics
        active_sites = [r['active_site_atoms'] for r in results if r['active_site_atoms']]
        if active_sites:
            print(f"Active site atom count range: {min(active_sites)} - {max(active_sites)}")
            print(f"Average active site atom count: {sum(active_sites)/len(active_sites):.1f}")
    
    # Quit PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()