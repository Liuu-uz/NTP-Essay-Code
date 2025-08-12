#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA File Generator
Function: Convert protein sequences to standard FASTA format files
"""

import os
import re
from datetime import datetime

def create_fasta_file(sequence, sequence_id=None, description="", output_dir="fasta_files"):
    """
    Create FASTA format file
    
    Parameters:
    - sequence: protein sequence
    - sequence_id: sequence ID (optional)
    - description: sequence description (optional)
    - output_dir: output directory
    
    Returns:
    - fasta_path: path of generated FASTA file
    - clean_seq_id: cleaned sequence ID
    """
    
    print("Starting FASTA file generation...")
    
    # 1. Clean sequence (remove spaces, newlines, non-amino acid characters)
    print("Cleaning sequence...")
    clean_sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    
    if len(clean_sequence) == 0:
        raise ValueError("Error: Sequence is empty or contains no valid amino acids")
    
    print(f"   Original sequence length: {len(sequence)} characters")
    print(f"   Cleaned sequence length: {len(clean_sequence)} amino acids")
    
    # 2. Generate sequence ID
    if sequence_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sequence_id = f"seq_{timestamp}"
        print(f"Auto-generated sequence ID: {sequence_id}")
    else:
        print(f"Using provided sequence ID: {sequence_id}")
    
    # 3. Clean sequence ID (keep only alphanumeric and underscores, superfamily.pl compliant)
    clean_seq_id = re.sub(r'[^a-zA-Z0-9_]', '_', sequence_id)
    if clean_seq_id != sequence_id:
        print(f"Sequence ID cleaned: {sequence_id} -> {clean_seq_id}")
    
    # 4. Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # 5. Create FASTA content
    fasta_content = f">{clean_seq_id}"
    if description:
        fasta_content += f" {description}"
    fasta_content += f"\n{clean_sequence}\n"
    
    # 6. Save to file (use .fa extension, superfamily.pl compliant)
    fasta_filename = f"{clean_seq_id}.fa"
    fasta_path = os.path.join(output_dir, fasta_filename)
    
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    print(f"FASTA file created: {fasta_path}")
    print(f"   Sequence ID: {clean_seq_id}")
    print(f"   Sequence length: {len(clean_sequence)} amino acids")
    print(f"   File size: {os.path.getsize(fasta_path)} bytes")
    
    return fasta_path, clean_seq_id

def batch_create_fasta(sequences_list, output_dir="fasta_files"):
    """
    Batch create FASTA files
    
    Parameters:
    - sequences_list: list of sequence info, each element contains sequence, sequence_id, description
    - output_dir: output directory
    
    Returns:
    - created_files: list of created files
    """
    
    print(f"Starting batch generation of {len(sequences_list)} FASTA files...")
    
    created_files = []
    
    for i, seq_info in enumerate(sequences_list, 1):
        print(f"\n[{i}/{len(sequences_list)}] Processing sequence: {seq_info.get('sequence_id', 'unknown')}")
        
        try:
            fasta_path, clean_seq_id = create_fasta_file(
                sequence=seq_info['sequence'],
                sequence_id=seq_info.get('sequence_id'),
                description=seq_info.get('description', ''),
                output_dir=output_dir
            )
            
            created_files.append({
                'sequence_id': clean_seq_id,
                'fasta_path': fasta_path,
                'original_info': seq_info
            })
            
        except Exception as e:
            print(f"Error processing sequence: {e}")
            continue
    
    print(f"\nBatch generation complete! Successfully created {len(created_files)} files")
    return created_files

def validate_fasta_file(fasta_path):
    """
    Validate FASTA file format
    """
    print(f"Validating FASTA file: {fasta_path}")
    
    if not os.path.exists(fasta_path):
        print("Error: File does not exist")
        return False
    
    try:
        with open(fasta_path, 'r') as f:
            content = f.read()
        
        lines = content.strip().split('\n')
        
        # Check if first line starts with >
        if not lines[0].startswith('>'):
            print("Error: First line is not a valid FASTA header")
            return False
        
        # Check sequence lines
        sequence_lines = lines[1:]
        sequence = ''.join(sequence_lines)
        
        # Check if contains valid amino acids
        valid_aa = re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence)
        if not valid_aa:
            print("Error: Sequence contains invalid characters")
            return False
        
        print(f"FASTA file validation passed")
        print(f"   Header: {lines[0]}")
        print(f"   Sequence length: {len(sequence)} amino acids")
        return True
        
    except Exception as e:
        print(f"Validation failed: {e}")
        return False

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA File Generator - Auto-run Version
Add your sequences to the SEQUENCES_TO_PROCESS list below and run directly
"""

import os
import re
from datetime import datetime

# ==================== Sequence Configuration Area ====================
# Add sequences you want to generate FASTA for here
SEQUENCES_TO_PROCESS = [
    {
        "sequence": """
        MMKKIDVKILDPRVGKEFPLPTYATSGSAGLDLRACLNDAVELAPGDTTLVPTGLAIHIA
        DPSLAAMMLPRSGLGHKHGIVLGNLVGLIDSDYQGQLMISVWNRGQDSFTIQPGERIAQM
        IFVPVVQAEFNLVEDFDFRQ
        """,
        "sequence_id": "test",
        "description": "test"
    },
    
    # Add your sequences here - copy the format above
    {
        "sequence": """
        Paste your protein sequence here
        Can include newlines and spaces, program will clean automatically
        """,
        "sequence_id": "my_protein_001",
        "description": "My test protein"
    }
    
    # Continue adding more sequences...
]

# Output directory setting
OUTPUT_DIR = "fasta_files"

# ==================== Processing Functions ====================

def create_fasta_file(sequence, sequence_id=None, description="", output_dir=OUTPUT_DIR):
    """
    Create FASTA format file
    """
    
    # 1. Clean sequence (remove spaces, newlines, non-amino acid characters)
    clean_sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    
    if len(clean_sequence) == 0:
        raise ValueError("Error: Sequence is empty or contains no valid amino acids")
    
    # 2. Generate sequence ID
    if sequence_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sequence_id = f"seq_{timestamp}"
    
    # 3. Clean sequence ID (keep only alphanumeric and underscores, superfamily.pl compliant)
    clean_seq_id = re.sub(r'[^a-zA-Z0-9_]', '_', sequence_id)
    
    # 4. Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 5. Create FASTA content
    fasta_content = f">{clean_seq_id}"
    if description:
        fasta_content += f" {description}"
    fasta_content += f"\n{clean_sequence}\n"
    
    # 6. Save to file (use .fa extension, superfamily.pl compliant)
    fasta_filename = f"{clean_seq_id}.fa"
    fasta_path = os.path.join(output_dir, fasta_filename)
    
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    return fasta_path, clean_seq_id, len(clean_sequence)

def auto_generate_fasta_files():
    """
    Auto-generate all configured FASTA files
    """
    print("FASTA File Auto-Generator")
    print("=" * 50)
    
    # Filter valid sequences
    valid_sequences = []
    for seq_info in SEQUENCES_TO_PROCESS:
        if seq_info['sequence'].strip() and "Paste your" not in seq_info['sequence']:
            valid_sequences.append(seq_info)
        else:
            print(f"Warning: Skipping empty sequence or template: {seq_info.get('sequence_id', 'unknown')}")
    
    if not valid_sequences:
        print("Error: No valid sequences found")
        print("Tip: Please add your sequences to SEQUENCES_TO_PROCESS at the top of the script")
        return []
    
    print(f"Starting generation of {len(valid_sequences)} FASTA files...")
    
    created_files = []
    
    for i, seq_info in enumerate(valid_sequences, 1):
        print(f"\n[{i}/{len(valid_sequences)}] Processing sequence: {seq_info['sequence_id']}")
        
        try:
            fasta_path, clean_seq_id, seq_length = create_fasta_file(
                sequence=seq_info['sequence'],
                sequence_id=seq_info['sequence_id'],
                description=seq_info.get('description', '')
            )
            
            print(f"Generation successful: {fasta_path}")
            print(f"   Sequence ID: {clean_seq_id}")
            print(f"   Sequence length: {seq_length} amino acids")
            
            created_files.append({
                'sequence_id': clean_seq_id,
                'fasta_path': fasta_path,
                'length': seq_length
            })
            
        except Exception as e:
            print(f"Generation failed: {e}")
            continue
    
    # Summary report
    print(f"\n{'='*50}")
    print("Generation Complete Statistics")
    print(f"{'='*50}")
    print(f"Configured sequences: {len(SEQUENCES_TO_PROCESS)}")
    print(f"Valid sequences: {len(valid_sequences)}")
    print(f"Successfully generated: {len(created_files)}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    if created_files:
        print(f"\nGenerated file list:")
        for file_info in created_files:
            print(f"  - {file_info['fasta_path']} ({file_info['length']} aa)")
        
        print(f"\nNext step:")
        print(f"   Use remote_superfamily_runner.py to upload and analyze these files")
    
    return created_files

def validate_all_files():
    """
    Validate all generated FASTA files
    """
    print(f"\nValidating FASTA files in {OUTPUT_DIR} directory...")
    
    if not os.path.exists(OUTPUT_DIR):
        print(f"Error: Directory does not exist: {OUTPUT_DIR}")
        return
    
    fasta_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.fa')]
    
    if not fasta_files:
        print(f"Error: No .fa files found in directory")
        return
    
    print(f"Found {len(fasta_files)} FASTA files:")
    
    for fasta_file in fasta_files:
        fasta_path = os.path.join(OUTPUT_DIR, fasta_file)
        
        try:
            with open(fasta_path, 'r') as f:
                content = f.read().strip()
            
            lines = content.split('\n')
            header = lines[0] if lines else ""
            sequence = ''.join(lines[1:]) if len(lines) > 1 else ""
            
            if header.startswith('>') and sequence:
                valid_aa = re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence)
                if valid_aa:
                    print(f"  Valid: {fasta_file} - {len(sequence)} aa")
                else:
                    print(f"  Warning: {fasta_file} - contains invalid characters")
            else:
                print(f"  Error: {fasta_file} - format error")
                
        except Exception as e:
            print(f"  Error: {fasta_file} - read failed: {e}")

def list_generated_files():
    """
    List generated files for other scripts to call
    """
    if not os.path.exists(OUTPUT_DIR):
        return []
    
    fasta_files = []
    for filename in os.listdir(OUTPUT_DIR):
        if filename.endswith('.fa'):
            fasta_path = os.path.join(OUTPUT_DIR, filename)
            sequence_id = filename.replace('.fa', '')
            fasta_files.append({
                'sequence_id': sequence_id,
                'fasta_path': fasta_path,
                'filename': filename
            })
    
    return fasta_files

def main():
    """
    Main function - auto-run
    """
    # Auto-generate FASTA files
    created_files = auto_generate_fasta_files()
    
    # Validate generated files
    if created_files:
        validate_all_files()
    
    return created_files

if __name__ == "__main__":
    # Run directly, no user interaction needed
    main()