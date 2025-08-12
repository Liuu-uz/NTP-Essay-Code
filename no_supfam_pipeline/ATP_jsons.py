#!/usr/bin/env python3
"""
NTP Processing Sequence Extractor - Ligand Version
Solve the overwrite problem of same sequence_id in different folders and add ATP ligand information
"""

import json
import sys
import os
import pandas as pd
from typing import Dict, List, Optional
from pathlib import Path
import re

class NTPProcessingSequenceExtractor:
    def __init__(self):
        # Fixed absolute paths
        self.excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report.xlsx"
        self.fasta_base_path = Path("/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files")
        
    def sanitize_filename(self, filename: str) -> str:
        """Clean special characters in filename"""
        # Replace characters not allowed by filesystem with underscore
        sanitized = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', filename)
        
        # Remove leading and trailing dots and spaces
        sanitized = sanitized.strip('. ')
        
        # Ensure not empty
        if not sanitized:
            sanitized = "unnamed"
        
        # Limit length
        if len(sanitized) > 200:
            sanitized = sanitized[:200]
        
        return sanitized
        
    def validate_sequence(self, sequence: str) -> bool:
        """Validate amino acid sequence - including extended IUPAC codes"""
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*')
        cleaned_sequence = sequence.replace(' ', '').replace('\n', '').upper()
        return all(aa in valid_amino_acids for aa in cleaned_sequence)
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean sequence string"""
        return sequence.replace(' ', '').replace('\n', '').replace('\r', '').upper().strip()
    
    def read_ntp_processing_records(self) -> List[Dict]:
        """Read NTP processing SF found records from Excel file"""
        try:
            print(f"Reading Excel file: {self.excel_path}")
            
            try:
                df = pd.read_excel(self.excel_path, sheet_name='Detailed_Results')
            except ValueError:
                df = pd.read_excel(self.excel_path, sheet_name=0)
                print("Using first sheet")
            
            print(f"Excel file has {len(df)} rows of data")
            print(f"Column names: {list(df.columns)}")
            
            # Filter records with Result_Category as 'NTP processing SF found'
            ntp_records = df[df['Result_Category'] == 'NTP processing SF found']
            print(f"Found {len(ntp_records)} 'NTP processing SF found' records")
            
            if len(ntp_records) == 0:
                unique_categories = df['Result_Category'].unique()
                print("Available Result_Category values:")
                for cat in unique_categories:
                    count = len(df[df['Result_Category'] == cat])
                    print(f"  - '{cat}': {count} records")
                
                print("\nTrying fuzzy match for records containing 'NTP'...")
                ntp_like_records = df[df['Result_Category'].str.contains('NTP', case=False, na=False)]
                print(f"Found {len(ntp_like_records)} records containing 'NTP'")
                
                if len(ntp_like_records) > 0:
                    ntp_records = ntp_like_records
                    print("Using records containing 'NTP'")
            
            # Check duplicate sequence_id
            sequence_ids = ntp_records['Sequence_ID'].astype(str).tolist()
            folders = ntp_records['Folder'].astype(str).tolist()
            
            print(f"Unique sequence_id count: {len(set(sequence_ids))}")
            print(f"Total records: {len(sequence_ids)}")
            
            # Count duplicates
            from collections import Counter
            id_counts = Counter(sequence_ids)
            duplicates = {seq_id: count for seq_id, count in id_counts.items() if count > 1}
            
            if duplicates:
                print(f"Found {len(duplicates)} duplicate sequence_ids:")
                for seq_id, count in list(duplicates.items())[:5]:  # Show only first 5
                    print(f"  - {seq_id}: {count} times")
                    # Show which folders these duplicate IDs come from
                    related_folders = [folders[i] for i, sid in enumerate(sequence_ids) if sid == seq_id]
                    print(f"    From folders: {set(related_folders)}")
                if len(duplicates) > 5:
                    print(f"  ... and {len(duplicates) - 5} more duplicate IDs")
            
            records = []
            for _, row in ntp_records.iterrows():
                records.append({
                    'folder': str(row['Folder']),
                    'sequence_id': str(row['Sequence_ID']),
                    'html_file': str(row['HTML_File']),
                    'found_superfamily': str(row['Found_Superfamily']) if pd.notna(row['Found_Superfamily']) else 'N/A',
                    'result_category': str(row['Result_Category'])
                })
            
            return records
            
        except FileNotFoundError:
            print(f"Error: Excel file not found - {self.excel_path}")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading Excel file: {e}")
            sys.exit(1)
    
    def find_fasta_file(self, folder: str, sequence_id: str) -> Optional[Path]:
        """Find corresponding FASTA file"""
        folder_path = self.fasta_base_path / folder
        
        if not folder_path.exists():
            print(f"Warning: Folder does not exist - {folder_path}")
            return None
        
        # Possible FASTA file extensions
        extensions = ['.fasta', '.fa', '.fas', '.seq', '.txt']
        
        # Exact match
        for ext in extensions:
            fasta_file = folder_path / f"{sequence_id}{ext}"
            if fasta_file.exists():
                return fasta_file
        
        # Fuzzy match
        for file in folder_path.iterdir():
            if file.is_file() and sequence_id in file.stem:
                return file
        
        print(f"Warning: FASTA file not found - {folder}/{sequence_id}")
        return None
    
    def read_fasta_sequence(self, fasta_file: Path) -> Optional[str]:
        """Read sequence from FASTA file"""
        try:
            with open(fasta_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            
            # Process FASTA format - skip lines starting with >
            lines = content.split('\n')
            sequence = ""
            
            for line in lines:
                line = line.strip()
                if not line.startswith('>') and line:
                    sequence += line
            
            return self.clean_sequence(sequence) if sequence else None
                
        except Exception as e:
            print(f"Error reading FASTA file {fasta_file}: {e}")
            return None
    
    def create_unique_name(self, folder: str, sequence_id: str) -> str:
        """Create unique sequence name combining sequence ID and folder"""
        # Scheme: sequenceID_folder
        unique_name = f"{sequence_id}_{folder}"
        
        # Clean special characters
        unique_name = self.sanitize_filename(unique_name)
        
        return unique_name
    
    def extract_ntp_processing_sequences(self) -> List[Dict]:
        """Extract all NTP processing SF found sequences including ATP ligand information"""
        ntp_records = self.read_ntp_processing_records()
        
        if not ntp_records:
            print("No NTP processing SF found records found")
            return []
        
        results = []
        success_count = 0
        seen_combinations = set()  # Track processed (folder, sequence_id) combinations
        duplicate_combinations = 0
        
        print(f"\nStarting to process {len(ntp_records)} NTP processing records...")
        
        for i, record in enumerate(ntp_records, 1):
            folder = record['folder']
            sequence_id = record['sequence_id']
            found_superfamily = record['found_superfamily']
            result_category = record['result_category']
            
            # Create unique identifier
            combination_key = (folder, sequence_id)
            
            print(f"[{i}/{len(ntp_records)}] Processing: {sequence_id} (Folder: {folder}, Category: {result_category})")
            
            # Check if this combination has been processed
            if combination_key in seen_combinations:
                duplicate_combinations += 1
                print(f"  -> Skipping duplicate combination: {folder}/{sequence_id}")
                continue
            
            seen_combinations.add(combination_key)
            
            # Find FASTA file
            fasta_file = self.find_fasta_file(folder, sequence_id)
            if not fasta_file:
                continue
            
            # Read sequence
            sequence = self.read_fasta_sequence(fasta_file)
            if not sequence:
                print(f"  -> Sequence is empty, skipping")
                continue
            
            # Validate sequence
            if not self.validate_sequence(sequence):
                print(f"  -> Sequence contains invalid characters, skipping")
                continue
            
            # Create unique name
            unique_name = self.create_unique_name(folder, sequence_id)
            
            # Create JSON entry including ligand information
            json_entry = {
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence,
                            "count": 1
                          }
                    },
                    {
                        "ligand": {
                            "ligand": "CCD_ATP",
                            "count": 1
                        }
                    }
                ],
                "name": unique_name,  # Use unique name
            }
            
            results.append(json_entry)
            success_count += 1
            print(f"  -> Success (Unique name: {unique_name}, Length: {len(sequence)}, Superfamily: {found_superfamily}, Ligand: ATP)")
        
        print(f"\nProcessing complete:")
        print(f"- Total records: {len(ntp_records)}")
        print(f"- Successfully processed: {success_count}")
        print(f"- Skipped duplicate combinations: {duplicate_combinations}")
        print(f"- Unique (folder, sequence_id) combinations: {len(seen_combinations)}")
        print(f"- All sequences include ATP ligand information")
        
        return results
    
    def save_sequences(self, sequences: List[Dict], output_path: str, single_file: bool = False):
        """Save sequences to JSON files with detailed tracking"""
        if not sequences:
            print("No sequences to save")
            return
        
        if single_file:
            # Save as single file
            output_file = Path(output_path)
            if not output_file.suffix:
                output_file = output_file.with_suffix('.json')
            
            try:
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(sequences, f, indent=2, ensure_ascii=False)
                print(f"Saved to single file: {output_file} ({len(sequences)} sequences)")
            except Exception as e:
                print(f"Error saving file: {e}")
        else:
            # Save as multiple files to folder
            output_dir = Path(output_path)
            output_dir.mkdir(exist_ok=True)
            
            saved_count = 0
            failed_count = 0
            failed_files = []
            seen_filenames = set()
            
            print(f"Starting to save {len(sequences)} sequences to {output_dir}")
            
            for i, sequence in enumerate(sequences):
                unique_name = sequence['name']
                original_id = sequence.get('original_sequence_id', 'unknown')
                folder = sequence.get('folder', 'unknown')
                
                # Ensure filename uniqueness
                base_filename = self.sanitize_filename(unique_name)
                filename = f"{base_filename}.json"
                
                # Add number if filename is duplicate
                counter = 1
                while filename in seen_filenames:
                    filename = f"{base_filename}_{counter}.json"
                    counter += 1
                
                seen_filenames.add(filename)
                filepath = output_dir / filename
                
                try:
                    with open(filepath, 'w', encoding='utf-8') as f:
                        json.dump([sequence], f, indent=2, ensure_ascii=False)
                    saved_count += 1
                    
                    # Report progress every 100 files
                    if (i + 1) % 100 == 0:
                        print(f"Progress: {i + 1}/{len(sequences)} saved")
                        
                except Exception as e:
                    failed_count += 1
                    failed_files.append({
                        'unique_name': unique_name, 
                        'original_id': original_id,
                        'folder': folder,
                        'error': str(e)
                    })
                    print(f"Error saving file {filename}: {e}")
            
            # Final report
            print(f"\nSaving complete statistics:")
            print(f"- Successfully saved: {saved_count} files")
            print(f"- Failed to save: {failed_count} files")
            print(f"- Save location: {output_dir}")
            print(f"- Each sequence includes ATP ligand information")
            
            # Show failed files
            if failed_files:
                print(f"\nFailed file details:")
                for fail_info in failed_files[:5]:  # Show only first 5
                    print(f"  - {fail_info['folder']}/{fail_info['original_id']}: {fail_info['error']}")
                if len(failed_files) > 5:
                    print(f"  ... and {len(failed_files) - 5} more failures")
            
            # Verify actual file count
            actual_files = list(output_dir.glob("*.json"))
            print(f"\nVerification: Directory actually contains {len(actual_files)} JSON files")
            
            if len(actual_files) != saved_count:
                print(f"Warning: Reported success count ({saved_count}) does not match actual files ({len(actual_files)})!")
            else:
                print("File saving verification passed!")

def main():
    # Fixed file paths
    excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report.xlsx"
    fasta_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
    output_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/ntp_processing_sequences_with_ligand"
    
    # Check if files exist
    if not os.path.exists(excel_path):
        print(f"Error: Excel file does not exist - {excel_path}")
        sys.exit(1)
    
    if not os.path.exists(fasta_path):
        print(f"Error: FASTA base path does not exist - {fasta_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("NTP Processing Sequence Extractor - Ligand Version")
    print("=" * 60)
    print(f"Excel file: {excel_path}")
    print(f"FASTA path: {fasta_path}")
    print(f"Output path: {output_path}")
    print(f"Ligand information: CCD_ATP (all sequences)")
    print("=" * 60)
    
    # Create extractor and process
    extractor = NTPProcessingSequenceExtractor()
    sequences = extractor.extract_ntp_processing_sequences()
    
    if sequences:
        extractor.save_sequences(sequences, output_path, single_file=False)
        print(f"\nProcessing complete! Successfully extracted {len(sequences)} NTP processing sequences")
        print(f"JSON files saved to: {output_path}")
        print(f"File naming format: sequenceID_folder.json")
        print(f"Each sequence includes ATP ligand information")
    else:
        print("\nNo sequences were successfully extracted")

if __name__ == "__main__":
    main()