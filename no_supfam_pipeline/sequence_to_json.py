#!/usr/bin/env python3
"""
Sequence to JSON generation script
Read sequence information from sequence_list.txt file and generate JSON files
"""

import json
import argparse
import sys
import os
from typing import Dict, List

class SequenceToJSONConverter:
    def __init__(self):
        pass
    
    def validate_sequence(self, sequence: str) -> bool:
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        cleaned_sequence = sequence.replace(' ', '').replace('\n', '').upper()
        return all(aa in valid_amino_acids for aa in cleaned_sequence)
    
    def clean_sequence(self, sequence: str) -> str:
        return sequence.replace(' ', '').replace('\n', '').replace('\r', '').upper().strip()
    
    def process_sequence_list_file(self, txt_file: str) -> List[Dict]:
        try:
            with open(txt_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            if not content:
                print("Error: File is empty")
                return []
            results = []
            if content.startswith('>'):
                results = self._process_fasta_format(content)
            elif '|' in content or ':' in content:
                results = self._process_named_sequences(content)
            else:
                results = self._process_simple_sequences(content)
            return results
        except FileNotFoundError:
            print(f"Error: File '{txt_file}' does not exist")
            return []
        except Exception as e:
            print(f"Error processing file: {e}")
            return []
    
    def _process_fasta_format(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        current_name = None
        current_sequence = ""
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_name and current_sequence:
                    results.append(self._create_sequence_entry(current_sequence, current_name))
                current_name = line[1:].strip()
                if not current_name:
                    current_name = f"sequence_{len(results) + 1}"
                current_sequence = ""
            elif line and not line.startswith('>'):
                current_sequence += line
        if current_name and current_sequence:
            results.append(self._create_sequence_entry(current_sequence, current_name))
        return results
    
    def _process_named_sequences(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line:
                continue
            name = None
            sequence = None
            if '|' in line:
                parts = line.split('|', 1)
                if len(parts) == 2:
                    name, sequence = parts
            elif ':' in line:
                parts = line.split(':', 1)
                if len(parts) == 2:
                    name, sequence = parts
            if name and sequence:
                name = name.strip()
                sequence = sequence.strip()
                if not name:
                    name = f"sequence_{i}"
                results.append(self._create_sequence_entry(sequence, name))
            else:
                print(f"Warning: Line {i} format incorrect, skipping: {line[:50]}...")
        return results
    
    def _process_simple_sequences(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line:
                continue
            name = f"sequence_{i}"
            results.append(self._create_sequence_entry(line, name))
        return results
    
    def _create_sequence_entry(self, sequence: str, name: str) -> Dict:
        cleaned_sequence = self.clean_sequence(sequence)
        if not cleaned_sequence:
            print(f"Warning: Sequence '{name}' is empty, skipping")
            return None
        if not self.validate_sequence(cleaned_sequence):
            print(f"Warning: Sequence '{name}' contains invalid characters, skipping")
            return None
        print(f"Processing sequence: {name} (length: {len(cleaned_sequence)})")
        return {
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": cleaned_sequence,
                        "count": 1
                    }
                }
            ],
            "name": name.replace(' ', '_').replace('/', '_')
        }
    
    def process_single_sequence(self, sequence: str, name: str = "user_sequence") -> Dict:
        return self._create_sequence_entry(sequence, name)
    
    def save_to_json_folder(self, data: List[Dict], output_folder: str = "sequences_json"):
        try:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
                print(f"Created folder: {output_folder}")
            saved_count = 0
            for entry in data:
                if entry is None:
                    continue
                filename = f"{entry['name']}.json"
                filepath = os.path.join(output_folder, filename)
                entry_data = [entry]
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(entry_data, f, indent=2, ensure_ascii=False)
                print(f"Saved file: {filepath}")
                saved_count += 1
            print(f"Successfully saved {saved_count} JSON files to folder: {output_folder}")
        except Exception as e:
            print(f"Error saving JSON files: {e}")
    
    def save_to_json(self, data: List[Dict], output_file: str):
        try:
            valid_data = [entry for entry in data if entry is not None]
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(valid_data, f, indent=2, ensure_ascii=False)
            print(f"Successfully saved to file: {output_file}")
        except Exception as e:
            print(f"Error saving JSON file: {e}")

def main():
    parser = argparse.ArgumentParser(description='Generate JSON files from sequence files')
    parser.add_argument('input_file', nargs='?', default='sequence_list.txt',
                       help='Input file path (default: sequence_list.txt)')
    parser.add_argument('-f', '--folder', default='sequences_json', 
                       help='Output JSON folder path (default: sequences_json)')
    parser.add_argument('-s', '--sequence', type=str,
                       help='Directly provide a sequence string')
    parser.add_argument('-n', '--name', type=str, default='user_sequence',
                       help='Sequence name (use with -s)')
    parser.add_argument('--single-file', action='store_true',
                       help='Save as a single JSON file instead of separate files in folder')
    
    args = parser.parse_args()
    
    converter = SequenceToJSONConverter()
    
    if args.sequence:
        print(f"Processing single sequence: {args.name}")
        sequence_data = converter.process_single_sequence(args.sequence, args.name)
        if sequence_data:
            if args.single_file:
                output_file = f"{args.folder}.json" if not args.folder.endswith('.json') else args.folder
                converter.save_to_json([sequence_data], output_file)
            else:
                converter.save_to_json_folder([sequence_data], args.folder)
            print("Processing completed!")
        else:
            print("Sequence processing failed")
    else:
        if not os.path.exists(args.input_file):
            print(f"Error: File '{args.input_file}' does not exist")
            print("\nSupported file formats:")
            print("1. One sequence per line:")
            print("   MKLAVIFG...")
            print("   MKLVQTSR...")
            print("\n2. name|sequence or name:sequence:")
            print("   protein1|MKLAVIFG...")
            print("   protein2:MKLVQTSR...")
            print("\n3. FASTA format:")
            print("   >protein1")
            print("   MKLAVIFG...")
            print("   >protein2") 
            print("   MKLVQTSR...")
            sys.exit(1)
        
        print(f"Starting to process file: {args.input_file}")
        sequences_data = converter.process_sequence_list_file(args.input_file)
        
        if sequences_data:
            if args.single_file:
                output_file = f"{args.folder}.json" if not args.folder.endswith('.json') else args.folder
                converter.save_to_json(sequences_data, output_file)
                print(f"Processing completed! Successfully processed {len(sequences_data)} sequences, saved to single file")
            else:
                converter.save_to_json_folder(sequences_data, args.folder)
                print(f"Processing completed! Successfully processed {len(sequences_data)} sequences, saved to separate files in folder")
        else:
            print("No valid sequences were processed")

if __name__ == "__main__":
    main()