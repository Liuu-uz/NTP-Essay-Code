#!/usr/bin/env python3
"""
FASTA File Splitter
Split multi-FASTA files into individual files organized by source filename
"""

import os
import sys
import re
from pathlib import Path

class FASTASplitter:
    def __init__(self, source_dir="/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_sequence", target_base_dir="/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"):
        self.source_dir = source_dir
        self.target_base_dir = target_base_dir
        self.stats = {
            "files_processed": 0,
            "empty_files": 0,
            "total_sequences": 0,
            "files_created": 0,
            "directories_created": 0
        }
    
    def clean_filename(self, header):
        """Extract clean filename from FASTA header"""
        # Remove '>' from header
        header = header.lstrip('>')
        
        # Try to extract protein ID or accession
        # Pattern 1: sp|ACCESSION|NAME or similar
        match = re.search(r'\|([A-Z0-9_]+)\|', header)
        if match:
            return match.group(1)
        
        # Pattern 2: First word before space
        parts = header.split()
        if parts:
            clean_name = re.sub(r'[^\w\-_.]', '_', parts[0])
            return clean_name
        
        # Fallback: use first 20 characters, cleaned
        clean_name = re.sub(r'[^\w\-_.]', '_', header[:20])
        return clean_name if clean_name else "sequence"
    
    def parse_fasta_file(self, file_path):
        """Parse FASTA file and return list of (header, sequence) tuples"""
        sequences = []
        current_header = None
        current_sequence = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    
                    if not line:  # Skip empty lines
                        continue
                    
                    if line.startswith('>'):  # Header line
                        # Save previous sequence if exists
                        if current_header is not None:
                            sequence_str = ''.join(current_sequence)
                            if sequence_str:  # Only save non-empty sequences
                                sequences.append((current_header, sequence_str))
                        
                        # Start new sequence
                        current_header = line
                        current_sequence = []
                    else:  # Sequence line
                        if current_header is not None:
                            current_sequence.append(line)
                
                # Don't forget the last sequence
                if current_header is not None:
                    sequence_str = ''.join(current_sequence)
                    if sequence_str:
                        sequences.append((current_header, sequence_str))
        
        except Exception as e:
            print(f"❌ Error reading {file_path}: {e}")
            return []
        
        return sequences
    
    def create_output_directory(self, source_filename):
        """Create output directory based on source filename"""
        # Extract base name without extension
        base_name = Path(source_filename).stem  # e.g., "1.16.99.1" from "1.16.99.1.fa"
        
        output_dir = os.path.join(self.target_base_dir, base_name)
        
        try:
            os.makedirs(output_dir, exist_ok=True)
            if not os.path.exists(output_dir):
                print(f"❌ Failed to create directory: {output_dir}")
                return None
            return output_dir
        except Exception as e:
            print(f"❌ Error creating directory {output_dir}: {e}")
            return None
    
    def write_individual_fasta(self, output_dir, filename, header, sequence):
        """Write individual FASTA file"""
        # Ensure filename has .fa extension
        if not filename.endswith('.fa'):
            filename += '.fa'
        
        output_path = os.path.join(output_dir, filename)
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(f"{header}\n")
                
                # Write sequence in lines of 80 characters (standard FASTA format)
                for i in range(0, len(sequence), 80):
                    f.write(f"{sequence[i:i+80]}\n")
            
            return output_path
        except Exception as e:
            print(f"❌ Error writing {output_path}: {e}")
            return None
    
    def process_single_file(self, file_path):
        """Process a single FASTA file"""
        filename = os.path.basename(file_path)
        print(f"\n🔍 Processing: {filename}")
        
        # Parse FASTA file
        sequences = self.parse_fasta_file(file_path)
        
        if not sequences:
            print(f"   ⚠️  Empty file or no valid sequences found")
            self.stats["empty_files"] += 1
            return False
        
        print(f"   📊 Found {len(sequences)} sequences")
        
        # Create output directory
        output_dir = self.create_output_directory(filename)
        if not output_dir:
            return False
        
        self.stats["directories_created"] += 1
        print(f"   📁 Output directory: {output_dir}")
        
        # Process each sequence
        created_files = []
        for i, (header, sequence) in enumerate(sequences, 1):
            # Generate filename from header
            clean_name = self.clean_filename(header)
            
            # Avoid filename conflicts
            base_filename = clean_name
            counter = 1
            while True:
                test_filename = f"{base_filename}.fa"
                test_path = os.path.join(output_dir, test_filename)
                if not os.path.exists(test_path):
                    break
                base_filename = f"{clean_name}_{counter}"
                counter += 1
            
            # Write individual file
            output_path = self.write_individual_fasta(output_dir, base_filename, header, sequence)
            
            if output_path:
                created_files.append(output_path)
                file_size = os.path.getsize(output_path)
                print(f"   ✅ Created: {os.path.basename(output_path)} ({file_size} bytes)")
            else:
                print(f"   ❌ Failed to create file for sequence {i}")
        
        self.stats["total_sequences"] += len(sequences)
        self.stats["files_created"] += len(created_files)
        
        print(f"   🎉 Successfully created {len(created_files)} files")
        return True
    
    def process_directory(self):
        """Process all FASTA files in source directory"""
        if not os.path.exists(self.source_dir):
            print(f"❌ Source directory does not exist: {self.source_dir}")
            return False
        
        # Find all .fa files
        fa_files = []
        for filename in os.listdir(self.source_dir):
            if filename.endswith('.fa'):
                file_path = os.path.join(self.source_dir, filename)
                fa_files.append(file_path)
        
        if not fa_files:
            print(f"❌ No .fa files found in {self.source_dir}")
            return False
        
        print(f"🧬 FASTA File Splitter")
        print("="*60)
        print(f"📂 Source directory: {self.source_dir}")
        print(f"📁 Target directory: {self.target_base_dir}")
        print(f"📋 Found {len(fa_files)} FASTA files to process")
        print("="*60)
        
        # Create base target directory
        try:
            os.makedirs(self.target_base_dir, exist_ok=True)
        except Exception as e:
            print(f"❌ Failed to create target directory {self.target_base_dir}: {e}")
            return False
        
        # Process each file
        for i, file_path in enumerate(fa_files, 1):
            print(f"\n[{i}/{len(fa_files)}] {'='*50}")
            success = self.process_single_file(file_path)
            if success:
                self.stats["files_processed"] += 1
        
        # Print summary
        self.print_summary()
        return True
    
    def process_single_file_direct(self, file_path):
        """Process a single file directly (for command line usage)"""
        if not os.path.exists(file_path):
            print(f"❌ File does not exist: {file_path}")
            return False
        
        print(f"🧬 FASTA File Splitter (Single File)")
        print("="*60)
        print(f"📄 Input file: {file_path}")
        print(f"📁 Target directory: {self.target_base_dir}")
        print("="*60)
        
        # Create base target directory
        try:
            os.makedirs(self.target_base_dir, exist_ok=True)
        except Exception as e:
            print(f"❌ Failed to create target directory {self.target_base_dir}: {e}")
            return False
        
        success = self.process_single_file(file_path)
        if success:
            self.stats["files_processed"] = 1
        
        self.print_summary()
        return success
    
    def print_summary(self):
        """Print processing summary"""
        print(f"\n{'='*60}")
        print("📊 Processing Summary")
        print(f"{'='*60}")
        print(f"Files processed:     {self.stats['files_processed']}")
        print(f"Empty files skipped: {self.stats['empty_files']}")
        print(f"Directories created: {self.stats['directories_created']}")
        print(f"Total sequences:     {self.stats['total_sequences']}")
        print(f"Files created:       {self.stats['files_created']}")
        print(f"Output location:     {self.target_base_dir}")
        
        if self.stats['files_created'] > 0:
            print(f"\n🎉 Successfully split {self.stats['total_sequences']} sequences into {self.stats['files_created']} individual files!")
        else:
            print(f"\n⚠️  No files were created. Check if input files contain valid FASTA sequences.")
        
        print("="*60)

def main():
    if len(sys.argv) == 1:
        # Process all files in default directory
        splitter = FASTASplitter()
        splitter.process_directory()
        
    elif len(sys.argv) == 2:
        arg = sys.argv[1]
        
        if arg == "--help":
            print("FASTA File Splitter")
            print("="*50)
            print("Split multi-FASTA files into individual files")
            print()
            print("Usage:")
            print("  python fasta_splitter.py                    # Process all .fa files in 'fasta_sequence' directory")
            print("  python fasta_splitter.py <file.fa>          # Process single file")
            print("  python fasta_splitter.py <directory>        # Process all .fa files in directory")
            print()
            print("Output structure:")
            print("  /Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files/")
            print("  ├── 1.16.99.1/")
            print("  │   ├── P65643_EUTT_ECOLI.fa")
            print("  │   ├── P65644_EUTT_ECOL6.fa")
            print("  │   └── ...")
            print("  ├── 2.3.1.193/")
            print("  │   └── ...")
            print("  └── ...")
            
        elif os.path.isfile(arg):
            # Process single file
            splitter = FASTASplitter()
            splitter.process_single_file_direct(arg)
            
        elif os.path.isdir(arg):
            # Process directory
            splitter = FASTASplitter(source_dir=arg)
            splitter.process_directory()
            
        else:
            print(f"❌ File or directory does not exist: {arg}")
    
    else:
        print("Usage: python fasta_splitter.py [file.fa|directory]")
        print("Use --help for more information")

if __name__ == "__main__":
    main()