#!/usr/bin/env python3
"""
Batch upload script - Upload all DALI generated JSON files to remote server
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import glob

class BatchFileUploader:
    def __init__(self):
        # Fixed path - DALI generated JSON folder
        self.json_folder = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/dali_sequences"
        
        # Remote server configuration
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_input_dir = "~/student/students_webserver/zhijing/input_jsons"
    
    def get_json_files(self):
        """Get all JSON files"""
        if not os.path.exists(self.json_folder):
            print(f"Error: JSON folder does not exist: {self.json_folder}")
            return []
        
        json_files = glob.glob(os.path.join(self.json_folder, "*.json"))
        
        if not json_files:
            print(f"Error: No JSON files found in folder: {self.json_folder}")
            return []
        
        print(f"Found {len(json_files)} JSON files")
        return json_files
    
    def verify_json_file(self, json_file):
        """Verify single JSON file"""
        if not os.path.exists(json_file):
            print(f"Error: File does not exist: {json_file}")
            return False
        
        if os.path.getsize(json_file) == 0:
            print(f"Error: File is empty: {json_file}")
            return False
        
        return True
    
    def upload_single_file(self, json_file):
        """Upload single JSON file"""
        filename = os.path.basename(json_file)
        
        print(f"Uploading: {filename}")
        
        cmd = f"scp '{json_file}' {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Upload failed for {filename}: {result.stderr}")
            return False
        
        print(f"Upload successful: {filename}")
        return True
    
    def upload_all_files(self):
        """Batch upload all JSON files"""
        print("Starting batch upload of DALI JSON files...")
        print(f"Source folder: {self.json_folder}")
        print(f"Target server: {self.remote_server}")
        print(f"Target path: {self.remote_input_dir}")
        print("="*60)
        
        # Get all JSON files
        json_files = self.get_json_files()
        if not json_files:
            return False
        
        # Statistics variables
        total_files = len(json_files)
        successful_uploads = 0
        failed_uploads = 0
        
        # Upload files one by one
        for i, json_file in enumerate(json_files, 1):
            filename = os.path.basename(json_file)
            print(f"[{i}/{total_files}] Processing: {filename}")
            
            # Verify file
            if not self.verify_json_file(json_file):
                failed_uploads += 1
                continue
            
            # Upload file
            if self.upload_single_file(json_file):
                successful_uploads += 1
            else:
                failed_uploads += 1
            
            print()  # Empty line separator
        
        # Output statistics
        print("="*60)
        print("Upload statistics:")
        print(f"   Total files: {total_files}")
        print(f"   Successful uploads: {successful_uploads}")
        print(f"   Failed uploads: {failed_uploads}")
        print("="*60)
        
        if successful_uploads > 0:
            print("Batch upload completed!")
            print(f"Next steps:")
            print(f"   1. SSH login to remote server: ssh {self.remote_server}")
            print(f"   2. Change to working directory: cd ~/student/students_webserver/zhijing")
            print(f"   3. Run remote scripts in batch (example):")
            
            # Show some example commands
            example_files = [os.path.basename(f) for f in json_files[:3]]
            for filename in example_files:
                print(f"      python remote_pipeline.py {filename}")
            
            if len(json_files) > 3:
                print(f"      ... (and {len(json_files)-3} more files)")
        
        return successful_uploads > 0

class SingleFileUploader:
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.json_filename = os.path.basename(json_file_path)
        self.output_name = os.path.splitext(self.json_filename)[0]
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_input_dir = "~/student/students_webserver/zhijing/input_jsons"
    
    def verify_json_file(self):
        """Verify if JSON file exists"""
        if not os.path.exists(self.json_file_path):
            print(f"Error: JSON file does not exist: {self.json_file_path}")
            return False
        
        print(f"JSON file verification successful: {self.json_file_path}")
        return True
    
    def upload_json_to_remote(self):
        """Upload JSON file to remote server"""
        print(f"Uploading JSON file to remote server...")
        
        cmd = f"scp '{self.json_file_path}' {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Upload failed: {result.stderr}")
            return False
        
        print(f"File upload successful: {self.json_filename}")
        return True
    
    def upload_file(self):
        """Main function for uploading single file"""
        print(f"Starting file upload...")
        print(f"JSON file: {self.json_file_path}")
        print(f"Filename: {self.json_filename}")
        print(f"Output name: {self.output_name}")
        print("="*50)
        
        # Verify file
        if not self.verify_json_file():
            return False
        
        # Upload file
        if not self.upload_json_to_remote():
            return False
        
        print("="*50)
        print("File upload completed!")
        print(f"Next steps:")
        print(f"   1. SSH login to remote server: ssh {self.remote_server}")
        print(f"   2. Change to working directory: cd ~/student/students_webserver/zhijing")
        print(f"   3. Run remote script: python remote_pipeline.py {self.json_filename}")
        print("="*50)
        
        return True

def main():
    parser = argparse.ArgumentParser(
        description='Upload JSON files to remote server',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  python upload_json.py --batch                    # Batch upload all DALI JSON files
  python upload_json.py single_file.json          # Upload single file
  python upload_json.py --file single_file.json   # Upload single file (parameter form)
        """
    )
    
    parser.add_argument('json_file', nargs='?', help='Single JSON file path')
    parser.add_argument('--file', help='Single JSON file path (parameter form)')
    parser.add_argument('--batch', action='store_true', 
                       help='Batch upload all DALI generated JSON files')
    
    args = parser.parse_args()
    
    if args.batch:
        # Batch upload mode
        uploader = BatchFileUploader()
        success = uploader.upload_all_files()
    else:
        # Single file upload mode
        json_file = args.json_file or args.file
        if not json_file:
            print("Error: Please specify JSON file path or use --batch for batch upload")
            print("Usage: python upload_json.py --batch")
            print("Or:    python upload_json.py your_file.json")
            sys.exit(1)
        
        uploader = SingleFileUploader(json_file)
        success = uploader.upload_file()
    
    if success:
        sys.exit(0)
    else:
        print("Upload failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()