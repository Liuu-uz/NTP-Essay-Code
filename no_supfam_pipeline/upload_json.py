#!/usr/bin/env python3
"""
Local upload script - Only responsible for uploading JSON files to remote server
"""

import os
import sys
import subprocess
import argparse

class FileUploader:
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.json_filename = os.path.basename(json_file_path)
        self.output_name = os.path.splitext(self.json_filename)[0]
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_input_dir = "~/student/students_webserver/zhijing/input_jsons"
        
    def verify_json_file(self):
        """Verify if JSON file exists"""
        if not os.path.exists(self.json_file_path):
            print(f"❌ JSON file does not exist: {self.json_file_path}")
            return False
        
        print(f"✅ JSON file verification successful: {self.json_file_path}")
        return True
    
    def upload_json_to_remote(self):
        """Upload JSON file to remote server"""
        print(f"📤 Uploading JSON file to remote server...")
        
        cmd = f"scp {self.json_file_path} {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ Upload failed: {result.stderr}")
            return False
        
        print(f"✅ File uploaded successfully: {self.json_filename}")
        return True
    
    def upload_file(self):
        """Main function for uploading files"""
        print(f"🚀 Starting file upload...")
        print(f"📋 JSON file: {self.json_file_path}")
        print(f"📋 Filename: {self.json_filename}")
        print(f"📋 Output name: {self.output_name}")
        print("="*50)
        
        # Verify file
        if not self.verify_json_file():
            return False
        
        # Upload file
        if not self.upload_json_to_remote():
            return False
        
        print("="*50)
        print("🎉 File upload completed!")
        print(f"📝 Next steps:")
        print(f"   1. SSH login to remote server: ssh {self.remote_server}")
        print(f"   2. Change to working directory: cd ~/student/students_webserver/zhijing")
        print(f"   3. Run remote script: python remote_pipeline.py {self.json_filename}")
        print("="*50)
        
        return True

def main():
    parser = argparse.ArgumentParser(description='Upload JSON files to remote server')
    parser.add_argument('json_file', help='JSON file path')
    
    args = parser.parse_args()
    
    uploader = FileUploader(args.json_file)
    success = uploader.upload_file()
    
    if success:
        sys.exit(0)
    else:
        print("💥 Upload failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()