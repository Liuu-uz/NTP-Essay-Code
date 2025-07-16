#!/usr/bin/env python3
"""
Interactive FASTA Uploader
Upload files and open interactive SSH session for manual analysis
"""

import os
import sys
import subprocess

class InteractiveUploader:
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path
        self.fasta_filename = os.path.basename(fasta_file_path)
        self.sequence_id = os.path.splitext(self.fasta_filename)[0]
        
        # Fixed configuration
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.upload_dir = "~/student/students_webserver/zhijing/fasta_sequence"
        self.superfamily_dir = "/mnt/data2/supfam/supfam"
        
    def verify_and_upload(self):
        """Verify and upload file"""
        if not os.path.exists(self.fasta_file_path):
            print(f"❌ FASTA file does not exist: {self.fasta_file_path}")
            return False
        
        if not self.fasta_filename.endswith('.fa'):
            print(f"❌ File must end with .fa: {self.fasta_filename}")
            return False
        
        print(f"✅ FASTA file verification successful: {self.fasta_file_path}")
        
        # Upload file
        print(f"📤 Uploading FASTA file to remote server...")
        cmd = f"scp {self.fasta_file_path} {self.remote_server}:{self.upload_dir}/"
        result = subprocess.run(cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Upload failed")
            return False
        
        print(f"✅ File uploaded successfully: {self.fasta_filename}")
        return True
    
    def copy_file_remote(self):
        """Copy file to SUPERFAMILY directory on remote server"""
        print(f"📋 Copying file to SUPERFAMILY directory...")
        
        copy_cmd = f"cp {self.upload_dir}/{self.fasta_filename} {self.superfamily_dir}/"
        ssh_cmd = f"ssh {self.remote_server} '{copy_cmd}'"
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Failed to copy file")
            return False
        
        print(f"✅ File copied to SUPERFAMILY directory")
        return True
    
    def start_interactive_session(self):
        """Start interactive SSH session"""
        print(f"\n🚀 Starting interactive SSH session...")
        print("="*60)
        print(f"📋 File is ready: {self.superfamily_dir}/{self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*60)
        print("🔄 SSH session will now open. Please run on the remote server:")
        print(f"   cd {self.superfamily_dir}")
        print(f"   ./superfamily.pl {self.fasta_filename}")
        print()
        print("⚠️  After analysis is complete, type 'exit' to quit SSH session")
        print("   Script will automatically download result files")
        print("="*60)
        
        input("Press Enter to start SSH session...")
        
        # Start interactive SSH, directly entering SUPERFAMILY directory
        ssh_cmd = f"ssh -t {self.remote_server} 'cd {self.superfamily_dir}; bash'"
        
        print(f"\n🔗 Connecting to remote server...")
        subprocess.run(ssh_cmd, shell=True)
        
        print(f"\n✅ SSH session ended")
        return True
    
    def download_results(self):
        """Download HTML result files"""
        print(f"\n📥 Starting to download HTML result files...")
        
        # Create HTML results directory
        html_dir = "html_results"
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
            print(f"📁 Created HTML results directory: {html_dir}")
        
        # Only download HTML files
        html_files_to_try = [
            f"{self.sequence_id}.html",
            f"{self.sequence_id}_html"
        ]
        
        downloaded_files = []
        
        for file_name in html_files_to_try:
            remote_path = f"{self.remote_server}:{self.superfamily_dir}/{file_name}"
            local_path = os.path.join(html_dir, f"{self.sequence_id}.html")  # Uniform naming as .html
            
            print(f"   Trying to download: {file_name}")
            cmd = f"scp {remote_path} {local_path}"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            
            if result.returncode == 0:
                downloaded_files.append(local_path)
                print(f"   ✅ HTML file downloaded successfully")
                print(f"      Saved to: {local_path}")
                
                # Check file size
                if os.path.exists(local_path):
                    file_size = os.path.getsize(local_path)
                    print(f"      File size: {file_size} bytes")
                    
                    if file_size > 500:
                        print(f"      ✅ HTML file appears complete")
                    else:
                        print(f"      ⚠️  HTML file may be incomplete")
                
                break  # Stop trying after successful download
            else:
                print(f"   ❌ {file_name} download failed (may not exist)")
        
        if not downloaded_files:
            print(f"   ⚠️  No HTML result files found")
        
        return downloaded_files
    
    def process_file(self):
        """Main function to process file"""
        print(f"🧬 Interactive FASTA Processor")
        print(f"📋 File: {self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*50)
        
        # 1. Verify and upload file
        if not self.verify_and_upload():
            return False
        
        # 2. Copy file to SUPERFAMILY directory
        if not self.copy_file_remote():
            return False
        
        # 3. Start interactive session
        if not self.start_interactive_session():
            return False
        
        # 4. Download results
        downloaded_files = self.download_results()
        
        # 5. Summary
        print(f"\n{'='*60}")
        print("🎉 Processing complete!")
        
        if downloaded_files:
            print(f"📋 HTML result files:")
            for file_path in downloaded_files:
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    print(f"   📄 {file_path} ({file_size} bytes)")
                    print(f"   🌐 Can be opened in browser to view results")
        else:
            print("⚠️  No HTML result files downloaded")
            print("💡 Please check if analysis completed successfully")
        
        print("="*60)
        return True

def process_all_files():
    """Batch process all files"""
    fasta_dir = "fasta_files"
    
    if not os.path.exists(fasta_dir):
        print(f"❌ Directory does not exist: {fasta_dir}")
        return False
    
    # Find all .fa files
    fasta_files = []
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fa'):
            fasta_path = os.path.join(fasta_dir, filename)
            fasta_files.append(fasta_path)
    
    if not fasta_files:
        print(f"❌ No .fa files found in {fasta_dir} directory")
        return False
    
    print(f"🧬 Batch Interactive Processor")
    print("="*50)
    print(f"📁 Found {len(fasta_files)} FASTA files:")
    for fasta_path in fasta_files:
        filename = os.path.basename(fasta_path)
        print(f"   - {filename}")
    
    # Confirm processing
    confirm = input(f"\nProcess these files one by one? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 Operation cancelled")
        return False
    
    # Process one by one
    for i, fasta_path in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {i}/{len(fasta_files)}")
        print(f"{'='*60}")
        
        uploader = InteractiveUploader(fasta_path)
        uploader.process_file()
        
        if i < len(fasta_files):
            input(f"\nPress Enter to continue to next file...")
    
    print(f"\n🎉 All files processed!")
    return True

def main():
    if len(sys.argv) == 1:
        # Batch mode
        process_all_files()
    elif len(sys.argv) == 2:
        # Single file mode
        fasta_file = sys.argv[1]
        uploader = InteractiveUploader(fasta_file)
        uploader.process_file()
    else:
        print("Usage:")
        print("  python interactive_uploader.py                # Batch process all files")
        print("  python interactive_uploader.py <file.fa>      # Process single file")
        print()
        print("Examples:")
        print("  python interactive_uploader.py fasta_files/test.fa")
        print("  python interactive_uploader.py")

if __name__ == "__main__":
    main()