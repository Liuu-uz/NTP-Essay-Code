#!/usr/bin/env python3
"""
Interactive FASTA Uploader
Upload files and open interactive SSH session for manual analysis
Support for symlink-based processing
"""

import os
import sys
import subprocess

class InteractiveUploader:
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path
        self.fasta_filename = os.path.basename(fasta_file_path)
        self.sequence_id = os.path.splitext(self.fasta_filename)[0]
        
        # Fixed configuration - upload directly to target directory
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.upload_dir = "/mnt/data2/supfam/zhijing/fasta_sequence"  # Direct upload to target
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
        
        # Upload file directly to target directory
        print(f"📤 Uploading FASTA file to remote server...")
        print(f"   Target: {self.upload_dir}/{self.fasta_filename}")
        cmd = f"scp {self.fasta_file_path} {self.remote_server}:{self.upload_dir}/"
        result = subprocess.run(cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Upload failed")
            return False
        
        print(f"✅ File uploaded successfully: {self.fasta_filename}")
        return True
    
    def create_symlink_and_run(self):
        """Create symlink in SUPERFAMILY directory and run analysis"""
        print(f"🔗 Creating symlink and running SUPERFAMILY analysis...")
        
        # Commands to run on remote server
        commands = [
            # Change to SUPERFAMILY directory
            f"cd {self.superfamily_dir}",
            
            # Remove existing symlink if it exists
            f"rm -f {self.fasta_filename}",
            
            # Create new symlink
            f"ln -s {self.upload_dir}/{self.fasta_filename} {self.fasta_filename}",
            
            # Verify symlink was created
            f"ls -la {self.fasta_filename}",
            
            # Run SUPERFAMILY analysis
            f"./superfamily.pl {self.fasta_filename}"
        ]
        
        # Combine commands with && to ensure they run in sequence
        full_command = " && ".join(commands)
        ssh_cmd = f"ssh {self.remote_server} '{full_command}'"
        
        print(f"📋 Running commands on remote server:")
        for i, cmd in enumerate(commands, 1):
            print(f"   {i}. {cmd}")
        print()
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Remote analysis failed")
            return False
        
        print(f"✅ SUPERFAMILY analysis completed")
        return True
    
    def start_interactive_session(self):
        """Start interactive SSH session for manual verification"""
        print(f"\n🚀 Opening interactive SSH session for verification...")
        print("="*60)
        print(f"📋 File location: {self.upload_dir}/{self.fasta_filename}")
        print(f"📋 Symlink location: {self.superfamily_dir}/{self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*60)
        print("🔍 You can verify the results or run additional commands")
        print("📋 Available commands:")
        print(f"   ls -la {self.fasta_filename}        # Check symlink")
        print(f"   ls -la *.html                     # Check HTML results")
        print(f"   ls -la {self.sequence_id}.*       # Check all result files")
        print()
        print("⚠️  Type 'exit' to quit SSH session and download results")
        print("="*60)
        
        input("Press Enter to start SSH session...")
        
        # Start interactive SSH, directly entering SUPERFAMILY directory
        ssh_cmd = f"ssh -t {self.remote_server} 'cd {self.superfamily_dir}; bash'"
        
        print(f"\n🔗 Connecting to remote server...")
        subprocess.run(ssh_cmd, shell=True)
        
        print(f"\n✅ SSH session ended")
        return True
    
    def cleanup_symlink(self):
        """Clean up symlink after processing"""
        print(f"🧹 Cleaning up symlink...")
        
        cleanup_cmd = f"rm -f {self.superfamily_dir}/{self.fasta_filename}"
        ssh_cmd = f"ssh {self.remote_server} '{cleanup_cmd}'"
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode == 0:
            print(f"✅ Symlink cleaned up")
        else:
            print(f"⚠️  Failed to clean up symlink (may not exist)")
        
        return True
    
    def download_results(self):
        """Download HTML result files"""
        print(f"\n📥 Starting to download HTML result files...")
        
        # Create HTML results directory
        html_dir = "html_results"
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
            print(f"📁 Created HTML results directory: {html_dir}")
        
        # Try downloading from SUPERFAMILY directory
        html_files_to_try = [
            f"{self.sequence_id}.html",
            f"{self.sequence_id}_html",
            f"{self.fasta_filename}.html"
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
            print(f"   💡 Check if analysis completed successfully")
        
        return downloaded_files
    
    def process_file_automatic(self):
        """Automatic processing without interactive session"""
        print(f"🧬 Automatic FASTA Processor (with Symlink)")
        print(f"📋 File: {self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*50)
        
        # 1. Verify and upload file
        if not self.verify_and_upload():
            return False
        
        # 2. Create symlink and run analysis
        if not self.create_symlink_and_run():
            return False
        
        # 3. Download results
        downloaded_files = self.download_results()
        
        # 4. Clean up symlink
        self.cleanup_symlink()
        
        # 5. Summary
        print(f"\n{'='*60}")
        print("🎉 Automatic processing complete!")
        
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
    
    def process_file_interactive(self):
        """Interactive processing with manual verification"""
        print(f"🧬 Interactive FASTA Processor (with Symlink)")
        print(f"📋 File: {self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*50)
        
        # 1. Verify and upload file
        if not self.verify_and_upload():
            return False
        
        # 2. Create symlink and run analysis
        if not self.create_symlink_and_run():
            return False
        
        # 3. Start interactive session for verification
        if not self.start_interactive_session():
            return False
        
        # 4. Download results
        downloaded_files = self.download_results()
        
        # 5. Clean up symlink
        self.cleanup_symlink()
        
        # 6. Summary
        print(f"\n{'='*60}")
        print("🎉 Interactive processing complete!")
        
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

def process_all_files(mode='interactive'):
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
    
    mode_name = "Interactive" if mode == 'interactive' else "Automatic"
    print(f"🧬 Batch {mode_name} Processor (with Symlink)")
    print("="*50)
    print(f"📁 Found {len(fasta_files)} FASTA files:")
    for fasta_path in fasta_files:
        filename = os.path.basename(fasta_path)
        print(f"   - {filename}")
    
    # Show processing mode
    if mode == 'interactive':
        print(f"\n🔄 Mode: Interactive (SSH session for each file)")
    else:
        print(f"\n🔄 Mode: Automatic (no SSH sessions)")
    
    # Confirm processing
    confirm = input(f"\nProcess these files in {mode} mode? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 Operation cancelled")
        return False
    
    # Process one by one
    successful = 0
    for i, fasta_path in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {i}/{len(fasta_files)}")
        print(f"{'='*60}")
        
        uploader = InteractiveUploader(fasta_path)
        
        if mode == 'interactive':
            result = uploader.process_file_interactive()
        else:
            result = uploader.process_file_automatic()
        
        if result:
            successful += 1
        
        if i < len(fasta_files):
            if mode == 'interactive':
                input(f"\nPress Enter to continue to next file...")
            else:
                print(f"\n⏳ Waiting 2 seconds before next file...")
                import time
                time.sleep(2)
    
    print(f"\n🎉 Batch processing complete!")
    print(f"📊 Successfully processed: {successful}/{len(fasta_files)} files")
    return True

def main():
    if len(sys.argv) == 1:
        # Batch mode (default: interactive)
        process_all_files(mode='interactive')
    elif len(sys.argv) == 2:
        arg = sys.argv[1]
        if arg == '--auto':
            # Batch automatic mode
            process_all_files(mode='auto')
        elif arg == '--interactive':
            # Batch interactive mode
            process_all_files(mode='interactive')
        elif arg.endswith('.fa'):
            # Single file mode (interactive)
            uploader = InteractiveUploader(arg)
            uploader.process_file_interactive()
        else:
            print(f"❌ Invalid argument: {arg}")
            print_usage()
    elif len(sys.argv) == 3 and sys.argv[1] == '--file':
        # Single file with mode
        fasta_file = sys.argv[2]
        if fasta_file.endswith('.fa'):
            uploader = InteractiveUploader(fasta_file)
            uploader.process_file_interactive()
        else:
            print(f"❌ File must end with .fa: {fasta_file}")
    elif len(sys.argv) == 4 and sys.argv[1] == '--file':
        # Single file with specified mode
        mode = sys.argv[2]
        fasta_file = sys.argv[3]
        
        if not fasta_file.endswith('.fa'):
            print(f"❌ File must end with .fa: {fasta_file}")
            return
        
        uploader = InteractiveUploader(fasta_file)
        if mode == 'auto':
            uploader.process_file_automatic()
        elif mode == 'interactive':
            uploader.process_file_interactive()
        else:
            print(f"❌ Invalid mode: {mode}")
            print_usage()
    else:
        print_usage()

def print_usage():
    print("Usage:")
    print("  python interactive_uploader.py                        # Batch interactive mode")
    print("  python interactive_uploader.py --auto                 # Batch automatic mode")
    print("  python interactive_uploader.py --interactive          # Batch interactive mode")
    print("  python interactive_uploader.py <file.fa>              # Single file interactive")
    print("  python interactive_uploader.py --file <file.fa>       # Single file interactive")
    print("  python interactive_uploader.py --file auto <file.fa>  # Single file automatic")
    print()
    print("Examples:")
    print("  python interactive_uploader.py fasta_files/test.fa")
    print("  python interactive_uploader.py --auto")
    print("  python interactive_uploader.py --file auto test.fa")
    print()
    print("Modes:")
    print("  interactive: Upload → Run analysis → SSH session → Download results")
    print("  auto:        Upload → Run analysis → Download results (no SSH)")

if __name__ == "__main__":
    main()