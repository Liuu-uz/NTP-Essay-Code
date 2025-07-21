#!/usr/bin/env python3
"""
Improved Interactive FASTA Uploader
Upload files to organized remote directories and run SUPERFAMILY analysis
"""

import os
import sys
import subprocess

class ImprovedUploader:
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path
        self.fasta_filename = os.path.basename(fasta_file_path)
        self.sequence_id = os.path.splitext(self.fasta_filename)[0]
        
        # Get the parent folder name (e.g., "1.16.99.1" from "fasta_files/1.16.99.1/protein.fa")
        self.local_parent_dir = os.path.basename(os.path.dirname(fasta_file_path))
        
        # Remote server configuration
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.base_upload_dir = "/mnt/data2/supfam/zhijing/fasta_sequence"
        self.superfamily_dir = "/mnt/data2/supfam/supfam"
        
        # Create remote directory path matching local structure
        self.remote_dir = f"{self.base_upload_dir}/{self.local_parent_dir}"
        self.remote_file_path = f"{self.remote_dir}/{self.fasta_filename}"
        
    def create_remote_directory(self):
        """Create remote directory structure matching local structure"""
        print(f"📁 Creating remote directory structure...")
        print(f"   Local parent folder: {self.local_parent_dir}")
        print(f"   Remote directory: {self.remote_dir}")
        
        # Create directory on remote server
        create_cmd = f"mkdir -p {self.remote_dir}"
        ssh_cmd = f"ssh {self.remote_server} '{create_cmd}'"
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Failed to create remote directory")
            return False
        
        print(f"✅ Remote directory created: {self.remote_dir}")
        return True
    
    def verify_and_upload(self):
        """Verify local file and upload to remote directory"""
        if not os.path.exists(self.fasta_file_path):
            print(f"❌ FASTA file does not exist: {self.fasta_file_path}")
            return False
        
        if not self.fasta_filename.endswith('.fa'):
            print(f"❌ File must end with .fa: {self.fasta_filename}")
            return False
        
        print(f"✅ FASTA file verification successful: {self.fasta_file_path}")
        
        # Upload file to the created directory
        print(f"📤 Uploading FASTA file to remote server...")
        print(f"   Target: {self.remote_file_path}")
        
        cmd = f"scp {self.fasta_file_path} {self.remote_server}:{self.remote_file_path}"
        result = subprocess.run(cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Upload failed")
            return False
        
        print(f"✅ File uploaded successfully to: {self.remote_file_path}")
        return True
    
    def create_symlink_and_manual_run(self):
        """Create symlink and provide instructions for manual SUPERFAMILY run"""
        print(f"🔗 Setting up symlink for SUPERFAMILY analysis...")
        
        # Step 1: Change to SUPERFAMILY directory and create symlink
        print(f"📋 Creating symlink...")
        setup_commands = [
            f"cd {self.superfamily_dir}",
            f"rm -f {self.fasta_filename}",  # Remove existing symlink if exists
            f"ln -s {self.remote_file_path} {self.fasta_filename}"  # Create symlink
        ]
        
        setup_command = " && ".join(setup_commands)
        ssh_setup_cmd = f"ssh {self.remote_server} '{setup_command}'"
        
        print(f"   Creating symlink: {self.fasta_filename} -> {self.remote_file_path}")
        result = subprocess.run(ssh_setup_cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ Symlink creation failed")
            return False
        
        print(f"✅ Symlink created successfully")
        return True
    
    def show_manual_instructions(self):
        """Show instructions for manual SUPERFAMILY analysis"""
        print(f"\n📋 Manual SUPERFAMILY Analysis Required")
        print("="*60)
        print(f"🔗 SSH Connection Command:")
        print(f"   ssh {self.remote_server}")
        print()
        print(f"📂 Commands to run on remote server:")
        print(f"   cd {self.superfamily_dir}")
        print(f"   ls -la {self.fasta_filename}  # Verify symlink exists")
        print(f"   ./superfamily.pl {self.fasta_filename}")
        print()
        print(f"📄 File being analyzed:")
        print(f"   Filename: {self.fasta_filename}")
        print(f"   Sequence ID: {self.sequence_id}")
        print(f"   Remote path: {self.remote_file_path}")
        print()
        print(f"⚠️  Expected result files after analysis:")
        print(f"   {self.sequence_id}.html")
        print(f"   {self.sequence_id}.ass")
        print(f"   {self.sequence_id}.*")
        print()
        print(f"🔄 After analysis completes:")
        print(f"   1. Check that {self.sequence_id}.html was created")
        print(f"   2. Return to this terminal")
        print(f"   3. Press Enter to continue with download")
        print("="*60)
        
        input("⏸️  Run the analysis manually, then press Enter to continue...")
        return True
    
    def download_results(self):
        """Download HTML result files to organized folders"""
        print(f"\n📥 Starting to download HTML result files...")
        
        # Create HTML results directory structure matching the folder structure
        base_html_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/html_results"
        folder_html_dir = os.path.join(base_html_dir, self.local_parent_dir)
        
        # Create the folder-specific directory
        if not os.path.exists(folder_html_dir):
            try:
                os.makedirs(folder_html_dir, exist_ok=True)
                print(f"📁 Created HTML results folder: {folder_html_dir}")
            except Exception as e:
                print(f"❌ Failed to create directory {folder_html_dir}: {e}")
                # Fallback to base html_results directory
                if not os.path.exists(base_html_dir):
                    os.makedirs(base_html_dir, exist_ok=True)
                folder_html_dir = base_html_dir
                print(f"📁 Using fallback directory: {folder_html_dir}")
        else:
            print(f"📁 Using existing HTML results folder: {folder_html_dir}")
        
        # Try downloading from SUPERFAMILY directory
        html_files_to_try = [
            f"{self.sequence_id}.html",
            f"{self.sequence_id}_html",
            f"{self.fasta_filename}.html"
        ]
        
        downloaded_files = []
        
        for file_name in html_files_to_try:
            remote_path = f"{self.remote_server}:{self.superfamily_dir}/{file_name}"
            local_path = os.path.join(folder_html_dir, f"{self.sequence_id}.html")  # Uniform naming
            
            print(f"   Trying to download: {file_name}")
            print(f"   Target path: {local_path}")
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
    
    def cleanup_remote_files(self):
        """Clean up only SUPERFAMILY result files, keep uploaded files"""
        print(f"🧹 Cleaning up SUPERFAMILY result files...")
        
        # Commands to clean up only SUPERFAMILY results
        cleanup_commands = [
            f"rm -f {self.superfamily_dir}/{self.fasta_filename}",  # Remove symlink
            f"rm -f {self.superfamily_dir}/{self.sequence_id}.*"    # Remove result files
        ]
        
        for i, cmd in enumerate(cleanup_commands, 1):
            desc = "symlink" if "ln -s" not in cmd and self.fasta_filename in cmd else "result files"
            print(f"   {i}. Removing {desc}")
            ssh_cmd = f"ssh {self.remote_server} '{cmd}'"
            result = subprocess.run(ssh_cmd, shell=True, capture_output=True)
            
            if result.returncode == 0:
                print(f"      ✅ Completed")
            else:
                print(f"      ⚠️  May not exist (ok)")
        
        print(f"✅ SUPERFAMILY cleanup completed")
        print(f"💾 Uploaded files preserved at: {self.remote_file_path}")
        return True
    
    def process_file_automatic(self):
        """Semi-automatic processing with manual SUPERFAMILY step"""
        print(f"🧬 Semi-Automatic FASTA Processor")
        print(f"📋 File: {self.fasta_filename}")
        print(f"📋 Local folder: {self.local_parent_dir}")
        print(f"📋 Remote directory: {self.remote_dir}")
        print("="*50)
        
        # 1. Create remote directory
        if not self.create_remote_directory():
            return False
        
        # 2. Verify and upload file
        if not self.verify_and_upload():
            return False
        
        # 3. Create symlink and show manual instructions
        if not self.create_symlink_and_manual_run():
            return False
        
        # 4. Show manual instructions
        if not self.show_manual_instructions():
            return False
        
        # 5. Download results
        downloaded_files = self.download_results()
        
        # 6. Clean up remote files
        self.cleanup_remote_files()
        
        # 7. Summary
        print(f"\n{'='*60}")
        print("🎉 Semi-automatic processing complete!")
        
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
    
    def start_interactive_session(self):
        """Start interactive SSH session for manual analysis"""
        print(f"\n🚀 Opening interactive SSH session...")
        print("="*60)
        print(f"📋 Remote file: {self.remote_file_path}")
        print(f"📋 Symlink location: {self.superfamily_dir}/{self.fasta_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print("="*60)
        print("🔄 Commands to run in SSH session:")
        print(f"   cd {self.superfamily_dir}")
        print(f"   ls -la {self.fasta_filename}        # Verify symlink")
        print(f"   ./superfamily.pl {self.fasta_filename}  # Run analysis")
        print(f"   ls -la {self.sequence_id}.*         # Check results")
        print()
        print("⚠️  After analysis completes, type 'exit' to quit SSH session")
        print("="*60)
        
        input("Press Enter to start SSH session...")
        
        # Start interactive SSH, directly entering SUPERFAMILY directory
        ssh_cmd = f"ssh -t {self.remote_server} 'cd {self.superfamily_dir}; bash'"
        
        print(f"\n🔗 Connecting to remote server...")
        subprocess.run(ssh_cmd, shell=True)
        
        print(f"\n✅ SSH session ended")
        return True
    
    def process_file_interactive(self):
        """Interactive processing with manual verification"""
        print(f"🧬 Improved Interactive FASTA Processor")
        print(f"📋 File: {self.fasta_filename}")
        print(f"📋 Local folder: {self.local_parent_dir}")
        print(f"📋 Remote directory: {self.remote_dir}")
        print("="*50)
        
        # 1. Create remote directory
        if not self.create_remote_directory():
            return False
        
        # 2. Verify and upload file
        if not self.verify_and_upload():
            return False
        
        # 3. Create symlink and show manual instructions
        if not self.create_symlink_and_manual_run():
            return False
        
        # 4. Show manual instructions
        if not self.show_manual_instructions():
            return False
        
        # 4. Start interactive session for verification
        if not self.start_interactive_session():
            return False
        
        # 5. Download results
        downloaded_files = self.download_results()
        
        # 6. Clean up remote files
        self.cleanup_remote_files()
        
        # 7. Summary
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

def select_and_process_folder(base_directory, mode='interactive'):
    """List folders and let user select one to process"""
    if not os.path.exists(base_directory):
        print(f"❌ Directory does not exist: {base_directory}")
        return False
    
    # Find all subdirectories that contain .fa files
    folders_with_fa = []
    
    for item in os.listdir(base_directory):
        item_path = os.path.join(base_directory, item)
        if os.path.isdir(item_path):
            # Check if this directory contains .fa files
            fa_files = [f for f in os.listdir(item_path) if f.endswith('.fa')]
            if fa_files:
                folders_with_fa.append({
                    'name': item,
                    'path': item_path,
                    'fa_count': len(fa_files)
                })
    
    if not folders_with_fa:
        print(f"❌ No folders containing .fa files found in {base_directory}")
        return False
    
    # Display folders for selection
    print(f"🧬 Folder Selection for FASTA Processing")
    print("="*60)
    print(f"📂 Base directory: {base_directory}")
    print(f"📁 Found {len(folders_with_fa)} folders containing .fa files:")
    print()
    
    for i, folder in enumerate(folders_with_fa, 1):
        print(f"  {i:2d}. {folder['name']:<20} ({folder['fa_count']} .fa files)")
    
    print()
    print("  0. Exit")
    print("="*60)
    
    # Get user selection
    while True:
        try:
            choice = input(f"Select folder to process (1-{len(folders_with_fa)}, 0 to exit): ").strip()
            
            if choice == '0':
                print("👋 Operation cancelled")
                return False
            
            choice_num = int(choice)
            if 1 <= choice_num <= len(folders_with_fa):
                selected_folder = folders_with_fa[choice_num - 1]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(folders_with_fa)}")
                
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Show selected folder details
    selected_path = selected_folder['path']
    fa_files = [f for f in os.listdir(selected_path) if f.endswith('.fa')]
    
    print(f"\n✅ Selected folder: {selected_folder['name']}")
    print(f"📂 Path: {selected_path}")
    print(f"📋 Contains {len(fa_files)} .fa files:")
    
    # Show first few files
    for i, fa_file in enumerate(fa_files[:5], 1):
        print(f"   {i}. {fa_file}")
    
    if len(fa_files) > 5:
        print(f"   ... and {len(fa_files) - 5} more files")
    
    mode_name = "Interactive" if mode == 'interactive' else "Automatic"
    print(f"\n🔄 Processing mode: {mode_name}")
    
    # Final confirmation
    confirm = input(f"\nProcess folder '{selected_folder['name']}' with {len(fa_files)} files? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 Operation cancelled")
        return False
    
    # Process the selected folder
    return process_directory_files(selected_path, mode)

def process_directory_files(directory_path, mode='interactive'):
    """Process .fa files in a directory (max 5 files)"""
    if not os.path.exists(directory_path):
        print(f"❌ Directory does not exist: {directory_path}")
        return False
    
    # Find all .fa files in the directory
    all_fasta_files = []
    for file in os.listdir(directory_path):
        if file.endswith('.fa'):
            file_path = os.path.join(directory_path, file)
            all_fasta_files.append(file_path)
    
    if not all_fasta_files:
        print(f"❌ No .fa files found in {directory_path}")
        return False
    
    # Sort files for consistent selection
    all_fasta_files.sort()
    
    # Limit to first 5 files
    fasta_files = all_fasta_files[:5]
    
    folder_name = os.path.basename(directory_path)
    mode_name = "Interactive" if mode == 'interactive' else "Automatic"
    
    print(f"\n🧬 Processing Folder: {folder_name}")
    print("="*50)
    print(f"📁 Directory: {directory_path}")
    print(f"📋 Found {len(all_fasta_files)} FASTA files total")
    
    if len(all_fasta_files) > 5:
        print(f"⚠️  Limiting to first 5 files (sorted alphabetically):")
    else:
        print(f"📋 Processing {len(fasta_files)} FASTA files:")
    
    for i, fasta_path in enumerate(fasta_files, 1):
        filename = os.path.basename(fasta_path)
        print(f"   {i:2d}. {filename}")
    
    if len(all_fasta_files) > 5:
        skipped_count = len(all_fasta_files) - 5
        print(f"   ... ({skipped_count} files skipped)")
    
    print(f"\n🔄 Mode: {mode_name}")
    print("="*50)
    
    # Show all SUPERFAMILY commands that will need to be run
    print(f"\n📋 SUPERFAMILY Commands to Run Manually:")
    print("="*50)
    print(f"🔗 SSH Connection: ssh webserver@coulomb.phys.ucl.ac.uk")
    print(f"📂 Change directory: cd /mnt/data2/supfam/supfam")
    print(f"🚀 Run these commands in sequence:")
    print()
    
    for i, fasta_path in enumerate(fasta_files, 1):
        filename = os.path.basename(fasta_path)
        sequence_id = os.path.splitext(filename)[0]
        print(f"   {i}. ./superfamily.pl {filename}    # For {sequence_id}")
    
    print()
    print(f"💡 Tip: You can run them one by one, or create a batch script:")
    batch_script = " && ".join([f"./superfamily.pl {os.path.basename(f)}" for f in fasta_files])
    print(f"   {batch_script}")
    print("="*50)
    
    # Confirm processing
    confirm = input(f"\nProcess these {len(fasta_files)} files? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 Operation cancelled")
        return False
    
    # Process files one by one
    successful = 0
    failed = 0
    uploaded_files = []  # Track uploaded files for batch command display
    
    for i, fasta_path in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {i}/{len(fasta_files)} from folder: {folder_name}")
        print(f"File: {os.path.basename(fasta_path)}")
        print(f"{'='*60}")
        
        try:
            uploader = ImprovedUploader(fasta_path)
            
            # For automatic mode, just do setup
            if mode == 'auto':
                # 1. Create remote directory
                if not uploader.create_remote_directory():
                    failed += 1
                    continue
                
                # 2. Upload file
                if not uploader.verify_and_upload():
                    failed += 1
                    continue
                
                # 3. Create symlink
                if not uploader.create_symlink_and_manual_run():
                    failed += 1
                    continue
                
                uploaded_files.append(uploader)
                successful += 1
                print(f"✅ File {i} uploaded and ready for analysis")
                
            else:
                # Interactive mode - full process
                result = uploader.process_file_interactive()
                if result:
                    successful += 1
                    print(f"✅ File {i} processed successfully")
                else:
                    failed += 1
                    print(f"❌ File {i} processing failed")
        
        except Exception as e:
            failed += 1
            print(f"❌ File {i} processing failed with error: {e}")
    
    # For automatic mode, show batch processing instructions
    if mode == 'auto' and uploaded_files:
        print(f"\n{'='*60}")
        print(f"🎯 BATCH ANALYSIS READY")
        print(f"{'='*60}")
        print(f"📋 {len(uploaded_files)} files uploaded and symlinks created")
        print(f"🔗 SSH to server: ssh webserver@coulomb.phys.ucl.ac.uk")
        print(f"📂 Change directory: cd /mnt/data2/supfam/supfam")
        print()
        print(f"🚀 Run these commands for batch analysis:")
        for j, uploader in enumerate(uploaded_files, 1):
            print(f"   {j}. ./superfamily.pl {uploader.fasta_filename}")
        
        print(f"\n💡 Or run all at once:")
        batch_cmd = " && ".join([f"./superfamily.pl {u.fasta_filename}" for u in uploaded_files])
        print(f"   {batch_cmd}")
        
        print(f"\n⏸️  After all analyses complete, continue for download and cleanup...")
        input("Press Enter when all analyses are finished...")
        
        # Download results and cleanup for all files
        download_successful = 0
        for j, uploader in enumerate(uploaded_files, 1):
            print(f"\n[{j}/{len(uploaded_files)}] Downloading results for {uploader.fasta_filename}")
            downloaded_files = uploader.download_results()
            uploader.cleanup_remote_files()
            
            if downloaded_files:
                download_successful += 1
        
        print(f"\n📊 Batch Download Results:")
        print(f"   ✅ Downloaded: {download_successful}")
        print(f"   ❌ Failed: {len(uploaded_files) - download_successful}")
    
    # Final summary
    print(f"\n🎉 Folder '{folder_name}' processing complete!")
    print(f"📊 Results:")
    print(f"   ✅ Successful: {successful}")
    print(f"   ❌ Failed: {failed}")
    print(f"   📁 Processed: {len(fasta_files)}")
    if len(all_fasta_files) > 5:
        print(f"   ⏩ Skipped: {len(all_fasta_files) - 5}")
    print(f"   📈 Success rate: {successful/len(fasta_files)*100:.1f}%")
    
    return successful > 0
    """Process all .fa files in a directory"""
    if not os.path.exists(directory_path):
        print(f"❌ Directory does not exist: {directory_path}")
        return False
    
    # Find all .fa files recursively
    fasta_files = []
    
    # Check if directory contains subdirectories with .fa files
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.fa'):
                file_path = os.path.join(root, file)
                fasta_files.append(file_path)
    
    if not fasta_files:
        print(f"❌ No .fa files found in {directory_path}")
        return False
    
    mode_name = "Interactive" if mode == 'interactive' else "Automatic"
    print(f"🧬 Batch {mode_name} Processor (Improved)")
    print("="*50)
    print(f"📁 Found {len(fasta_files)} FASTA files in {directory_path}:")
    for fasta_path in fasta_files[:10]:  # Show first 10 files
        rel_path = os.path.relpath(fasta_path, directory_path)
        print(f"   - {rel_path}")
    
    if len(fasta_files) > 10:
        print(f"   ... and {len(fasta_files) - 10} more files")
    
    # Show processing mode
    if mode == 'interactive':
        print(f"\n🔄 Mode: Interactive (SSH session for each file)")
    else:
        print(f"\n🔄 Mode: Automatic (no SSH sessions)")
    
    # Confirm processing
    confirm = input(f"\nProcess these {len(fasta_files)} files in {mode} mode? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 Operation cancelled")
        return False
    
    # Process files one by one
    successful = 0
    failed = 0
    
    for i, fasta_path in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {i}/{len(fasta_files)}")
        print(f"File: {os.path.basename(fasta_path)}")
        print(f"{'='*60}")
        
        try:
            uploader = ImprovedUploader(fasta_path)
            
            if mode == 'interactive':
                result = uploader.process_file_interactive()
            else:
                result = uploader.process_file_automatic()
            
            if result:
                successful += 1
                print(f"✅ File {i} processed successfully")
            else:
                failed += 1
                print(f"❌ File {i} processing failed")
        
        except Exception as e:
            failed += 1
            print(f"❌ File {i} processing failed with error: {e}")
        
        # Wait between files
        if i < len(fasta_files):
            if mode == 'interactive':
                input(f"\nPress Enter to continue to next file...")
            else:
                print(f"\n⏳ Waiting 3 seconds before next file...")
                import time
                time.sleep(3)
    
    # Final summary
    print(f"\n🎉 Batch processing complete!")
    print(f"📊 Results:")
    print(f"   ✅ Successful: {successful}")
    print(f"   ❌ Failed: {failed}")
    print(f"   📁 Total: {len(fasta_files)}")
    print(f"   📈 Success rate: {successful/len(fasta_files)*100:.1f}%")
    
    return successful > 0

def main():
    if len(sys.argv) == 1:
        # Default: show folder selection menu
        directory_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
        select_and_process_folder(directory_path, mode='interactive')
        
    elif len(sys.argv) == 2:
        arg = sys.argv[1]
        
        if arg == '--auto':
            # Show folder selection menu for automatic mode
            directory_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
            select_and_process_folder(directory_path, mode='auto')
            
        elif arg == '--interactive':
            # Show folder selection menu for interactive mode
            directory_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
            select_and_process_folder(directory_path, mode='interactive')
            
        elif arg == '--help':
            print_usage()
            
        elif os.path.isfile(arg) and arg.endswith('.fa'):
            # Single file mode (interactive)
            uploader = ImprovedUploader(arg)
            uploader.process_file_interactive()
            
        elif os.path.isdir(arg):
            # Show folder selection for specified directory
            select_and_process_folder(arg, mode='interactive')
            
        else:
            print(f"❌ Invalid argument: {arg}")
            print_usage()
            
    elif len(sys.argv) == 3:
        mode_arg = sys.argv[1]
        target_arg = sys.argv[2]
        
        if mode_arg == '--auto' and os.path.isdir(target_arg):
            # Show folder selection for specified directory in automatic mode
            select_and_process_folder(target_arg, mode='auto')
            
        elif mode_arg == '--interactive' and os.path.isdir(target_arg):
            # Show folder selection for specified directory in interactive mode
            select_and_process_folder(target_arg, mode='interactive')
            
        elif mode_arg == '--folder' and os.path.isdir(target_arg):
            # Process specific folder directly (interactive)
            process_directory_files(target_arg, mode='interactive')
            
        elif mode_arg == '--file' and os.path.isfile(target_arg) and target_arg.endswith('.fa'):
            # Single file interactive
            uploader = ImprovedUploader(target_arg)
            uploader.process_file_interactive()
            
        else:
            print(f"❌ Invalid arguments")
            print_usage()
            
    elif len(sys.argv) == 4:
        if sys.argv[1] == '--folder':
            # Process specific folder with mode
            mode = sys.argv[2]
            folder_path = sys.argv[3]
            
            if os.path.isdir(folder_path):
                if mode == 'auto':
                    process_directory_files(folder_path, mode='auto')
                elif mode == 'interactive':
                    process_directory_files(folder_path, mode='interactive')
                else:
                    print(f"❌ Invalid mode: {mode}")
            else:
                print(f"❌ Invalid folder: {folder_path}")
                
        elif sys.argv[1] == '--file':
            # Single file with specified mode
            mode = sys.argv[2]
            fasta_file = sys.argv[3]
            
            if not (os.path.isfile(fasta_file) and fasta_file.endswith('.fa')):
                print(f"❌ Invalid file: {fasta_file}")
                return
            
            uploader = ImprovedUploader(fasta_file)
            if mode == 'auto':
                uploader.process_file_automatic()
            elif mode == 'interactive':
                uploader.process_file_interactive()
            else:
                print(f"❌ Invalid mode: {mode}")
                print_usage()
        else:
            print_usage()
    else:
        print_usage()

def print_usage():
    print("Improved Interactive FASTA Uploader")
    print("="*50)
    print("Usage:")
    print("  python improved_uploader.py                           # Show folder selection menu (interactive)")
    print("  python improved_uploader.py --auto                    # Show folder selection menu (automatic)")
    print("  python improved_uploader.py --interactive             # Show folder selection menu (interactive)")
    print("  python improved_uploader.py <directory>               # Show folder selection for directory")
    print("  python improved_uploader.py --auto <directory>        # Show folder selection for directory (auto)")
    print("  python improved_uploader.py --folder <folder_path>    # Process specific folder (interactive)")
    print("  python improved_uploader.py --folder auto <folder>    # Process specific folder (automatic)")
    print("  python improved_uploader.py <file.fa>                 # Process single file (interactive)")
    print("  python improved_uploader.py --file auto <file.fa>     # Process single file (automatic)")
    print()
    print("Examples:")
    print("  python improved_uploader.py                           # Show folder menu")
    print("  python improved_uploader.py --auto                    # Show folder menu (auto mode)")
    print("  python improved_uploader.py --folder auto 1.16.99.1  # Process folder directly")
    print("  python improved_uploader.py protein1.fa               # Single file")
    print()
    print("Features:")
    print("  - Maximum 5 files per folder (sorted alphabetically)")
    print("  - Organized HTML results in matching folder structure")
    print("  - Batch command generation for multiple files")
    print("  - Preserves uploaded files on remote server")
    print()
    print("Directory Structure:")
    print("  Local Input:   fasta_files/1.16.99.1/protein.fa")
    print("  Remote Upload: /mnt/data2/supfam/zhijing/fasta_sequence/1.16.99.1/protein.fa")
    print("  Local Results: html_results/1.16.99.1/protein.html")
    print()
    print("Workflow:")
    print("  1. Create remote directory: /mnt/data2/supfam/zhijing/fasta_sequence/{folder_name}/")
    print("  2. Upload files to remote directory (files are preserved)")
    print("  3. Create symlinks in SUPERFAMILY directory")
    print("  4. Manual SUPERFAMILY analysis (batch commands provided)")
    print("  5. Download HTML results to organized folders")
    print("  6. Clean up only SUPERFAMILY result files (uploaded files kept)")

if __name__ == "__main__":
    main()