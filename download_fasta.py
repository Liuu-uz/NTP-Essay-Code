#!/usr/bin/env python3
"""
FASTA Files Downloader
Download all .fa files from remote server /mnt/data2/supfam/expasy/
to local directory /Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_sequence/
"""

import os
import sys
import subprocess

class FASTADownloader:
    def __init__(self):
        # Remote server configuration
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_dir = "/mnt/data2/supfam/expasy"
        
        # Local directory configuration
        self.local_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_sequence"
        
        # Create local directory if not exists
        self.setup_local_directory()
    
    def setup_local_directory(self):
        """Create local directory if it doesn't exist"""
        if not os.path.exists(self.local_dir):
            try:
                os.makedirs(self.local_dir)
                print(f"✅ Created local directory: {self.local_dir}")
            except Exception as e:
                print(f"❌ Failed to create directory {self.local_dir}: {e}")
                return False
        else:
            print(f"✅ Local directory exists: {self.local_dir}")
        return True
    
    def list_remote_fasta_files(self):
        """List all .fa files in remote directory"""
        print(f"🔍 Searching for .fa files in remote directory...")
        
        # Command to list .fa files on remote server - use basename to get only filenames
        list_cmd = f"ssh {self.remote_server} 'cd {self.remote_dir} && ls -1 *.fa 2>/dev/null'"
        
        try:
            result = subprocess.run(list_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"❌ Failed to access remote directory or no .fa files found")
                print(f"   Error: {result.stderr.strip()}")
                return []
            
            # Parse the ls output to get filenames
            fasta_files = []
            lines = result.stdout.strip().split('\n')
            
            for line in lines:
                filename = line.strip()
                if filename and filename.endswith('.fa'):
                    fasta_files.append(filename)
            
            if fasta_files:
                print(f"✅ Found {len(fasta_files)} .fa files:")
                for i, filename in enumerate(fasta_files, 1):
                    print(f"   {i:2d}. {filename}")
            else:
                print(f"⚠️  No .fa files found in {self.remote_dir}")
            
            return fasta_files
            
        except Exception as e:
            print(f"❌ Error listing remote files: {e}")
            return []
    
    def get_file_sizes_batch(self, filenames):
        """Get file sizes for all files in one SSH session to avoid multiple password prompts"""
        if not filenames:
            return {}
        
        # Build a single command to get all file sizes at once
        size_commands = []
        for filename in filenames:
            size_commands.append(f"ls -lh {self.remote_dir}/{filename} 2>/dev/null | awk '{{print \"{filename}:\" $5}}' || echo \"{filename}:unknown\"")
        
        # Combine all commands with semicolon
        batch_cmd = f"ssh {self.remote_server} '{'; '.join(size_commands)}'"
        
        sizes = {}
        try:
            result = subprocess.run(batch_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if ':' in line:
                        filename, size = line.split(':', 1)
                        sizes[filename] = size.strip()
            
            # Fill in unknown for any missing files
            for filename in filenames:
                if filename not in sizes:
                    sizes[filename] = "unknown"
                    
        except Exception as e:
            print(f"⚠️  Could not get file sizes: {e}")
            # Return unknown for all files
            for filename in filenames:
                sizes[filename] = "unknown"
        
        return sizes
    
    def download_single_file(self, filename, file_size="unknown"):
        """Download a single .fa file"""
        remote_path = f"{self.remote_server}:{self.remote_dir}/{filename}"
        local_path = os.path.join(self.local_dir, filename)
        
        # Check if file already exists locally
        if os.path.exists(local_path):
            local_size = os.path.getsize(local_path)
            print(f"   ⚠️  File already exists locally ({local_size} bytes)")
            
            choice = input(f"   Overwrite {filename}? (y/N): ").strip().lower()
            if choice not in ['y', 'yes']:
                print(f"   ⏩ Skipping {filename}")
                return True
        
        # Use provided file size to avoid extra SSH call
        print(f"   📥 Downloading {filename} (size: {file_size})...")
        
        # Download using scp
        scp_cmd = f"scp {remote_path} {local_path}"
        
        try:
            result = subprocess.run(scp_cmd, shell=True, capture_output=True)
            
            if result.returncode == 0:
                # Verify download
                if os.path.exists(local_path):
                    local_size = os.path.getsize(local_path)
                    print(f"   ✅ Downloaded successfully ({local_size} bytes)")
                    print(f"      Saved to: {local_path}")
                    return True
                else:
                    print(f"   ❌ Download failed - file not found locally")
                    return False
            else:
                print(f"   ❌ Download failed")
                if result.stderr:
                    print(f"      Error: {result.stderr.decode().strip()}")
                return False
                
        except Exception as e:
            print(f"   ❌ Download error: {e}")
            return False
    
    def download_all_files(self):
        """Download all .fa files"""
        print(f"🧬 FASTA Files Downloader")
        print("="*60)
        print(f"📡 Remote: {self.remote_server}:{self.remote_dir}")
        print(f"💾 Local:  {self.local_dir}")
        print("="*60)
        
        # List remote files
        fasta_files = self.list_remote_fasta_files()
        
        if not fasta_files:
            print(f"❌ No files to download")
            return False
        
        # Confirm download
        print(f"\n📋 Getting file sizes...")
        file_sizes = self.get_file_sizes_batch(fasta_files)
        
        print(f"\n📋 Files to download:")
        total_files = len(fasta_files)
        for i, filename in enumerate(fasta_files, 1):
            size = file_sizes.get(filename, "unknown")
            print(f"   {i:2d}. {filename:30s} ({size})")
        
        print(f"\n📥 Total files: {total_files}")
        confirm = input(f"Download all files? (y/N): ").strip().lower()
        
        if confirm not in ['y', 'yes']:
            print("👋 Download cancelled")
            return False
        
        # Download files
        print(f"\n🚀 Starting download...")
        print("-" * 60)
        
        successful = 0
        failed = 0
        
        for i, filename in enumerate(fasta_files, 1):
            print(f"\n[{i}/{len(fasta_files)}] {filename}")
            file_size = file_sizes.get(filename, "unknown")
            
            if self.download_single_file(filename, file_size):
                successful += 1
            else:
                failed += 1
        
        # Summary
        print(f"\n{'='*60}")
        print("📊 Download Summary")
        print(f"{'='*60}")
        print(f"✅ Successful: {successful}")
        print(f"❌ Failed:     {failed}")
        print(f"📁 Total:      {len(fasta_files)}")
        print(f"💾 Local dir:  {self.local_dir}")
        
        if successful > 0:
            print(f"\n🎉 {successful} files downloaded successfully!")
            
            # List downloaded files
            print(f"\n📋 Downloaded files:")
            try:
                local_files = [f for f in os.listdir(self.local_dir) if f.endswith('.fa')]
                for filename in sorted(local_files):
                    file_path = os.path.join(self.local_dir, filename)
                    file_size = os.path.getsize(file_path)
                    print(f"   📄 {filename} ({file_size} bytes)")
            except:
                pass
        
        print("="*60)
        return successful > 0
    
    def download_specific_files(self, file_patterns):
        """Download specific files based on patterns"""
        print(f"🧬 FASTA Files Downloader (Selective)")
        print("="*60)
        print(f"📡 Remote: {self.remote_server}:{self.remote_dir}")
        print(f"💾 Local:  {self.local_dir}")
        print(f"🔍 Patterns: {', '.join(file_patterns)}")
        print("="*60)
        
        # List all remote files
        all_files = self.list_remote_fasta_files()
        
        if not all_files:
            print(f"❌ No .fa files found on remote server")
            return False
        
        # Filter files based on patterns
        matching_files = []
        for pattern in file_patterns:
            for filename in all_files:
                if pattern.lower() in filename.lower() and filename not in matching_files:
                    matching_files.append(filename)
        
        if not matching_files:
            print(f"❌ No files match the given patterns")
            return False
        
        print(f"\n✅ Found {len(matching_files)} matching files:")
        file_sizes = self.get_file_sizes_batch(matching_files)
        
        for i, filename in enumerate(matching_files, 1):
            size = file_sizes.get(filename, "unknown")
            print(f"   {i:2d}. {filename} ({size})")
        
        # Download matching files
        confirm = input(f"\nDownload these {len(matching_files)} files? (y/N): ").strip().lower()
        
        if confirm not in ['y', 'yes']:
            print("👋 Download cancelled")
            return False
        
        print(f"\n🚀 Starting selective download...")
        print("-" * 60)
        
        successful = 0
        for i, filename in enumerate(matching_files, 1):
            print(f"\n[{i}/{len(matching_files)}] {filename}")
            file_size = file_sizes.get(filename, "unknown")
            if self.download_single_file(filename, file_size):
                successful += 1
        
        print(f"\n🎉 {successful}/{len(matching_files)} files downloaded!")
        return successful > 0

def main():
    downloader = FASTADownloader()
    
    if len(sys.argv) == 1:
        # Download all .fa files
        downloader.download_all_files()
        
    elif len(sys.argv) >= 2:
        if sys.argv[1] == "--help":
            print("Usage:")
            print("  python fasta_downloader.py                    # Download all .fa files")
            print("  python fasta_downloader.py pattern1 pattern2  # Download files matching patterns")
            print()
            print("Examples:")
            print("  python fasta_downloader.py                    # Download all files")
            print("  python fasta_downloader.py protein kinase     # Download files containing 'protein' or 'kinase'")
            print("  python fasta_downloader.py 1abc 2def          # Download files containing '1abc' or '2def'")
            
        else:
            # Download specific files based on patterns
            patterns = sys.argv[1:]
            downloader.download_specific_files(patterns)
    
    else:
        print("Usage: python fasta_downloader.py [pattern1] [pattern2] ...")

if __name__ == "__main__":
    main()