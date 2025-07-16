#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HTML and NTP Superfamily Comparator
Compare SUPERFAMILY generated HTML files with NTP-processing superfamilies in Excel table
"""

import os
import sys
import pandas as pd
from bs4 import BeautifulSoup

class NTPChecker:
    def __init__(self, excel_path=None):
        self.excel_path = excel_path
        self.ntp_sf_dict = {}
        self.sf_names = []
        self.find_excel_file()
        self.load_ntp_data()
    
    def find_excel_file(self):
        """Find Excel file"""
        if self.excel_path and os.path.exists(self.excel_path):
            return
        
        # Possible Excel file locations
        possible_paths = [
            "SUPFAM based NTP processing.xlsx",  # Current directory
            "../SUPFAM based NTP processing.xlsx",  # Parent directory
            "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx",  # Specified directory
            os.path.expanduser("~/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx"),  # User directory
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                self.excel_path = path
                print(f"✅ Found Excel file: {path}")
                return
        
        print("❌ Cannot find Excel file, please check the following locations:")
        for path in possible_paths:
            print(f"   {path}")
        self.excel_path = None
    
    def load_ntp_data(self):
        """Load NTP-processing superfamily data"""
        if not self.excel_path or not os.path.exists(self.excel_path):
            print(f"❌ Excel file does not exist or not found")
            print("💡 Please ensure the file exists in one of the following locations:")
            print("   - Current directory: SUPFAM based NTP processing.xlsx")
            print("   - Specified directory: /Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/")
            return False
        
        try:
            df = pd.read_excel(self.excel_path)
            
            # Store complete row information, including Representative (PDB ID)
            self.ntp_data = {}  # Store complete information
            
            for i, row in df.iterrows():
                sf_id = row['SF ID'] if pd.notna(row['SF ID']) else None
                sf_name = row['SF'] if pd.notna(row['SF']) else None
                representative = row['Representative'] if pd.notna(row['Representative']) else None
                reaction = row['Reaction'] if pd.notna(row['Reaction']) else None
                coordination = row['Coordination'] if pd.notna(row['Coordination']) else None
                
                # Create entry information
                entry_info = {
                    'sf_name': sf_name,
                    'sf_id': sf_id,
                    'representative': representative,
                    'reaction': reaction,
                    'coordination': coordination
                }
                
                if sf_id and sf_id != 'N/A':
                    self.ntp_sf_dict[sf_id] = sf_name
                    self.ntp_data[sf_id] = entry_info
                
                if sf_name:
                    self.sf_names.append(sf_name)
                    if sf_id and sf_id != 'N/A':
                        self.ntp_sf_dict[sf_name] = sf_id
                        self.ntp_data[sf_name] = entry_info
                    else:
                        self.ntp_sf_dict[sf_name] = 'N/A'
                        self.ntp_data[sf_name] = entry_info
            
            sf_with_id = len([k for k in self.ntp_sf_dict.keys() if '.' in str(k)])
            print(f"✅ Loaded {len(df)} NTP-processing superfamilies")
            print(f"   {sf_with_id} with SF ID, {len(df) - sf_with_id} with name only")
            return True
            
        except Exception as e:
            print(f"❌ Failed to read Excel file: {e}")
            return False
    
    def parse_html_file(self, html_path):
        """Parse SUPERFAMILY HTML result file"""
        if not os.path.exists(html_path):
            print(f"❌ HTML file does not exist: {html_path}")
            return []
        
        try:
            # Read file content
            with open(html_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            print(f"   Debug: File size: {len(content)} characters")
            
            # Fix HTML format issues: change self-closing table tags to normal format
            # Original: <table border="1" />
            # Fixed: <table border="1">
            content = content.replace('<table border="1" />', '<table border="1">')
            
            # Ensure table tags are properly closed
            if '</table>' not in content and '<table' in content:
                # Add </table> before body end
                content = content.replace('</body>', '</table></body>')
            
            print(f"   Debug: Fixed HTML snippet:")
            # Find table section and display
            table_start = content.find('<table')
            table_end = content.find('</table>') + 8
            if table_start >= 0 and table_end > table_start:
                table_html = content[table_start:table_end]
                print(f"   {table_html[:300]}...")
            
            soup = BeautifulSoup(content, 'html.parser')
            found_superfamilies = []
            
            # Find tables
            tables = soup.find_all('table')
            print(f"   Debug: Found {len(tables)} tables after fixing")
            
            for table_idx, table in enumerate(tables):
                print(f"   Debug: Processing table {table_idx + 1}")
                
                rows = table.find_all('tr')
                print(f"   Debug: Found {len(rows)} tr tags")
                
                # Process each row
                for row_idx, row in enumerate(rows):
                    cells = row.find_all(['td', 'th'])
                    print(f"   Debug: Row {row_idx + 1} has {len(cells)} cells")
                    
                    # Display row content
                    if len(cells) > 0:
                        print(f"   Debug: Row {row_idx + 1} content:")
                        for i, cell in enumerate(cells):
                            cell_text = cell.get_text(strip=True)
                            print(f"     Column{i+1}: '{cell_text}'")
                            # Display links
                            link = cell.find('a')
                            if link:
                                link_text = link.get_text(strip=True)
                                print(f"       Link: '{link_text}'")
                    
                    # Skip header row
                    if row_idx == 0:
                        print("   Debug: This is header row, skipping")
                        continue
                    
                    # Process data rows
                    if len(cells) >= 4:
                        try:
                            seq_id = cells[0].get_text(strip=True)
                            
                            # Extract SCOP superfamily from column 4
                            scop_superfamily_cell = cells[3]
                            sf_link = scop_superfamily_cell.find('a')
                            
                            if sf_link:
                                sf_name = sf_link.get_text(strip=True)
                                print(f"   Debug: ✅ Extracted superfamily from column 4 link: '{sf_name}'")
                            else:
                                sf_name = scop_superfamily_cell.get_text(strip=True)
                                print(f"   Debug: Extracted superfamily from column 4 text: '{sf_name}'")
                            
                            # Determine SCOP ID
                            scop_id = 'unknown'
                            if sf_name == 'HIT-like':
                                scop_id = 'd.13.1'  # Known SCOP ID
                            
                            if sf_name and len(sf_name) > 2 and seq_id:
                                print(f"   Debug: ✅ Adding superfamily: name='{sf_name}', SCOP_ID='{scop_id}', sequence='{seq_id}'")
                                found_superfamilies.append({
                                    'name': sf_name,
                                    'scop_id': scop_id,
                                    'seq_id': seq_id
                                })
                            else:
                                print(f"   Debug: ❌ Skipping invalid entry")
                                
                        except Exception as e:
                            print(f"   Debug: Error processing row: {e}")
                            continue
            
            print(f"   Debug: Finally found {len(found_superfamilies)} superfamilies")
            return found_superfamilies
            
        except Exception as e:
            print(f"❌ Failed to parse HTML file: {e}")
            import traceback
            traceback.print_exc()
            return []
    
    def check_ntp_matches(self, found_superfamilies):
        """Check if found superfamilies match with NTP-processing"""
        if not found_superfamilies:
            return []
        
        matches = []
        
        for sf in found_superfamilies:
            sf_name = sf['name']
            scop_id = sf['scop_id']
            seq_id = sf.get('seq_id', 'unknown')
            
            matched_entry = None
            match_type = None
            
            # 1. SCOP ID exact match (primary matching method)
            if scop_id and scop_id in self.ntp_data:
                matched_entry = self.ntp_data[scop_id]
                match_type = 'SCOP_ID_exact'
            
            # 2. Superfamily name exact match
            elif sf_name and sf_name in self.ntp_data:
                matched_entry = self.ntp_data[sf_name]
                match_type = 'SF_name_exact'
            
            # 3. Look up corresponding superfamily name through SCOP ID in dictionary
            elif scop_id:
                # Check if Excel has superfamily corresponding to this SCOP ID
                for key, entry in self.ntp_data.items():
                    if entry['sf_id'] == scop_id:
                        matched_entry = entry
                        match_type = 'SCOP_ID_lookup'
                        break
            
            # 4. Smart fuzzy matching (only as last resort)
            if not matched_entry and sf_name:
                for ntp_name in self.sf_names:
                    if len(ntp_name) > 15:  # Only match longer names
                        sf_lower = sf_name.lower()
                        ntp_lower = ntp_name.lower()
                        
                        # Stricter matching conditions
                        if (ntp_lower == sf_lower or  # Exact match
                            (len(ntp_lower) > 10 and ntp_lower in sf_lower) or  # Long substring match
                            (len(sf_lower) > 10 and sf_lower in ntp_lower)):
                            
                            matched_entry = self.ntp_data[ntp_name]
                            match_type = 'fuzzy'
                            break
            
            if matched_entry and match_type:
                matches.append({
                    'found_name': sf_name,
                    'found_scop': scop_id,
                    'found_seq_id': seq_id,
                    'matched_entry': matched_entry,
                    'match_type': match_type
                })
        
        return matches
    
    def check_single_html(self, html_path):
        """Check single HTML file"""
        print(f"🔍 Checking file: {os.path.basename(html_path)}")
        print("-" * 50)
        
        # Parse HTML file
        found_superfamilies = self.parse_html_file(html_path)
        
        if not found_superfamilies:
            print("⚠️  No superfamily information found in HTML")
            print("❌ NO NTP-processing family found")
            return False
        
        print(f"📊 Found {len(found_superfamilies)} superfamilies in HTML:")
        for i, sf in enumerate(found_superfamilies, 1):
            print(f"   {i}. {sf['name']}")
            if sf['scop_id'] and sf['scop_id'] != 'unknown':
                print(f"      SCOP ID: {sf['scop_id']}")
            if sf.get('seq_id') and sf['seq_id'] != 'unknown':
                print(f"      Sequence ID: {sf['seq_id']}")
        
        # Check NTP matches
        matches = self.check_ntp_matches(found_superfamilies)
        
        print(f"\n🧬 NTP-processing match results:")
        print("=" * 50)
        
        if matches:
            print(f"✅ Found {len(matches)} NTP-processing related superfamilies!")
            for i, match in enumerate(matches, 1):
                entry = match['matched_entry']
                print(f"\n   {i}. Found superfamily: {match['found_name']}")
                if match['found_scop'] and match['found_scop'] != 'unknown':
                    print(f"      Found SCOP ID: {match['found_scop']}")
                
                print(f"      ✅ Matched NTP superfamily: {entry['sf_name']}")
                if entry['sf_id'] and entry['sf_id'] != 'N/A':
                    print(f"      ✅ Matched SCOP ID: {entry['sf_id']}")
                if entry['representative']:
                    print(f"      🧬 PDB ID: {entry['representative']}")
                if entry['reaction']:
                    print(f"      ⚗️  Reaction type: {entry['reaction']}")
                if entry['coordination']:
                    print(f"      🔗 Coordination: {entry['coordination']}")
                print(f"      📊 Match type: {match['match_type']}")
            
            print(f"\n🔬 Suggestion: Please perform AlphaFold prediction to find active sites!")
            return True
        else:
            print("❌ NO NTP-processing family found")
            return False
    
    def check_directory(self, directory_path):
        """Check all HTML files in directory"""
        if not os.path.exists(directory_path):
            print(f"❌ Directory does not exist: {directory_path}")
            return
        
        # Find all HTML files
        html_files = []
        for file in os.listdir(directory_path):
            if file.endswith('.html') or file.endswith('.htm'):
                html_files.append(os.path.join(directory_path, file))
        
        if not html_files:
            print(f"❌ No HTML files found in directory {directory_path}")
            return
        
        print(f"🧬 Batch NTP-processing check")
        print("=" * 60)
        print(f"📁 Found {len(html_files)} HTML files in {directory_path}")
        print("=" * 60)
        
        # Statistics results
        results = {"total": len(html_files), "ntp_found": 0, "no_result": 0}
        
        for i, html_file in enumerate(html_files, 1):
            print(f"\n[{i}/{len(html_files)}] " + "="*60)
            
            result = self.check_single_html(html_file)
            if result:
                results["ntp_found"] += 1
            elif not result:
                results["no_result"] += 1
        
        # Generate summary report
        print(f"\n{'='*60}")
        print("📊 Batch check summary report")
        print(f"{'='*60}")
        print(f"Total files: {results['total']}")
        print(f"NTP-processing found: {results['ntp_found']}")
        print(f"NTP-processing not found: {results['total'] - results['ntp_found']}")
        print(f"NTP match rate: {results['ntp_found']/results['total']*100:.1f}%")
        
        if results['ntp_found'] > 0:
            print(f"\n🔬 Suggestion: {results['ntp_found']} proteins need AlphaFold prediction!")
        
        print("=" * 60)

def main():
    # Create checker instance
    checker = NTPChecker()
    
    # Check if Excel file was successfully loaded
    if not checker.excel_path:
        print("\n💡 If Excel file is in another location, you can specify the path:")
        print("python html_ntp_checker.py --excel /path/to/excel/file.xlsx file.html")
        return
    
    if len(sys.argv) == 1:
        # No parameters, check current directory or html_results directory
        possible_dirs = ["html_results", ".", "downloaded_results"]
        
        for dir_path in possible_dirs:
            if os.path.exists(dir_path):
                html_files = [f for f in os.listdir(dir_path) if f.endswith('.html')]
                if html_files:
                    checker.check_directory(dir_path)
                    return
        
        print("❌ No HTML files found")
        print("💡 Usage:")
        print("   python html_ntp_checker.py                    # Auto-find HTML files")
        print("   python html_ntp_checker.py file.html          # Check single file")
        print("   python html_ntp_checker.py html_results/      # Check directory")
        
    elif len(sys.argv) >= 2:
        # Process command line arguments
        if sys.argv[1] == "--excel" and len(sys.argv) >= 4:
            # Specify Excel file path
            excel_path = sys.argv[2]
            target = sys.argv[3]
            checker = NTPChecker(excel_path)
            
            if not checker.excel_path:
                print(f"❌ Specified Excel file does not exist: {excel_path}")
                return
        else:
            target = sys.argv[1]
        
        if os.path.isfile(target):
            # Check single file
            print(f"🧬 Single file NTP-processing check")
            print("=" * 50)
            checker.check_single_html(target)
            
        elif os.path.isdir(target):
            # Check directory
            checker.check_directory(target)
            
        else:
            print(f"❌ File or directory does not exist: {target}")
    else:
        print("Usage:")
        print("  python html_ntp_checker.py                              # Auto-find HTML files") 
        print("  python html_ntp_checker.py file.html                    # Check single file")
        print("  python html_ntp_checker.py html_results/                # Check directory")
        print("  python html_ntp_checker.py --excel /path/to/excel.xlsx file.html  # Specify Excel path")

if __name__ == "__main__":
    # Check dependencies
    try:
        import pandas as pd
        from bs4 import BeautifulSoup
    except ImportError as e:
        print("❌ Missing required Python packages, please install:")
        print("pip install pandas openpyxl beautifulsoup4 lxml")
        print(f"Error details: {e}")
        sys.exit(1)
    
    main()