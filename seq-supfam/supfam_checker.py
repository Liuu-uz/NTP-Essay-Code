#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch HTML and NTP Superfamily Comparator
Process organized HTML folders and generate comprehensive Excel report
"""

import os
import sys
import pandas as pd
from bs4 import BeautifulSoup
from datetime import datetime

class BatchNTPChecker:
    def __init__(self, excel_path=None):
        self.excel_path = excel_path
        self.ntp_sf_dict = {}
        self.sf_names = []
        self.ntp_data = {}
        self.find_excel_file()
        self.load_ntp_data()
        
        # Statistics for final report
        self.total_files = 0
        self.total_folders = 0
        self.ntp_matches = 0
        self.results_data = []  # Store all results for Excel export
    
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
            
            # Store complete row information
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
        """Parse SUPERFAMILY HTML result file (simplified version)"""
        if not os.path.exists(html_path):
            return []
        
        try:
            with open(html_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Fix HTML format issues
            content = content.replace('<table border="1" />', '<table border="1">')
            if '</table>' not in content and '<table' in content:
                content = content.replace('</body>', '</table></body>')
            
            soup = BeautifulSoup(content, 'html.parser')
            found_superfamilies = []
            
            # Find tables and process rows
            tables = soup.find_all('table')
            
            for table in tables:
                rows = table.find_all('tr')
                
                for row_idx, row in enumerate(rows):
                    if row_idx == 0:  # Skip header
                        continue
                    
                    cells = row.find_all(['td', 'th'])
                    if len(cells) >= 4:
                        try:
                            seq_id = cells[0].get_text(strip=True)
                            scop_superfamily_cell = cells[3]
                            
                            # Extract superfamily name
                            sf_link = scop_superfamily_cell.find('a')
                            if sf_link:
                                sf_name = sf_link.get_text(strip=True)
                            else:
                                sf_name = scop_superfamily_cell.get_text(strip=True)
                            
                            if sf_name and len(sf_name) > 2 and seq_id:
                                found_superfamilies.append({
                                    'name': sf_name,
                                    'seq_id': seq_id
                                })
                                
                        except Exception:
                            continue
            
            return found_superfamilies
            
        except Exception:
            return []
    
    def check_ntp_matches(self, found_superfamilies):
        """Check if found superfamilies match with NTP-processing"""
        matches = []
        
        for sf in found_superfamilies:
            sf_name = sf['name']
            seq_id = sf.get('seq_id', 'unknown')
            
            matched_entry = None
            match_type = None
            
            # 1. Superfamily name exact match
            if sf_name and sf_name in self.ntp_data:
                matched_entry = self.ntp_data[sf_name]
                match_type = 'SF_name_exact'
            
            # 2. Fuzzy matching
            elif sf_name:
                for ntp_name in self.sf_names:
                    if len(ntp_name) > 10:
                        sf_lower = sf_name.lower()
                        ntp_lower = ntp_name.lower()
                        
                        if (ntp_lower == sf_lower or 
                            (len(ntp_lower) > 10 and ntp_lower in sf_lower) or
                            (len(sf_lower) > 10 and sf_lower in ntp_lower)):
                            
                            matched_entry = self.ntp_data[ntp_name]
                            match_type = 'fuzzy'
                            break
            
            if matched_entry and match_type:
                matches.append({
                    'found_name': sf_name,
                    'found_seq_id': seq_id,
                    'matched_entry': matched_entry,
                    'match_type': match_type
                })
        
        return matches
    
    def process_single_html(self, html_path, folder_name):
        """Process single HTML file and store results"""
        html_filename = os.path.basename(html_path)
        sequence_id = os.path.splitext(html_filename)[0]

        # Parse HTML file
        found_superfamilies = self.parse_html_file(html_path)

        if not found_superfamilies:
            self.results_data.append({
                'Folder': folder_name,
                'HTML_File': html_filename,
                'Sequence_ID': sequence_id,
                'Found_Superfamily': 'N/A',
                'NTP_Match': 'No',
                'NTP_Superfamily': 'N/A',
                'PDB_ID': 'N/A',
                'Result_Category': 'DALI'
            })
            return False

        # Check for NTP matches
        matches = self.check_ntp_matches(found_superfamilies)

        if matches:
            match = matches[0]
            entry = match['matched_entry']
            self.results_data.append({
                'Folder': folder_name,
                'HTML_File': html_filename,
                'Sequence_ID': sequence_id,
                'Found_Superfamily': match['found_name'],
                'NTP_Match': 'Yes',
                'NTP_Superfamily': entry['sf_name'] if entry['sf_name'] else 'N/A',
                'PDB_ID': entry['representative'] if entry['representative'] else 'N/A',
                'Result_Category': 'NTP processing SF found'
            })
            self.ntp_matches += 1
            return True
        else:
            sf = found_superfamilies[0]
            self.results_data.append({
                'Folder': folder_name,
                'HTML_File': html_filename,
                'Sequence_ID': sequence_id,
                'Found_Superfamily': sf['name'],
                'NTP_Match': 'No',
                'NTP_Superfamily': 'N/A',
                'PDB_ID': 'N/A',
                'Result_Category': 'NO NTP processing SF found'
            })
            return False
    
    def process_html_folders(self, base_directory):
        """Process all folders containing HTML files"""
        if not os.path.exists(base_directory):
            print(f"❌ Directory does not exist: {base_directory}")
            return False
        
        print(f"🧬 Batch NTP Processing Analysis")
        print("="*60)
        print(f"📂 Base directory: {base_directory}")
        
        # Find all folders containing HTML files
        folders_with_html = []
        
        for item in os.listdir(base_directory):
            item_path = os.path.join(base_directory, item)
            if os.path.isdir(item_path):
                html_files = [f for f in os.listdir(item_path) if f.endswith('.html')]
                if html_files:
                    folders_with_html.append({
                        'name': item,
                        'path': item_path,
                        'html_count': len(html_files)
                    })
        
        if not folders_with_html:
            print("❌ No folders containing HTML files found")
            return False
        
        self.total_folders = len(folders_with_html)
        
        print(f"📁 Found {len(folders_with_html)} folders with HTML files:")
        for folder in folders_with_html:
            print(f"   📂 {folder['name']:<20} ({folder['html_count']} HTML files)")
        
        print(f"\n🔄 Processing all folders...")
        print("="*60)
        
        # Process each folder
        for i, folder in enumerate(folders_with_html, 1):
            print(f"\n[{i}/{len(folders_with_html)}] Processing folder: {folder['name']}")
            print("-" * 50)
            
            html_files = []
            for file in os.listdir(folder['path']):
                if file.endswith('.html'):
                    html_files.append(os.path.join(folder['path'], file))
            
            folder_ntp_count = 0
            for j, html_file in enumerate(html_files, 1):
                html_name = os.path.basename(html_file)
                print(f"   [{j}/{len(html_files)}] {html_name}... ", end="")
                
                result = self.process_single_html(html_file, folder['name'])
                if result:
                    folder_ntp_count += 1
                    print("✅ NTP match found")
                else:
                    print("❌ No NTP match")
                
                self.total_files += 1
            
            print(f"   📊 Folder summary: {folder_ntp_count}/{len(html_files)} files with NTP matches")
        
        return True
    
    def generate_excel_report(self, output_path=None):
        """Generate comprehensive Excel report"""
        if not self.results_data:
            print("❌ No data to export")
            return False
        
        # Generate output filename if not provided
        if not output_path:
            output_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam"
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_filename = f"NTP_Analysis_Report_{timestamp}.xlsx"
            
            # Ensure output directory exists
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    print(f"📁 Created output directory: {output_dir}")
                except Exception as e:
                    print(f"⚠️  Could not create directory {output_dir}, using current directory: {e}")
                    output_dir = "."
            
            output_path = os.path.join(output_dir, output_filename)
        
        try:
            # Create DataFrame
            df = pd.DataFrame(self.results_data)
            
            # Create Excel writer with multiple sheets
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                # Main results sheet
                df.to_excel(writer, sheet_name='Detailed_Results', index=False)
                
                # Summary by folder
                folder_summary = df.groupby('Folder').agg({
                    'HTML_File': 'count',
                    'NTP_Match': lambda x: (x == 'Yes').sum()
                }).rename(columns={
                    'HTML_File': 'Total_Files',
                    'NTP_Match': 'NTP_Matches'
                })
                folder_summary['Match_Rate'] = (folder_summary['NTP_Matches'] / folder_summary['Total_Files'] * 100).round(2)
                folder_summary.to_excel(writer, sheet_name='Folder_Summary')
                
                # NTP families found
                ntp_families = df[df['NTP_Match'] == 'Yes']['NTP_Superfamily'].value_counts()
                ntp_families.to_excel(writer, sheet_name='NTP_Families_Found')
                
                # Overall statistics
                stats_data = {
                    'Metric': [
                        'Total Folders',
                        'Total HTML Files',
                        'Files with NTP Matches',
                        'Files without NTP Matches',
                        'Overall NTP Match Rate (%)',
                        'Unique NTP Families Found',
                        'Analysis Date',
                        'Output Directory'
                    ],
                    'Value': [
                        self.total_folders,
                        self.total_files,
                        len(df[df['NTP_Match'] == 'Yes']),
                        len(df[df['NTP_Match'] == 'No']),
                        f"{len(df[df['NTP_Match'] == 'Yes'])/len(df)*100:.2f}%",
                        len([x for x in df[df['NTP_Match'] == 'Yes']['NTP_Superfamily'].unique() if x != 'N/A']),
                        datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        os.path.abspath(output_path)
                    ]
                }
                stats_df = pd.DataFrame(stats_data)
                stats_df.to_excel(writer, sheet_name='Statistics', index=False)
            
            print(f"✅ Excel report generated: {output_path}")
            print(f"📊 Report contains {len(df)} entries across {self.total_folders} folders")
            print(f"📁 Saved to: {os.path.dirname(output_path)}")
            return output_path
            
        except Exception as e:
            print(f"❌ Failed to generate Excel report: {e}")
            return False
    
    def print_summary(self):
        """Print processing summary"""
        print(f"\n{'='*60}")
        print("📊 Final Processing Summary")
        print(f"{'='*60}")
        print(f"📂 Folders processed: {self.total_folders}")
        print(f"📄 HTML files processed: {self.total_files}")
        print(f"✅ Files with NTP matches: {len([r for r in self.results_data if r['NTP_Match'] == 'Yes'])}")
        print(f"❌ Files without NTP matches: {len([r for r in self.results_data if r['NTP_Match'] == 'No'])}")
        
        if self.total_files > 0:
            ntp_rate = len([r for r in self.results_data if r['NTP_Match'] == 'Yes']) / self.total_files * 100
            print(f"📈 Overall NTP match rate: {ntp_rate:.2f}%")
        
        # Show unique NTP families found
        ntp_families = set([r['NTP_Superfamily'] for r in self.results_data if r['NTP_Match'] == 'Yes' and r['NTP_Superfamily'] != 'N/A'])
        if ntp_families:
            print(f"\n🧬 Unique NTP families found ({len(ntp_families)}):")
            for family in sorted(ntp_families):
                print(f"   - {family}")
        
        print("="*60)

def main():
    # Create checker instance
    checker = BatchNTPChecker()
    
    # Check if Excel file was successfully loaded
    if not checker.excel_path:
        print("\n💡 If Excel file is in another location, you can specify the path:")
        print("python batch_ntp_checker.py --excel /path/to/excel/file.xlsx html_results/")
        return
    
    if len(sys.argv) == 1:
        # Default: check html_results directory
        html_results_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/html_results"
        if os.path.exists(html_results_dir):
            print(f"📂 Using default directory: {html_results_dir}")
            if checker.process_html_folders(html_results_dir):
                checker.print_summary()
                output_path = checker.generate_excel_report()
                if output_path:
                    print(f"📄 Excel report saved: {output_path}")
        else:
            print(f"❌ Default html_results directory not found: {html_results_dir}")
            print("💡 Usage: python batch_ntp_checker.py [directory_path]")
        
    elif len(sys.argv) == 2:
        if sys.argv[1] == "--help":
            print("Batch HTML NTP Checker")
            print("="*50)
            print("Usage:")
            print("  python batch_ntp_checker.py                    # Use default html_results directory")
            print("  python batch_ntp_checker.py <directory>        # Process specified directory")
            print("  python batch_ntp_checker.py --excel <path> <dir>  # Specify Excel file")
            print()
            print("Output:")
            print("  - Comprehensive Excel report with multiple sheets")
            print("  - Detailed results, folder summaries, and statistics")
            print("  - NTP family matches with PDB IDs and reaction types")
            
        elif os.path.isdir(sys.argv[1]):
            directory = sys.argv[1]
            print(f"📂 Processing directory: {directory}")
            if checker.process_html_folders(directory):
                checker.print_summary()
                output_path = checker.generate_excel_report()
                if output_path:
                    print(f"📄 Excel report saved: {output_path}")
        else:
            print(f"❌ Directory does not exist: {sys.argv[1]}")
    
    elif len(sys.argv) == 4 and sys.argv[1] == "--excel":
        excel_path = sys.argv[2]
        directory = sys.argv[3]
        
        if not os.path.exists(excel_path):
            print(f"❌ Excel file does not exist: {excel_path}")
            return
        
        if not os.path.isdir(directory):
            print(f"❌ Directory does not exist: {directory}")
            return
        
        checker = BatchNTPChecker(excel_path)
        if checker.process_html_folders(directory):
            checker.print_summary()
            output_path = checker.generate_excel_report()
            if output_path:
                print(f"📄 Excel report saved: {output_path}")
    
    else:
        print("Usage:")
        print("  python batch_ntp_checker.py                    # Use default html_results directory")
        print("  python batch_ntp_checker.py <directory>        # Process specified directory")
        print("  python batch_ntp_checker.py --excel <excel_path> <directory>  # Specify Excel file")

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