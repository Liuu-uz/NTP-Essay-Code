#!/usr/bin/env python3
"""
Script to filter CSV files based on protein IDs from Excel file.
Keeps only CSV rows where the filename (without extension and suffix) matches 
protein IDs in the first column of the Excel file.
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

def extract_protein_id_from_filename(filename):
    """
    Extract protein ID from CSV filename.
    Example: 'a0a1A_vs_3amtA.txt' -> 'a0a1A'
    """
    basename = os.path.basename(filename)
    # Remove extension
    name_without_ext = os.path.splitext(basename)[0]
    # Extract the part before '_vs_'
    if '_vs_' in name_without_ext:
        protein_id = name_without_ext.split('_vs_')[0]
        return protein_id
    return name_without_ext

def load_excel_protein_ids(excel_path):
    """
    Load protein IDs from the first column of Excel file.
    """
    try:
        df = pd.read_excel(excel_path)
        # Get first column values, skip header, remove NaN values
        protein_ids = df.iloc[1:, 0].dropna().tolist()
        # Convert to set for faster lookup
        return set(str(pid).strip() for pid in protein_ids)
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        return set()

def process_csv_files(csv_directory, excel_path, output_directory=None):
    """
    Process all CSV files in directory and filter based on Excel protein IDs.
    """
    # Load protein IDs from Excel
    valid_protein_ids = load_excel_protein_ids(excel_path)
    print(f"Loaded {len(valid_protein_ids)} protein IDs from Excel file")
    print(f"Sample protein IDs: {list(valid_protein_ids)[:10]}")
    
    # Create output directory if specified
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)
    
    # Find all CSV files
    csv_pattern = os.path.join(csv_directory, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    # Also look for .txt files as your example shows .txt extension
    txt_pattern = os.path.join(csv_directory, "*.txt")
    txt_files = glob.glob(txt_pattern)
    
    all_files = csv_files + txt_files
    print(f"Found {len(all_files)} files to process")
    
    processed_count = 0
    filtered_count = 0
    
    for file_path in all_files:
        try:
            # Extract protein ID from filename
            protein_id = extract_protein_id_from_filename(file_path)
            
            # Check if this protein ID is in our valid set
            if protein_id in valid_protein_ids:
                # Read the CSV file
                df = pd.read_csv(file_path)
                
                # Determine output path
                if output_directory:
                    output_path = os.path.join(output_directory, os.path.basename(file_path))
                else:
                    # Overwrite original file
                    output_path = file_path
                
                # Save filtered file
                df.to_csv(output_path, index=False)
                print(f"Kept: {os.path.basename(file_path)} (protein ID: {protein_id})")
                filtered_count += 1
            else:
                print(f"Removed: {os.path.basename(file_path)} (protein ID: {protein_id} not in Excel)")
                if not output_directory:
                    # Delete the file if not using output directory
                    os.remove(file_path)
            
            processed_count += 1
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    print(f"\nProcessing complete!")
    print(f"Total files processed: {processed_count}")
    print(f"Files kept: {filtered_count}")
    print(f"Files removed: {processed_count - filtered_count}")

def main():
    # Configuration - Update these paths
    csv_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/batch_results"
    excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx"
    
    # Optional: specify output directory to keep original files intact
    # output_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/filtered_results"
    output_directory = None  # Set to None to overwrite original files
    
    # Verify paths exist
    if not os.path.exists(csv_directory):
        print(f"Error: CSV directory not found: {csv_directory}")
        return
    
    if not os.path.exists(excel_path):
        print(f"Error: Excel file not found: {excel_path}")
        return
    
    print(f"CSV directory: {csv_directory}")
    print(f"Excel file: {excel_path}")
    if output_directory:
        print(f"Output directory: {output_directory}")
    else:
        print("Will overwrite original files")
    
    # Process files
    process_csv_files(csv_directory, excel_path, output_directory)

if __name__ == "__main__":
    main()