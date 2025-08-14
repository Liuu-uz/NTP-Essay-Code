import os
import csv
import glob

def extract_first_data_row_from_csv(csv_file_path):
    """Extract the first data row from CSV file (skip header row)"""
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as file:
            reader = csv.reader(file)
            
            # Skip header row
            try:
                header = next(reader)
                print(f"  Header row: {header}")
            except StopIteration:
                print(f"  Empty file: {csv_file_path}")
                return None, None
            
            # Read first data row
            try:
                first_data_row = next(reader)
                return header, first_data_row
            except StopIteration:
                print(f"  Header only, no data: {csv_file_path}")
                return header, None
                
    except Exception as e:
        print(f"Error reading file {csv_file_path}: {e}")
        return None, None

def extract_first_rows_batch(input_directory, output_file):
    """Batch extract first data rows from all CSV files and merge"""
    
    # Find all CSV files
    csv_pattern = os.path.join(input_directory, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"No CSV files found in directory {input_directory}")
        return
    
    print(f"Found {len(csv_files)} CSV files")
    
    all_headers = set()
    successful_extractions = []
    failed_files = []
    
    # First collect all possible column headers
    for csv_file in csv_files:
        header, _ = extract_first_data_row_from_csv(csv_file)
        if header:
            all_headers.update(header)
    
    # Create unified column header list
    unified_headers = ['Source_File'] + sorted(list(all_headers))
    
    # Process each file
    for csv_file in csv_files:
        print(f"\nProcessing file: {os.path.basename(csv_file)}")
        
        header, first_data_row = extract_first_data_row_from_csv(csv_file)
        
        if first_data_row:
            # Create dictionary to store this file's data
            row_data = {'Source_File': os.path.basename(csv_file)}
            
            # Map data according to column headers
            for i, value in enumerate(first_data_row):
                if i < len(header):
                    row_data[header[i]] = value.strip() if value else ''
            
            successful_extractions.append(row_data)
            print(f"  Successfully extracted first data row: {first_data_row}")
        else:
            failed_files.append(os.path.basename(csv_file))
    
    # Save results
    if successful_extractions:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=unified_headers)
            
            # Write header row
            writer.writeheader()
            
            # Write all data
            for row_data in successful_extractions:
                writer.writerow(row_data)
        
        print(f"\nExtraction completed!")
        print(f"Successfully processed {len(successful_extractions)} files")
        print(f"Failed files: {len(failed_files)}")
        if failed_files:
            print(f"Failed file list: {', '.join(failed_files)}")
        print(f"Results saved to: {output_file}")
        
    else:
        print("No data was successfully extracted")

def simple_extract_first_rows(input_directory, output_file):
    """Simplified version: directly extract first data row from each file"""
    
    csv_files = glob.glob(os.path.join(input_directory, "*.csv"))
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # Write header
        first_file = True
        
        for csv_file in csv_files:
            file_name = os.path.basename(csv_file)
            print(f"Processing file: {file_name}")
            
            try:
                with open(csv_file, 'r', encoding='utf-8') as infile:
                    reader = csv.reader(infile)
                    
                    # Skip header row
                    try:
                        header = next(reader)
                        if first_file:
                            # Write column headers for first file
                            writer.writerow(['Source_File'] + header)
                            first_file = False
                    except StopIteration:
                        print(f"  Empty file: {file_name}")
                        continue
                    
                    # Read first data row
                    try:
                        first_data_row = next(reader)
                        writer.writerow([file_name] + first_data_row)
                        print(f"  Extracted data: {first_data_row}")
                    except StopIteration:
                        print(f"  Header only, no data: {file_name}")
                        
            except Exception as e:
                print(f"  Error processing {file_name}: {e}")
    
    print(f"\nExtraction completed, results saved to: {output_file}")

def extract_max_score_rows(input_directory, output_file, score_column_index=1):
    """Extract rows with highest score from each file (assumes second column is score)"""
    
    csv_files = glob.glob(os.path.join(input_directory, "*.csv"))
    
    results = []
    
    for csv_file in csv_files:
        file_name = os.path.basename(csv_file)
        print(f"Processing file: {file_name}")
        
        try:
            with open(csv_file, 'r', encoding='utf-8') as infile:
                reader = csv.reader(infile)
                
                # Read header row
                try:
                    header = next(reader)
                except StopIteration:
                    print(f"  Empty file: {file_name}")
                    continue
                
                # Find row with highest score
                max_score = float('-inf')
                max_row = None
                
                for row in reader:
                    if len(row) > score_column_index:
                        try:
                            score = float(row[score_column_index])
                            if score > max_score:
                                max_score = score
                                max_row = row
                        except (ValueError, IndexError):
                            continue
                
                if max_row:
                    result = {
                        'Source_File': file_name,
                        'Max_Score': max_score
                    }
                    
                    # Add all column data
                    for i, value in enumerate(max_row):
                        if i < len(header):
                            result[header[i]] = value.strip() if value else ''
                    
                    results.append(result)
                    print(f"  Highest score: {max_score}, data: {max_row}")
                else:
                    print(f"  No valid score data found")
                    
        except Exception as e:
            print(f"  Error processing {file_name}: {e}")
    
    # Save results
    if results:
        # Get all possible column names
        all_columns = set()
        for result in results:
            all_columns.update(result.keys())
        
        fieldnames = ['Source_File', 'Max_Score'] + sorted([col for col in all_columns if col not in ['Source_File', 'Max_Score']])
        
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        print(f"\nExtraction completed! Found {len(results)} highest score rows")
        print(f"Results saved to: {output_file}")
    else:
        print("No valid data found")

def main():
    """Main function"""
    # Set input and output paths
    input_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/no_supfam_pipeline/filtered_results"
    output_file = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/no_supfam_pipeline/extracted_first_rows.csv"
    
    # Check if input directory exists
    if not os.path.exists(input_directory):
        print(f"Error: Input directory does not exist: {input_directory}")
        return
    
    print("Extracting first data rows from CSV files (skip headers)...")
    print(f"Input directory: {input_directory}")
    print(f"Output file: {output_file}")
    
    # Choose extraction method
    print("\nChoose extraction method:")
    print("1. Extract first data row from each file (skip headers)")
    print("2. Extract highest score row from each file (assumes second column is score)")
    print("3. Simplified extraction (fast)")
    
    choice = input("Please enter your choice (1/2/3): ").strip()
    
    if choice == '1':
        extract_first_rows_batch(input_directory, output_file)
    elif choice == '2':
        extract_max_score_rows(input_directory, output_file)
    elif choice == '3':
        simple_extract_first_rows(input_directory, output_file)
    else:
        print("Invalid choice, defaulting to method 1")
        extract_first_rows_batch(input_directory, output_file)

if __name__ == "__main__":
    main()