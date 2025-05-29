#!/usr/bin/env python3
"""
Convert all JSON files in a folder to text files.
Each text file contains one JSON tuple per line with 'chunk' key removed.
"""

import json
import os
from pathlib import Path

def process_json_file(json_file_path, output_dir):
    """Process a single JSON file and create corresponding text file."""
    try:
        # Read the JSON file
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Create output file path (same name but .txt extension)
        json_filename = Path(json_file_path).stem
        output_file = os.path.join(output_dir, f"{json_filename}.txt")
        
        # Process and write to text file
        with open(output_file, 'w', encoding='utf-8') as f:
            processed_count = 0
            
            # Handle both single objects and arrays
            if isinstance(data, list):
                items = data
            else:
                items = [data]
            
            for item in items:
                # Remove the 'chunk' key if it exists
                if isinstance(item, dict) and 'chunk' in item:
                    del item['chunk']

                if isinstance(item, dict) and 'inferred' in item:
                    del item['inferred']
                
                # Convert to JSON string and write to file
                json_line = json.dumps(item, ensure_ascii=False)
                f.write(json_line + '\n')
                processed_count += 1
        
        print(f"✓ Processed '{json_file_path}' -> '{output_file}' ({processed_count} tuples)")
        return True
        
    except json.JSONDecodeError as e:
        print(f"✗ Error in '{json_file_path}': Invalid JSON format - {e}")
        return False
    except Exception as e:
        print(f"✗ Error processing '{json_file_path}': {e}")
        return False

def convert_folder():
    # Hardcoded folder paths - modify these as needed
    input_folder = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/output/data_output_exp_5_1"      # Folder containing JSON files
    output_folder = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/output/data_output_exp_5_1/txt"      # Folder to save text files
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Find all JSON files in input folder
    json_files = []
    if os.path.exists(input_folder):
        for filename in os.listdir(input_folder):
            if filename.lower().endswith('.json'):
                json_files.append(os.path.join(input_folder, filename))
    else:
        print(f"Error: Input folder '{input_folder}' not found.")
        return
    
    if not json_files:
        print(f"No JSON files found in '{input_folder}' folder.")
        return
    
    print(f"Found {len(json_files)} JSON file(s) to process...")
    print(f"Input folder: {input_folder}")
    print(f"Output folder: {output_folder}")
    print("-" * 50)
    
    # Process each JSON file
    successful = 0
    failed = 0
    
    for json_file in json_files:
        if process_json_file(json_file, output_folder):
            successful += 1
        else:
            failed += 1
    
    print("-" * 50)
    print(f"Conversion complete!")
    print(f"Successfully processed: {successful} files")
    if failed > 0:
        print(f"Failed to process: {failed} files")
    print(f"Output files saved in: {output_folder}")

if __name__ == "__main__":
    convert_folder()