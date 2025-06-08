#!/bin/bash

echo "Hello, World!"
echo "This is a bash script."
echo "This script is located in the src/knowledge_graph directory."
echo "The current directory is: $(pwd)"

# Have all files in one variable
files=$(ls ./data/input/data_input_sclc_exp_1)

# Loop through each file
for file in $files
# python generate_knowledge_graph.py --input_file $file
do
  # Check if the file is a .txt file
  if [[ $file == *.txt ]]; then
    # Print the file name
    echo "Processing file: $file"
    # remove txt from the file name
    output_file=${file%.txt}
    # Call the Python script with the file as an argument
    python generate-graph.py --input "./data/input/data_input_sclc_exp_1/$file" --output "./data/output/data_output_sclc_exp_1/$output_file.html" > ./logs/sclc_exp_1/$output_file.log 2>&1
    # Check if the Python script was successful
    if [ $? -eq 0 ]; then
      echo "Successfully processed $file"
    else
      echo "Failed to process $file"
    fi
  fi
done