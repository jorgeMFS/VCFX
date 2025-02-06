#!/bin/bash

# Output files
output1="first_20_directories_code.cpp"
output2="next_20_directories_code.cpp"
output3="remaining_directories_code.cpp"

# Ensure the output files are empty
> "$output1"
> "$output2"
> "$output3"

# Directory counter
counter=0

# Loop through the directories matching VCFX_*
for dir in $(find . -maxdepth 1 -type d -name "VCFX_*" | sort); do
    counter=$((counter + 1))

    # Select the appropriate output file based on the current directory count
    if [ $counter -le 20 ]; then
        output="$output1"
    elif [ $counter -le 40 ]; then
        output="$output2"
    else
        output="$output3"
    fi

    # Find and process .cpp and .h files in the current directory in sorted order
    for file in $(find "$dir" -type f \( -name "*.cpp" -o -name "*.h" \) | sort); do
        # Output the filename
        echo "$file" >> "$output"
        # Then output its content
        cat "$file" >> "$output"
        # Optionally add a newline or separator between files
        echo -e "\n" >> "$output"
    done

done

# Provide feedback
echo "Code files (with filenames) have been concatenated into $output1, $output2, and $output3."
