#!/bin/bash

# Output files
output1="first_20_directories_code.cpp"
output2="next_20_directories_code.cpp"
output3="remaining_directories_code.cpp"

# Ensure the output files are empty
> $output1
> $output2
> $output3

# Directory counter
counter=0

# Loop through the directories
for dir in $(find . -maxdepth 1 -type d -name "VCFX_*" | sort); do
    counter=$((counter + 1))

    # Concatenate files based on the range of directories
    if [ $counter -le 20 ]; then
        find "$dir" -type f -name "*.cpp" -exec cat {} + >> $output1
        find "$dir" -type f -name "*.h" -exec cat {} + >> $output1
    elif [ $counter -le 40 ]; then
        find "$dir" -type f -name "*.cpp" -exec cat {} + >> $output2
        find "$dir" -type f -name "*.h" -exec cat {} + >> $output2
    else
        find "$dir" -type f -name "*.cpp" -exec cat {} + >> $output3
        find "$dir" -type f -name "*.h" -exec cat {} + >> $output3
    fi

done

# Provide feedback
echo "Code files have been concatenated into $output1, $output2, and $output3."
