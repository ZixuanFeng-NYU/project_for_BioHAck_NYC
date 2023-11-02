#!/bin/bash

# Iterate over all subdirectories within "refined-set"
for dir in refined-set/*;
do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        # Iterate over all MOL2 files in the subdirectory
        for mol2_file in "$dir"/*.mol2;
        do
            # Extract the base filename without the extension
            base=$(basename "$mol2_file" .mol2)

            # Change the working directory to the subdirectory
            cd "$dir"

            # Use Open Babel to convert MOL2 to SMILES
            obabel -i mol2 "$base".mol2 -o smi -O "$base".smi

            # Change back to the previous working directory
            cd -
        done
    fi
