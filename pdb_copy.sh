#!/bin/bash

# Define the source directory
SOURCE_DIR="/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/RNAnneal_LiLab_2026-03-27/sequence/reps"

# Define the destination (relative to your current location or use an absolute path)
DEST_DIR="./pdb"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop from 0 to 19
for i in {0..150}
do
    # Construct the filename
    FILENAME="sgRNA_${i}_rank01.pdb"
    
    # Check if the file exists before copying
    if [ -f "$SOURCE_DIR/$FILENAME" ]; then
        cp "$SOURCE_DIR/$FILENAME" "$DEST_DIR/"
        echo "Copied $FILENAME to $DEST_DIR"
    else
        echo "Warning: $FILENAME not found in source directory."
    fi
done

echo "Transfer complete."
