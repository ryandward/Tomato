#!/bin/bash

# Array of dataset identifiers
datasets=("SRR23702320" "SRR23702319" "SRR23702316" "SRR23702315" "SRR23702314" "SRR23702313" "SRR23702312" "SRR23702311" "SRR23702310" "SRR23702309" "SRR23702318" "SRR23702317")

# Loop through the array and download each dataset
for id in "${datasets[@]}"
do
   echo "Downloading $id"
   prefetch $id
done

echo "All downloads completed"