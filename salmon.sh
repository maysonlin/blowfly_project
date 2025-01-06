module load jdk-11.0.1-gcc-8.2.0-ref6fpq

module load miniconda3-4.5.11-gcc-8.2.0-oqs2mbg
  
source activate salmon

# Set paths
INDEX_DIR="~/salmon/salmon_index"
FASTQ_DIR="~/trimmed"
OUTPUT_DIR="~/salmon_quan"

# Create the Salmon index (only run this if not already created)
if [ ! -d "$INDEX_DIR" ]; then
    salmon index \
        -t ~/rna.fna \
        -i "$INDEX_DIR" \
        -k 31
fi

# Loop over all _1P.fastq.gz files to find pairs and run Salmon quantification
for R1 in "$FASTQ_DIR"/*_1P.fastq.gz; do
    # Define the matching R2 file by replacing _1P.fastq.gz with _2P.fastq.gz
    R2="${R1/_1P.fastq.gz/_2P.fastq.gz}"
    
    # Check if the R2 file exists
    if [ -f "$R2" ]; then
        # Extract the sample name (remove path and suffix)
        SAMPLE_NAME=$(basename "$R1" "_1P.fastq.gz")
        
        # Define the output directory name following the desired format
        OUTPUT_SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}_quant"
        
        # Run Salmon quantification
        salmon quant -i "$INDEX_DIR" -l A \
            -1 "$R1" -2 "$R2" \
            -p 16 --validateMappings \
            -o "$OUTPUT_SAMPLE_DIR"
    else
        echo "Warning: Could not find matching R2 file for $R1. Skipping."
    fi
done
