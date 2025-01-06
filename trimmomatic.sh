# Load required modules
module load trimmomatic-0.36-gcc-8.2.0-siurwco
module load jdk-11.0.1-gcc-8.2.0-ref6fpq

# Define adapter file path
ADAPTER_PATH=" ~/mRNA_adapters.fa"

# Set the input and output directories
INPUT_DIR="~/libraries_folders"
OUTPUT_DIR="${INPUT_DIR}/trimmed"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each R1 .fastq.gz file
for R1_file in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    # Define corresponding R2 filename
    R2_file="${R1_file/_R1_/_R2_}"

    # Define base name for output files
    base_name=$(basename "$R1_file" _R1_001.fastq.gz)

    # Define output filenames
    filtered_R1="${OUTPUT_DIR}/${base_name}_Filtered_1P.fastq.gz"
    unfiltered_R1="${OUTPUT_DIR}/${base_name}_UnFiltered_1U.fastq.gz"
    filtered_R2="${OUTPUT_DIR}/${base_name}_Filtered_2P.fastq.gz"
    unfiltered_R2="${OUTPUT_DIR}/${base_name}_UnFiltered_2U.fastq.gz"

    # Run Trimmomatic with the paired-end command
    java -jar /home/applications/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 -phred33 \
        "$R1_file" "$R2_file" \
        "$filtered_R1" "$unfiltered_R1" "$filtered_R2" "$unfiltered_R2" \
        ILLUMINACLIP:"${ADAPTER_PATH}:2:30:10" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
