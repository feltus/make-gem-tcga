#!/bin/bash

# Define parameters
MANIFEST_FILE="gdc_manifest.2025-02-26.113058.txt"  # Your downloaded manifest file
OUTPUT_FILE="blca_gene_expression_matrix.tsv"
GDC_CLIENT="gdc-client"  # Path to GDC client executable
DOWNLOAD_DIR="gdc_downloads"

# Check if manifest file exists
if [ ! -f "$MANIFEST_FILE" ]; then
    echo "Error: Manifest file $MANIFEST_FILE not found."
    exit 1
fi

# Create download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# Download files using GDC client
echo "Downloading files from GDC using manifest..."
$GDC_CLIENT download -m "$MANIFEST_FILE" -d "$DOWNLOAD_DIR"

# Find all downloaded RNA-seq files
echo "Locating RNA-seq files..."
find "$DOWNLOAD_DIR" -name "*.htseq.counts" -o -name "*.FPKM.txt" -o -name "*.rsem.genes.results" > temp_file_list.txt

if [ ! -s temp_file_list.txt ]; then
    echo "Error: No RNA-seq files found in downloaded data."
    exit 1
fi

# Create a temporary directory for processing
mkdir -p temp_gem

# Process each RNA-seq file
echo "Processing RNA-seq files..."
while read expr_file; do
    # Extract sample ID from the file path
    # GDC files typically have TCGA barcodes in their directory structure
    # Example: TCGA-BL-A0C8-01A-11R-A10J-07
    
    # Extract TCGA barcode from file path
    tcga_barcode=$(basename "$(dirname "$expr_file")" | grep -o "TCGA-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+")
    
    # If no barcode found, use directory name
    if [ -z "$tcga_barcode" ]; then
        tcga_barcode=$(basename "$(dirname "$expr_file")")
    fi
    
    echo "Processing $tcga_barcode: $expr_file"
    
    # Skip first 6 lines and extract gene_name and stranded_second data
    # Adjust column numbers based on your file format
    tail -n +7 "$expr_file" | awk '{print $1"\t"$2}' > "temp_gem/${tcga_barcode}.processed.tsv"
done < temp_file_list.txt

# Check if any files were processed
if [ ! "$(ls -A temp_gem)" ]; then
    echo "Error: No files were successfully processed. Check file formats and paths."
    rm -rf temp_gem
    exit 1
fi

# Create the gene list (first column of the matrix)
echo "Creating gene list..."
cat temp_gem/*.processed.tsv | cut -f1 | sort | uniq > temp_gem/gene_list.txt

# Create header row with sample IDs
echo -n "gene_name" > temp_gem/header.txt
for sample_file in temp_gem/*.processed.tsv; do
    sample_id=$(basename "$sample_file" .processed.tsv)
    echo -ne "\t$sample_id" >> temp_gem/header.txt
done
echo "" >> temp_gem/header.txt

# Build the matrix row by row
echo "Building expression matrix..."
while read gene_id; do
    echo -n "$gene_id" > temp_gem/current_row.txt
    
    for sample_file in temp_gem/*.processed.tsv; do
        # Find expression value for this gene in this sample
        expr_value=$(grep -w "^$gene_id" "$sample_file" | cut -f2)
        if [ -z "$expr_value" ]; then
            expr_value="NA"  # Use NA for missing values
        fi
        echo -ne "\t$expr_value" >> temp_gem/current_row.txt
    done
    
    echo "" >> temp_gem/current_row.txt
    cat temp_gem/current_row.txt >> temp_gem/matrix_body.txt
done < temp_gem/gene_list.txt

# Combine header and data rows
cat temp_gem/header.txt temp_gem/matrix_body.txt > "$OUTPUT_FILE"

# Clean up
rm -rf temp_gem
rm -f temp_file_list.txt

echo "Gene Expression Matrix created: $OUTPUT_FILE"
echo "Matrix dimensions: $(wc -l < "$OUTPUT_FILE") rows × $(head -1 "$OUTPUT_FILE" | tr '\t' '\n' | wc -l) columns"