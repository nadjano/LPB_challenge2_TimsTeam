#!/bin/bash

# function for info messeages
info() {
    echo "$@" >&2
}

# Display header if no arguments
if [ $# -lt 1 ]; then
cat << "EOF" >&2
  _   _  _____ ____ _____    _____                 _    _                           _            
 | \ | |/ ____|  _ \_   _|  / ____|               | |  | |                         | |           
 |  \| | |    | |_) || |   | |  __  ___ _ __   ___| |__| | __ _ _ ____   _____  ___| |_ ___ _ __ 
 | . ` | |    |  _ < | |   | | |_ |/ _ \ '_ \ / _ \  __  |/ _` | '__\ \ / / _ \/ __| __/ _ \ '__|
 | |\  | |____| |_) || |_  | |__| |  __/ | | |  __/ |  | | (_| | |   \ V /  __/\__ \ ||  __/ |   
 |_| \_|\_____|____/_____|  \_____|\___|_| |_|\___|_|  |_|\__,_|_|    \_/ \___||___/\__\___|_|   
                                                                                                 
                                                                                                 
EOF
info "=============================================================================="
info "          NCBI GeneHarvester v0.3 - A Slow NCBI Gene Name Extractor           "
info "=============================================================================="
info "Author: Tim's Team"
info "Created: March 6, 2025"
info ""
info "DESCRIPTION:"
info "  Extracts gene names from NCBI for any organism"
info ""
info "USAGE:"
info "  $0 \"Organism Name\" your@email.com [-o output_file]"
info ""
info "EXAMPLES:"
info "  $0 \"Gorilla gorilla\" your@email.com > gorilla_genes.txt"
info "  $0 \"Solanum tuberosum\" your@email.com -o potato_genes.txt"
info ""
info "DEPENDENCIES:"
info "  - curl (to talk to NCBI)"
info "  - jq (because parsing JSON with grep is for masochists)"
info "  - patience (we need to be gentle with NCBI servers)"
info "You can install curl and jq with apt (\$ apt install curl jq) or conda (\$ conda install -c conda-forge curl jq)"
info ""
exit 0
fi

ORGANISM_NAME="$1"
EMAIL="$2"
OUTPUT_FILE=""

# Check for optional output file parameter
if [ "$3" = "-o" ] && [ -n "$4" ]; then
    OUTPUT_FILE="$4"
fi

# URL encode organism name
ENCODED_ORGANISM=$(echo "$ORGANISM_NAME" | sed 's/ /%20/g')

# Clear output file if it exists
if [ -n "$OUTPUT_FILE" ]; then
    > "$OUTPUT_FILE"
fi

# Search for the organism's Taxonomy ID
info "Asking NCBI if they know the organism $ORGANISM_NAME..."
TAXID=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=${ENCODED_ORGANISM}&format=json" | jq -r '.esearchresult.idlist[0]')

if [ -z "$TAXID" ] || [ "$TAXID" = "null" ]; then
    info "NCBI does not know what is '$ORGANISM_NAME'"
    info "Maybe they lack funding. You could try EBI instead..."
    exit 1
fi

info "'$ORGANISM_NAME' goes by ID ${TAXID} in NCBI"
info "Searcing for genes..."

# Function to save gene names to file
output_genes() {
    local names="$1"
    if [ -n "$OUTPUT_FILE" ]; then
        echo "$names" >> "$OUTPUT_FILE"
    else
        echo "$names"
    fi
}

# Get gene IDs for the organism with pagination
ENCODED_TERM="txid${TAXID}%5BOrganism%3Aexp%5D"
RETMAX=5000  # Number of results per page
RETSTART=0   # Starting position
TOTAL_COUNT=0
PROCESSED=0
MORE_RESULTS=true

# First, get the total count of genes
COUNT_INFO=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${ENCODED_TERM}&retmax=0&email=${EMAIL}&format=json")
TOTAL_COUNT=$(echo "$COUNT_INFO" | jq -r '.esearchresult.count')

if [ -z "$TOTAL_COUNT" ] || [ "$TOTAL_COUNT" = "null" ] || [ "$TOTAL_COUNT" -eq 0 ]; then
    info "No genes found for this organism."
    info "Check your spelling or try another organism."
    exit 1
fi

info "Found $TOTAL_COUNT genes! This might take a while for large genomes..."

# Process genes in batches
while [ $PROCESSED -lt $TOTAL_COUNT ]; do
    info "Fetching genes $PROCESSED to $((PROCESSED + RETMAX)) of $TOTAL_COUNT..."
    
    # Get batch of gene IDs
    GENE_BATCH=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${ENCODED_TERM}&retmax=${RETMAX}&retstart=${RETSTART}&email=${EMAIL}&format=json" | jq -r '.esearchresult.idlist[]')
    
    # Check if we got any gene IDs in this batch
    if [ -z "$GENE_BATCH" ]; then
        info "No more genes found or API error. Stopping."
        break
    fi
    
    # Get gene names in chunks to avoid overwhelming the API
    echo "$GENE_BATCH" | tr ' ' '\n' | xargs -n 100 | while read -r CHUNK; do
        IDS=$(echo "$CHUNK" | tr ' ' ',')
        
        # Retry mechanism for API call
        MAX_RETRIES=3
        RETRY=0
        SUCCESS=false
        
        while [ $RETRY -lt $MAX_RETRIES ] && [ "$SUCCESS" = "false" ]; do
            GENE_INFO=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=${IDS}&email=${EMAIL}&format=json")
            
            # Check if the API call was successful
            if echo "$GENE_INFO" | jq -e '.result' > /dev/null 2>&1; then
                SUCCESS=true
            else
                RETRY=$((RETRY + 1))
                info "API call failed, retrying ($RETRY/$MAX_RETRIES)..."
                sleep 2  # Wait before retrying
            fi
        done
        
        if [ "$SUCCESS" = "true" ]; then
            # Extract gene names using jq
            NAMES=$(echo "$GENE_INFO" | jq -r '.result | del(.uids) | .[] | .name')
            
            # Output to file or stdout
            output_genes "$NAMES"
        else
            info "Failed to fetch gene info for batch after $MAX_RETRIES attempts."
        fi
        
        # Sleep to avoid overwhelming the NCBI API
        sleep 0.4
    done
    
    # Update for next batch
    PROCESSED=$((PROCESSED + RETMAX))
    RETSTART=$((RETSTART + RETMAX))
    
    # Progress indicator
    info "Processed $PROCESSED/$TOTAL_COUNT genes ($(printf "%.1f" $((PROCESSED * 100 / TOTAL_COUNT)))%)..."
done

if [ -n "$OUTPUT_FILE" ]; then
    GENE_COUNT=$(wc -l < "$OUTPUT_FILE")
    info "Gene extraction complete! $GENE_COUNT gene names saved to $OUTPUT_FILE"
else
    info "Gene extraction complete!"
fi