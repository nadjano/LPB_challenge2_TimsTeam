#!/usr/bin/env python3

from Bio import Entrez
import json
import sys
import time

def fetch_gene_names(organism_name, email, output_file):
    # Set email for NCBI
    Entrez.email = email

    # Search for taxid for species
    print(f"Searching for species {organism_name} in NCBI database...")
    handle = Entrez.esearch(db="taxonomy", term=organism_name)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        print(f"I could not find {organism_name} in NCBI. Are you sure this species exists?")
        return

    taxid = record["IdList"][0]
    print(f"Found your species! I goes by the name {taxid} in NCBI database.")

    # Search for genes associated with taxid
    print(f"Searcing for {organism_name} genes...")
    gene_ids = []
    retstart = 0
    retmax = 10000  # max records
    while True:
        handle = Entrez.esearch(db="gene", term=f"txid{taxid}[Organism]", retstart=retstart, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()

        gene_ids.extend(record["IdList"])
        retstart += retmax

        if len(record["IdList"]) < retmax:
            break

    num_genes = len(gene_ids)
    print(f"NCBI stores {num_genes} genes for {organism_name}. Searching for their names...")

    # Fetch gene summaries and parse out gene name
    gene_names = set()
    batch_size = 500  # Number of gene IDs to fetch per API call
    delay = 0.34  # Delay between API calls not to be blocked by NCBI

    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]  # Get a batch of gene IDs
        percent_complete = (i + len(batch)) / num_genes * 100  # Calculate percentage completed
        print(f"Fetching genes {i}-{i+batch_size} ({percent_complete:.1f}% processed)")

        # Fetch summary
        handle = Entrez.esummary(db="gene", id=",".join(batch), format="json")
        result = json.loads(handle.read())
        handle.close()

        # Extract gene names from summary
        for gene_id in batch:
            gene_info = result["result"].get(gene_id, {})
            gene_name = gene_info.get("name", "Unknown")  # Use "Unknown" if name is missing
            gene_names.add(gene_name)

        # Add a delay not to be blocked by NCBI
        time.sleep(delay)

    # Write gene names to file
    print(f"Writing {len(gene_names)} gene names to {output_file}...")
    with open(output_file, "w") as f:
        for gene_name in sorted(gene_names):
            f.write(f"{gene_name}\n")

    print("Done!")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python fetch_genes.py <organism_name> <email> <output_file>")
        sys.exit(1)

    organism_name = sys.argv[1]
    email = sys.argv[2]
    output_file = sys.argv[3]

    fetch_gene_names(organism_name, email, output_file)