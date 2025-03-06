#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Identify amino acid mutations unique to pathogenic strains."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file path")
    parser.add_argument("-o", "--output", help="Output file path (default: stdout)")
    return parser.parse_args()

def fasta_to_df(file_path):
    """
    Parse a FASTA file and return a DataFrame where rows represent sequences
    and columns represent positions.
    
    Args:
        file_path (str): Path to the FASTA file
    
    Returns:
        pd.DataFrame: DataFrame with sequence IDs as index and positions as columns.
    """
    sequences = {}
    sequence_lengths = set()
    
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str = str(record.seq)
                sequences[record.id] = seq_str
                sequence_lengths.add(len(seq_str))
        
        # Check if all sequences have the same length
        if len(sequence_lengths) > 1:
            print("Error: Not all sequences have the same length. Ensure input is pre-aligned.")
            return pd.DataFrame()
        
        max_length = max(sequence_lengths)
        
        # Convert sequences to a DataFrame, ensuring uniform length
        df = pd.DataFrame({seq_id: list(seq.ljust(max_length, '-')) for seq_id, seq in sequences.items()}).T
        
        # Rename columns to indicate positions (1-based index)
        df.columns = [str(i + 1) for i in range(max_length)]
        
        # Extract pathogenicity labels from sequence headers
        df['pathogenic'] = df.index.str.split('|').str[0]
        
        return df
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return pd.DataFrame()

def identify_unique_mutations(df):
    """
    Identify amino acid mutations exclusive to pathogenic strains.
    
    Args:
        df (pd.DataFrame): DataFrame with sequence data and pathogenicity labels.
    
    Returns:
        dict: Dictionary with positions as keys and sets of amino acids for pathogenic and non-pathogenic strains.
    """
    unique_mutations = {}
    
    # Ensure at least one pathogenic and one non-pathogenic sequence exists
    if not ({'pathogenic', 'non_pathogenic'} <= set(df['pathogenic'])):
        print("Error: Input file must contain at least one pathogenic and one non-pathogenic sequence.")
        return unique_mutations
    
    for col in df.columns[:-1]:  # Exclude the 'pathogenic' column
        amino_acids = df.groupby('pathogenic')[col].apply(set)
        pathogenic_aa = amino_acids.get('pathogenic', set())
        non_pathogenic_aa = amino_acids.get('non_pathogenic', set())

        if not pathogenic_aa.intersection(non_pathogenic_aa):
            unique_mutations[col] = (pathogenic_aa, non_pathogenic_aa)
    
    return unique_mutations

def write_output(output_path, unique_mutations):
    """Write mutation results to file or stdout."""
    output_lines = [
        f"Position {pos}: Pathogenic -> {p_aa}, Non-Pathogenic -> {np_aa}"
        for pos, (p_aa, np_aa) in unique_mutations.items()
    ]
    
    if output_path:
        with open(output_path, "w") as f:
            f.write("\n".join(output_lines) + "\n")
    else:
        print("\n".join(output_lines))

def main():
    args = parse_args()
    sequences = fasta_to_df(args.input)
    
    if not sequences.empty:
        unique_mutations = identify_unique_mutations(sequences)
        if unique_mutations:
            write_output(args.output, unique_mutations)
    else:
        print("No valid sequences found.")

if __name__ == "__main__":
    main()