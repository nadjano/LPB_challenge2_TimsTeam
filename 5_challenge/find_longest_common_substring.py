import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze aptamer sequences for GC content and conserved motifs.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing two aptamer sequences.")
    parser.add_argument("-o", "--output", required=True, help="Output file to store results.")
    return parser.parse_args()

def gc_content(sequence):
    """Calculate the GC content percentage of a sequence."""
    gc_count = sum(1 for base in sequence if base in "GC")
    return round((gc_count / len(sequence)) * 100, 1)

def longest_common_substring(seq1, seq2):
    """Find the longest common substring between two sequences."""
    m, n = len(seq1), len(seq2)
    matrix = [[0] * (n + 1) for _ in range(m + 1)]
    longest, end_pos = 0, 0
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1] + 1
                if matrix[i][j] > longest:
                    longest = matrix[i][j]
                    end_pos = i
    
    return seq1[end_pos - longest:end_pos]

def process_sequences(input_file, output_file):
    with open(input_file, "r") as f:
        lines = f.readlines()
    # check that there are two lines
    if len(lines) != 2:
        raise ValueError("Input file must contain exactly two lines.")
    # TODO: get the seq names from the lines
    seq1 = lines[0].strip().split()[1]
    seq2 = lines[1].strip().split()[1]

    results = [
        f"seq1 length: {len(seq1)}",
        f"seq2 length: {len(seq2)}",
        f"seq1 GC content: {gc_content(seq1)}%",
        f"seq2 GC content: {gc_content(seq2)}%",
        f"Longest common substring: {longest_common_substring(seq1, seq2)}"
    ]
    
    with open(output_file, "w") as f:
        f.write("\n".join(results) + "\n")

def main():
    args = parse_args()

    process_sequences(args.input, args.output)

if __name__ == "__main__":
    main()
