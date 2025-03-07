#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate mutation rate with bootstrap confidence intervals."
    )
    parser.add_argument("-i", "--input", required=True, help="Input mutation file")
    parser.add_argument("-o", "--output", help="Output file path (default: stdout)")
    return parser.parse_args()

def parse_input(file_path):
    """Parse input file and extract required parameters."""
    mutations = None
    genome_length = None
    generations = None
    bootstraps = None
    
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith('mut'):
                # Extract the part after '=' and split by commas
                mutation_values = line.split('=')[1].strip()
                mutations = [int(x) for x in mutation_values.split(',')]
            elif line.startswith('l'):
                genome_length = int(line.split('=')[1].strip())
            elif line.startswith('g'):
                generations = int(line.split('=')[1].strip())
            elif line.startswith('b'):
                bootstraps = int(line.split('=')[1].strip())
    
    # Verify all required parameters were found
    if any(param is None for param in [mutations, genome_length, generations, bootstraps]):
        missing_params = []
        if mutations is None: missing_params.append("mut")
        if genome_length is None: missing_params.append("l")
        if generations is None: missing_params.append("g")
        if bootstraps is None: missing_params.append("b")
        raise ValueError(f"Missing required parameters in input file: {', '.join(missing_params)}")
        
    return mutations, genome_length, generations, bootstraps

def calculate_mutation_rate(mutations, genome_length, generations):
    """Calculate mutation rate for a single observation."""
    return mutations / (genome_length * generations)
    
def calculate_mean_mutation_rate(mutation_rates):
    """Calculate the mean mutation rate across all observations."""
    return np.mean(mutation_rates)

def bootstrap_confidence_interval(mutations, genome_length, generations, n_resamples, confidence=0.95):
    """
    Perform bootstrap resampling to estimate confidence interval of mutation rate.
    
    Args:
        mutations: List of mutation counts
        genome_length: Genome length in base pairs
        generations: Number of generations
        n_resamples: Number of bootstrap resamples
        confidence: Confidence level (default: 0.95 for 95% CI)
        
    Returns:
        tuple: (mean_mutation_rate, (lower_bound, upper_bound))
    """
    # Array to store bootstrap means
    bootstrap_means = np.zeros(n_resamples)
    
    # Perform bootstrap resampling
    for i in range(n_resamples):
        # Sample with replacement from the original mutation counts
        resample = np.random.choice(mutations, size=len(mutations), replace=True)
        # Calculate mutation rates for each resampled mutation count
        resampled_rates = [calculate_mutation_rate(mut, genome_length, generations) for mut in resample]
        # Calculate and store mean of resampled rates
        bootstrap_means[i] = np.mean(resampled_rates)
    
    # Calculate confidence interval
    lower = np.percentile(bootstrap_means, (1-confidence)*100/2)
    upper = np.percentile(bootstrap_means, 100-(1-confidence)*100/2)
    
    # Observed mean mutation rate from original data
    observed_rates = [calculate_mutation_rate(mut, genome_length, generations) for mut in mutations]
    observed_mean = np.mean(observed_rates)
    
    return observed_mean, (lower, upper)

def write_output(output_path, mean_rate, confidence_interval):
    """Write results to output file or stdout."""
    output_text = f"Mean mutation rate observed: {mean_rate:.2e}\n95% Confidence Interval: [{confidence_interval[0]:.1e}, {confidence_interval[1]:.1e}]"
    
    if output_path:
        with open(output_path, "w") as f:
            f.write(output_text)
    else:
        print(output_text)

def main():
    args = parse_args()
    
    try:
        # Parse input file
        mutation_list, genome_length, generations, bootstraps = parse_input(args.input)
    
        # Calculate mean mutation rate and confidence interval
        mean_mutation_rate, confidence_interval = bootstrap_confidence_interval(
            mutation_list, genome_length, generations, bootstraps
        )
        
        # Write output to file or stdout
        write_output(args.output, mean_mutation_rate, confidence_interval)
        
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    main()