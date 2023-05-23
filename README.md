# BioPython
This is a Python library for performing sequence analysis using the BioPython package. The library includes various tools for performing common sequence analysis tasks, such as parsing FASTA files, calculating sequence statistics, and aligning sequences.
## Installation

To install the package, run the following command:

pip install biopython-sequence-analysis

## Usage

### Parsing FASTA files

The `parse_fasta` function can be used to parse a FASTA file and return a dictionary of sequence IDs and their corresponding sequences. 

### Calculating sequence statistics
The calculate_statistics function can be used to calculate various statistics for a sequence, such as its length, GC content, and molecular weight. Here's an example of how to use it:

from biopython_sequence_analysis import calculate_statistics

sequence = "ATCGATCGATCG"
stats = calculate_statistics(sequence)
print(stats)
