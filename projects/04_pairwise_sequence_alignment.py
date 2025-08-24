#!/usr/bin/env python3
"""
Pairwise Sequence Alignment
Simple sequence alignment analysis - based on your original work!
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np

def create_sample_sequences():
    """Create sample sequences for alignment"""
    # Sample sequences (you can replace with your own data)
    seq1 = "ATGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGC"
    seq2 = "ATGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGC"
    seq3 = "ATGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGC"
    
    print("Sample Sequences:")
    print(f"Sequence 1: {seq1[:50]}...")
    print(f"Sequence 2: {seq2[:50]}...")
    print(f"Sequence 3: {seq3[:50]}...")
    print(f"\nLengths: {len(seq1)}, {len(seq2)}, {len(seq3)} bp")
    
    return seq1, seq2, seq3

def perform_global_alignment(seq1, seq2):
    """Perform global alignment between two sequences"""
    print("\nGlobal Alignment Results:")
    print("=" * 50)
    
    # Align sequence 1 vs 2
    alignments_1_2 = pairwise2.align.globalxx(seq1, seq2)
    best_alignment_1_2 = alignments_1_2[0]
    
    print("\nAlignment 1 vs 2:")
    print(format_alignment(*best_alignment_1_2))
    print(f"Score: {best_alignment_1_2.score}")
    print(f"Start: {best_alignment_1_2.start}, End: {best_alignment_1_2.end}")
    
    return best_alignment_1_2

def perform_local_alignment(seq1, seq3):
    """Perform local alignment between two sequences"""
    print("\nLocal Alignment Results:")
    print("=" * 50)
    
    # Align sequence 1 vs 3
    alignments_1_3 = pairwise2.align.localxx(seq1, seq3)
    best_alignment_1_3 = alignments_1_3[0]
    
    print("\nAlignment 1 vs 3 (Local):")
    print(format_alignment(*best_alignment_1_3))
    print(f"Score: {best_alignment_1_3.score}")
    print(f"Start: {best_alignment_1_3.start}, End: {best_alignment_1_3.end}")
    
    return best_alignment_1_3

def calculate_similarity(seq1, seq2):
    """Calculate simple similarity between two sequences"""
    min_length = min(len(seq1), len(seq2))
    matches = 0
    
    for i in range(min_length):
        if seq1[i] == seq2[i]:
            matches += 1
    
    similarity = (matches / min_length) * 100
    return similarity, matches, min_length

def analyze_similarities(seq1, seq2, seq3):
    """Calculate similarities between all sequence pairs"""
    print("\nSimilarity Analysis:")
    print("=" * 30)
    
    # Calculate similarities
    sim_1_2, matches_1_2, length_1_2 = calculate_similarity(seq1, seq2)
    sim_1_3, matches_1_3, length_1_3 = calculate_similarity(seq1, seq3)
    sim_2_3, matches_2_3, length_2_3 = calculate_similarity(seq2, seq3)
    
    print(f"Sequence 1 vs 2: {sim_1_2:.2f}% ({matches_1_2}/{length_1_2})")
    print(f"Sequence 1 vs 3: {sim_1_3:.2f}% ({matches_1_3}/{length_1_3})")
    print(f"Sequence 2 vs 3: {sim_2_3:.2f}% ({matches_2_3}/{length_2_3})")
    
    return sim_1_2, sim_1_3, sim_2_3

def visualize_similarity(sim_1_2, sim_1_3, sim_2_3):
    """Create similarity matrix visualization"""
    print("\nCreating similarity visualization...")
    
    # Create similarity matrix
    similarities = [
        [100.0, sim_1_2, sim_1_3],
        [sim_1_2, 100.0, sim_2_3],
        [sim_1_3, sim_2_3, 100.0]
    ]
    
    sequence_names = ["Seq 1", "Seq 2", "Seq 3"]
    
    plt.figure(figsize=(8, 6))
    plt.imshow(similarities, cmap='Blues', interpolation='nearest')
    plt.colorbar(label='Similarity (%)')
    plt.xticks(range(3), sequence_names)
    plt.yticks(range(3), sequence_names)
    plt.title('Sequence Similarity Matrix')
    
    # Add text annotations
    for i in range(3):
        for j in range(3):
            plt.text(j, i, f'{similarities[i][j]:.1f}%', 
                    ha='center', va='center', color='black', fontweight='bold')
    
    plt.tight_layout()
    plt.show()

def demonstrate_alignment_types(seq1, seq2):
    """Demonstrate different types of alignments"""
    print("\nDemonstrating Different Alignment Types:")
    print("=" * 50)
    
    # Global alignment
    print("\n1. Global Alignment (full sequence):")
    global_align = pairwise2.align.globalxx(seq1, seq2)[0]
    print(f"Score: {global_align.score}")
    
    # Local alignment
    print("\n2. Local Alignment (best matching region):")
    local_align = pairwise2.align.localxx(seq1, seq2)[0]
    print(f"Score: {local_align.score}")
    
    # Global with gap penalty
    print("\n3. Global Alignment with Gap Penalty:")
    gap_align = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)[0]
    print(f"Score: {gap_align.score}")
    
    return global_align, local_align, gap_align

def main():
    """Main function"""
    print("=" * 60)
    print("Pairwise Sequence Alignment")
    print("=" * 60)
    
    # Create sample sequences
    seq1, seq2, seq3 = create_sample_sequences()
    
    # Perform global alignment
    best_alignment_1_2 = perform_global_alignment(seq1, seq2)
    
    # Perform local alignment
    best_alignment_1_3 = perform_local_alignment(seq1, seq3)
    
    # Calculate similarities
    sim_1_2, sim_1_3, sim_2_3 = analyze_similarities(seq1, seq2, seq3)
    
    # Demonstrate different alignment types
    global_align, local_align, gap_align = demonstrate_alignment_types(seq1, seq2)
    
    # Visualize results
    visualize_similarity(sim_1_2, sim_1_3, sim_2_3)
    
    print("\nSequence alignment analysis complete!")

if __name__ == "__main__":
    main()
