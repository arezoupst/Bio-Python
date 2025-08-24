#!/usr/bin/env python3
"""
Metagenomics Analysis
Simple metagenomics analysis - based on your original work!
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def create_sample_sequences():
    """Create sample metagenomic sequences"""
    # Sample metagenomic sequences (you can replace with your own data)
    sample_sequences = [
        "ATGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGCGTACGTAGC",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    ]
    
    print(f"Loaded {len(sample_sequences)} sample sequences")
    for i, seq in enumerate(sample_sequences):
        print(f"Sequence {i+1}: {len(seq)} bp")
    
    return sample_sequences

def quality_control(sequences, min_length=50):
    """Simple quality control for sequences"""
    print("\nQuality Control Results:")
    print("=" * 30)
    
    quality_results = []
    for i, seq in enumerate(sequences):
        if len(seq) >= min_length:
            passed = True
            reason = "Passed"
        else:
            passed = False
            reason = "Too short"
        
        quality_results.append({
            'sequence_id': i+1,
            'length': len(seq),
            'passed_qc': passed,
            'reason': reason
        })
        
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"Sequence {i+1}: {status} - {reason} ({len(seq)} bp)")
    
    return quality_results

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def analyze_gc_content(sequences):
    """Analyze GC content for all sequences"""
    print("\nGC Content Analysis:")
    print("=" * 30)
    
    gc_contents = []
    for i, seq in enumerate(sequences):
        gc_content = calculate_gc_content(seq)
        gc_contents.append(gc_content)
        print(f"Sequence {i+1}: GC Content: {gc_content:.2f}%")
    
    print(f"\nAverage GC Content: {np.mean(gc_contents):.2f}%")
    print(f"GC Content Range: {np.min(gc_contents):.2f}% - {np.max(gc_contents):.2f}%")
    
    return gc_contents

def simple_classification(gc_content):
    """Simple classification based on GC content"""
    if gc_content < 40:
        return "Low GC bacteria"
    elif gc_content > 60:
        return "High GC bacteria"
    else:
        return "Medium GC bacteria"

def classify_sequences(gc_contents):
    """Classify all sequences based on GC content"""
    print("\nSequence Classification:")
    print("=" * 30)
    
    classifications = []
    for i, gc_content in enumerate(gc_contents):
        classification = simple_classification(gc_content)
        classifications.append(classification)
        print(f"Sequence {i+1}: {classification} (GC: {gc_content:.2f}%)")
    
    # Count classifications
    classification_counts = Counter(classifications)
    print(f"\nClassification Summary:")
    for classification, count in classification_counts.items():
        print(f"{classification}: {count} sequences")
    
    return classifications, classification_counts

def create_visualizations(gc_contents, classification_counts):
    """Create visualizations of the results"""
    print("\nCreating visualizations...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # GC content distribution
    ax1.hist(gc_contents, bins=10, alpha=0.7, color='lightblue', edgecolor='black')
    ax1.axvline(np.mean(gc_contents), color='red', linestyle='--', 
               label=f'Mean: {np.mean(gc_contents):.2f}%')
    ax1.set_xlabel('GC Content (%)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('GC Content Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Classification distribution
    if classification_counts:
        class_names = list(classification_counts.keys())
        class_counts_list = list(classification_counts.values())
        
        bars = ax2.bar(class_names, class_counts_list, 
                       color=['lightcoral', 'lightblue', 'lightgreen'])
        ax2.set_title('Sequence Classification Distribution')
        ax2.set_ylabel('Number of Sequences')
        ax2.set_xlabel('Classification')
        
        # Add value labels on bars
        for bar, count in zip(bars, class_counts_list):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.show()

def additional_analysis(sequences, gc_contents):
    """Perform additional simple analyses"""
    print("\nAdditional Analysis:")
    print("=" * 30)
    
    # Length analysis
    lengths = [len(seq) for seq in sequences]
    print(f"Sequence Length Analysis:")
    print(f"  Average length: {np.mean(lengths):.2f} bp")
    print(f"  Min length: {np.min(lengths)} bp")
    print(f"  Max length: {np.max(lengths)} bp")
    
    # GC content vs length correlation
    if len(gc_contents) > 1:
        correlation = np.corrcoef(lengths, gc_contents)[0, 1]
        print(f"  GC-Length correlation: {correlation:.3f}")
    
    # Nucleotide composition
    all_sequences = ''.join(sequences)
    total_bases = len(all_sequences)
    a_count = all_sequences.count('A')
    t_count = all_sequences.count('T')
    g_count = all_sequences.count('G')
    c_count = all_sequences.count('C')
    
    print(f"\nOverall Nucleotide Composition:")
    print(f"  A: {a_count} ({a_count/total_bases*100:.1f}%)")
    print(f"  T: {t_count} ({t_count/total_bases*100:.1f}%)")
    print(f"  G: {g_count} ({g_count/total_bases*100:.1f}%)")
    print(f"  C: {c_count} ({c_count/total_bases*100:.1f}%)")

def main():
    """Main function"""
    print("=" * 60)
    print("Metagenomics Analysis")
    print("=" * 60)
    
    # Load sample sequences
    sequences = create_sample_sequences()
    
    # Quality control
    quality_results = quality_control(sequences)
    
    # GC content analysis
    gc_contents = analyze_gc_content(sequences)
    
    # Sequence classification
    classifications, classification_counts = classify_sequences(gc_contents)
    
    # Additional analysis
    additional_analysis(sequences, gc_contents)
    
    # Create visualizations
    create_visualizations(gc_contents, classification_counts)
    
    print("\nMetagenomics analysis complete!")

if __name__ == "__main__":
    main()
