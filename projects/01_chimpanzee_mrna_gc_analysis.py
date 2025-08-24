#!/usr/bin/env python3
"""
Chimpanzee mRNA GC Content Analysis
Simple analysis of chimpanzee mRNA sequences - based on your original work!
"""

from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def setup_entrez():
    """Setup Entrez for NCBI access"""
    Entrez.email = "your.email@example.com"
    print("Entrez configured for NCBI access")

def search_chimpanzee_mrna():
    """Search for chimpanzee mRNA sequences"""
    search_term = "Pan troglodytes[Organism] AND mRNA[Title]"
    print(f"Searching for: {search_term}")
    
    # Perform search
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    
    print(f"Found {len(record['IdList'])} chimpanzee mRNA sequences")
    return record['IdList']

def retrieve_sequences(seq_ids, max_sequences=5):
    """Retrieve sequence information"""
    sequences = []
    for i, seq_id in enumerate(seq_ids[:max_sequences]):
        try:
            print(f"Retrieving sequence {i+1}: {seq_id}")
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            seq_record = SeqIO.read(handle, "genbank")
            handle.close()
            
            sequence_info = {
                'id': seq_record.id,
                'description': seq_record.description,
                'length': len(seq_record.seq),
                'gc_content': GC(str(seq_record.seq))
            }
            sequences.append(sequence_info)
            
            print(f"  Length: {len(seq_record.seq)} bp, GC: {sequence_info['gc_content']:.2f}%")
            
        except Exception as e:
            print(f"Error retrieving sequence {seq_id}: {e}")
            continue
    
    return sequences

def analyze_gc_content(sequences):
    """Analyze GC content"""
    if not sequences:
        print("No sequences to analyze")
        return
    
    gc_values = [seq['gc_content'] for seq in sequences]
    lengths = [seq['length'] for seq in sequences]
    
    print(f"\nGC Content Analysis:")
    print(f"Mean GC: {np.mean(gc_values):.2f}%")
    print(f"Median GC: {np.median(gc_values):.2f}%")
    print(f"Min GC: {np.min(gc_values):.2f}%")
    print(f"Max GC: {np.max(gc_values):.2f}%")
    print(f"Total sequences: {len(sequences)}")
    
    return gc_values, lengths

def visualize_results(gc_values, lengths):
    """Create visualization"""
    if not gc_values:
        print("No data to visualize")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # GC content histogram
    ax1.hist(gc_values, bins=10, alpha=0.7, color='lightblue', edgecolor='black')
    ax1.axvline(np.mean(gc_values), color='red', linestyle='--', 
               label=f'Mean: {np.mean(gc_values):.2f}%')
    ax1.set_xlabel('GC Content (%)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Chimpanzee mRNA GC Content Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Length vs GC content
    ax2.scatter(lengths, gc_values, alpha=0.6, color='orange', s=50)
    ax2.set_xlabel('Sequence Length (bp)')
    ax2.set_ylabel('GC Content (%)')
    ax2.set_title('Length vs GC Content')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def main():
    """Main function"""
    print("=" * 60)
    print("Chimpanzee mRNA GC Content Analysis")
    print("=" * 60)
    
    # Setup
    setup_entrez()
    
    # Search for sequences
    seq_ids = search_chimpanzee_mrna()
    if not seq_ids:
        print("No sequences found. Exiting.")
        return
    
    # Retrieve sequences
    sequences = retrieve_sequences(seq_ids)
    if not sequences:
        print("No sequences retrieved. Exiting.")
        return
    
    # Analyze GC content
    gc_values, lengths = analyze_gc_content(sequences)
    
    # Visualize results
    visualize_results(gc_values, lengths)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
