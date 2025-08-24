#!/usr/bin/env python3
"""
Nucleotide Search and Retrieval
Simple search and retrieval from biological databases - based on your original work!
"""

from Bio import Entrez, SeqIO
import requests
import time

def setup_entrez():
    """Setup Entrez for NCBI access"""
    Entrez.email = "your.email@example.com"
    print("Entrez configured for NCBI access")

def search_ncbi(search_term, max_results=5):
    """Search NCBI for nucleotide sequences"""
    print(f"Searching NCBI for: {search_term}")
    
    # Perform search
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    
    print(f"Found {len(record['IdList'])} results")
    return record['IdList']

def retrieve_sequence_details(seq_id):
    """Retrieve detailed information for a sequence"""
    try:
        print(f"Retrieving details for: {seq_id}")
        
        # Fetch sequence record
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()
        
        # Display information
        print(f"\nSequence Information:")
        print(f"ID: {seq_record.id}")
        print(f"Name: {seq_record.name}")
        print(f"Description: {seq_record.description}")
        print(f"Length: {len(seq_record.seq)} bp")
        print(f"Organism: {seq_record.annotations.get('organism', 'Unknown')}")
        
        # Show first 100 bases
        print(f"\nFirst 100 bases: {str(seq_record.seq)[:100]}...")
        
        return seq_record
        
    except Exception as e:
        print(f"Error retrieving sequence {seq_id}: {e}")
        return None

def save_to_fasta(seq_record, filename):
    """Save sequence to FASTA format"""
    if seq_record:
        SeqIO.write(seq_record, filename, "fasta")
        print(f"\nSequence saved to: {filename}")
    else:
        print("\nNo sequence to save.")

def search_with_delay(search_term, max_results=3):
    """Search NCBI with a delay to respect rate limits"""
    print(f"Searching for: {search_term}")
    
    # Search
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    
    # Wait before next request
    time.sleep(1)
    
    return record['IdList']

def demonstrate_search():
    """Demonstrate different search capabilities"""
    print("=" * 60)
    print("Demonstrating NCBI Search Capabilities")
    print("=" * 60)
    
    # Search for ubiquitin gene
    print("\n1. Searching for ubiquitin gene:")
    ubiquitin_results = search_ncbi("ubiquitin gene human", max_results=3)
    for i, seq_id in enumerate(ubiquitin_results):
        print(f"  Result {i+1}: {seq_id}")
    
    # Search for actin gene
    print("\n2. Searching for actin gene:")
    actin_results = search_with_delay("actin gene mouse", max_results=3)
    for i, seq_id in enumerate(actin_results):
        print(f"  Result {i+1}: {seq_id}")
    
    # Search for specific organism
    print("\n3. Searching for specific organism:")
    organism_results = search_with_delay("Escherichia coli plasmid", max_results=3)
    for i, seq_id in enumerate(organism_results):
        print(f"  Result {i+1}: {seq_id}")
    
    return ubiquitin_results, actin_results, organism_results

def main():
    """Main function"""
    print("=" * 60)
    print("Nucleotide Search and Retrieval")
    print("=" * 60)
    
    # Setup
    setup_entrez()
    
    # Demonstrate searches
    ubiquitin_results, actin_results, organism_results = demonstrate_search()
    
    # Retrieve details for first ubiquitin result
    if ubiquitin_results:
        print("\n" + "=" * 60)
        print("Retrieving Detailed Information")
        print("=" * 60)
        
        first_id = ubiquitin_results[0]
        seq_record = retrieve_sequence_details(first_id)
        
        # Save to FASTA
        if seq_record:
            save_to_fasta(seq_record, "retrieved_sequence.fasta")
    
    print("\nSearch and retrieval demonstration complete!")

if __name__ == "__main__":
    main()
