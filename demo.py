#!/usr/bin/env python3
"""
Biopython Tools Demo
Simple demonstration of the repository structure and projects
"""

import os
import sys

def print_header():
    """Print a nice header"""
    print("=" * 70)
    print("🐍 Biopython Tools - Clean & Simple Repository")
    print("=" * 70)
    print("📚 Perfect for Bioinformatics PhD Applications!")
    print("=" * 70)

def show_repository_structure():
    """Show the repository structure"""
    print("\n📁 Repository Structure:")
    print("-" * 40)
    
    structure = {
        "📂 Tools & Examples/": "Your original bioinformatics tools and projects",
        "📂 projects/": "5 simple Python script projects (based on your work)",
        "📂 data/": "Example datasets and FASTA files",
        "📂 notebooks/": "Jupyter notebooks for demonstrations"
    }
    
    for folder, description in structure.items():
        print(f"{folder:<25} {description}")

def show_original_projects():
    """Show your original projects (preserved as tutorials)"""
    print("\n🔬 Your Original Projects (Preserved as Tutorials):")
    print("-" * 50)
    
    original_projects = [
        "1. Chimpanzee mRNA GC Content Analysis",
        "2. Nucleotide Search and Retrieval (Entrez, Expassy, Swissprot)",
        "3. Protein Structure Analysis (PDB)",
        "4. Pairwise Sequence Alignment"
    ]
    
    for project in original_projects:
        print(f"   {project}")

def show_simple_projects():
    """Show the 5 simple Python script projects"""
    print("\n🚀 Simple Python Script Projects:")
    print("-" * 40)
    
    projects = [
        "01_chimpanzee_mrna_gc_analysis.py - Simple NCBI search & GC analysis",
        "02_nucleotide_search_retrieval.py - Database search & retrieval",
        "03_protein_structure_analysis.py - PDB structure analysis",
        "04_pairwise_sequence_alignment.py - Sequence alignment & similarity",
        "05_metagenomics_analysis.py - Quality control & classification"
    ]
    
    for project in projects:
        print(f"   {project}")

def show_usage_instructions():
    """Show how to use the projects"""
    print("\n💡 How to Use:")
    print("-" * 20)
    
    print("\n🔍 Your Original Tutorials:")
    print("   • Use your existing notebooks in Tools & Examples/")
    print("   • They work exactly as before, just cleaner")
    
    print("\n🐍 Simple Python Projects:")
    print("   • Run any script in projects/ directory:")
    print("   • python projects/01_chimpanzee_mrna_gc_analysis.py")
    print("   • python projects/02_nucleotide_search_retrieval.py")
    print("   • python projects/03_protein_structure_analysis.py")
    print("   • python projects/04_pairwise_sequence_alignment.py")
    print("   • python projects/05_metagenomics_analysis.py")
    
    print("\n📊 Each script is:")
    print("   • Self-contained and easy to understand")
    print("   • Based on your original work but simplified")
    print("   • Ready to run with clear output")
    print("   • Perfect for demonstrating skills")

def show_features():
    """Show repository features"""
    print("\n✨ What Makes This Clean & Simple:")
    print("-" * 40)
    
    features = [
        "✅ Preserves Your Original Work - All projects maintained as tutorials",
        "✅ Simple Python Scripts - 5 easy-to-understand scripts based on your work",
        "✅ Clean Code - Well-organized but not over-engineered",
        "✅ Easy to Understand - Simple functions and clear documentation",
        "✅ Academic Ready - Suitable for PhD applications without being overwhelming"
    ]
    
    for feature in features:
        print(f"   {feature}")

def main():
    """Main function"""
    print_header()
    
    show_repository_structure()
    show_original_projects()
    show_simple_projects()
    show_usage_instructions()
    show_features()
    
    print("\n" + "=" * 70)
    print("🎯 Ready for your PhD application!")
    print("   Simple, clean, and professional!")
    print("=" * 70)

if __name__ == "__main__":
    main()
