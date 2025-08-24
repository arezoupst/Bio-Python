#!/usr/bin/env python3
"""
Protein Structure Analysis
Simple protein structure analysis from PDB - based on your original work!
"""

from Bio import PDB
import matplotlib.pyplot as plt
import numpy as np

def download_protein_structure(pdb_id):
    """Download protein structure from PDB"""
    pdb_list = PDB.PDBList()
    pdb_parser = PDB.PDBParser()
    
    print(f"Downloading protein structure: {pdb_id}")
    
    try:
        # Download the structure
        pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
        print(f"Downloaded to: {pdb_file}")
        return pdb_file, pdb_parser
    except Exception as e:
        print(f"Error downloading structure: {e}")
        return None, None

def analyze_structure(structure, pdb_id):
    """Analyze the protein structure"""
    print(f"\nStructure Analysis: {pdb_id}")
    print("=" * 50)
    
    print(f"Number of models: {len(structure)}")
    print(f"Number of chains: {len(structure[0])}")
    print(f"Number of residues: {len(list(structure[0].get_residues()))}")
    print(f"Number of atoms: {len(list(structure[0].get_atoms()))}")
    
    return structure

def analyze_chains(structure):
    """Analyze each chain in the structure"""
    print("\nChain Analysis:")
    print("=" * 30)
    
    chain_info = []
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            residues = list(chain.get_residues())
            atoms = list(chain.get_atoms())
            
            chain_data = {
                'id': chain_id,
                'residues': len(residues),
                'atoms': len(atoms)
            }
            chain_info.append(chain_data)
            
            print(f"\nChain {chain_id}:")
            print(f"  Residues: {len(residues)}")
            print(f"  Atoms: {len(atoms)}")
            
            # Show first few residues
            print(f"  First 5 residues:")
            for i, residue in enumerate(residues[:5]):
                print(f"    {i+1}: {residue.get_resname()} {residue.get_id()[1]}")
    
    return chain_info

def calculate_statistics(structure):
    """Calculate basic statistics for the structure"""
    print("\nStatistical Analysis:")
    print("=" * 30)
    
    all_atoms = list(structure[0].get_atoms())
    all_residues = list(structure[0].get_residues())
    
    # Count atom types
    atom_types = {}
    for atom in all_atoms:
        atom_name = atom.get_name()
        atom_types[atom_name] = atom_types.get(atom_name, 0) + 1
    
    print("Atom Type Counts:")
    for atom_type, count in sorted(atom_types.items()):
        print(f"  {atom_type}: {count}")
    
    print(f"\nTotal Statistics:")
    print(f"  Total atoms: {len(all_atoms)}")
    print(f"  Total residues: {len(all_residues)}")
    print(f"  Unique atom types: {len(atom_types)}")
    
    return atom_types, all_atoms, all_residues

def visualize_structure(atom_types, chain_info):
    """Create simple visualizations"""
    print("\nCreating visualizations...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Atom type distribution
    if atom_types:
        atom_names = list(atom_types.keys())
        atom_counts = list(atom_types.values())
        
        ax1.bar(range(len(atom_names)), atom_counts, color='lightblue')
        ax1.set_xlabel('Atom Types')
        ax1.set_ylabel('Count')
        ax1.set_title('Atom Type Distribution')
        ax1.set_xticks(range(len(atom_names)))
        ax1.set_xticklabels(atom_names, rotation=45)
    
    # Chain comparison
    if chain_info:
        chain_ids = [chain['id'] for chain in chain_info]
        chain_residues = [chain['residues'] for chain in chain_info]
        
        ax2.bar(chain_ids, chain_residues, color='lightcoral')
        ax2.set_xlabel('Chain ID')
        ax2.set_ylabel('Number of Residues')
        ax2.set_title('Residues per Chain')
    
    plt.tight_layout()
    plt.show()

def main():
    """Main function"""
    print("=" * 60)
    print("Protein Structure Analysis")
    print("=" * 60)
    
    # Choose a simple protein (insulin)
    pdb_id = "1ZNI"  # Insulin structure
    print(f"Analyzing protein: {pdb_id}")
    
    # Download structure
    pdb_file, pdb_parser = download_protein_structure(pdb_id)
    if not pdb_file:
        print("Failed to download structure. Exiting.")
        return
    
    # Parse the structure
    try:
        structure = pdb_parser.get_structure(pdb_id, pdb_file)
        print(f"Structure loaded successfully!")
    except Exception as e:
        print(f"Error parsing structure: {e}")
        return
    
    # Analyze structure
    structure = analyze_structure(structure, pdb_id)
    
    # Analyze chains
    chain_info = analyze_chains(structure)
    
    # Calculate statistics
    atom_types, all_atoms, all_residues = calculate_statistics(structure)
    
    # Visualize results
    visualize_structure(atom_types, chain_info)
    
    print("\nProtein structure analysis complete!")

if __name__ == "__main__":
    main()
