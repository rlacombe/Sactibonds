import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import numpy as np
import requests
import tempfile
import os
from typing import List, Tuple, Dict

class SactibondEvaluator:
    def __init__(self):
        self.models = {
            'esmfold': self.predict_esmfold,
            # Add other models as we implement them
        }
        
    def predict_esmfold(self, sequence: str) -> str:
        """
        Predict structure using ESMFold API
        Returns path to PDB file
        """
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        headers = {'Content-Type': 'application/json'}
        response = requests.post(url, headers=headers, data=sequence)
        
        if response.status_code != 200:
            raise Exception(f"ESMFold prediction failed: {response.status_code}")
            
        # Save PDB to temporary file
        temp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        temp_pdb.write(response.content)
        temp_pdb.close()
        return temp_pdb.name

    def calculate_distances(self, pdb_path: str, donor_positions: List[int], 
                          acceptor_positions: List[int]) -> List[float]:
        """
        Calculate distances between Cys-S and acceptor C-alpha atoms
        """
        # Load structure
        structure = pdb.PDBFile.read(pdb_path)
        structure = structure.get_structure()
        
        # If we have multiple models, take the first one
        if isinstance(structure, struc.AtomArrayStack):
            structure = structure[0]
        
        # Debug info
        print(f"Available residue IDs: {sorted(set(structure.res_id))}")
        print(f"Looking for donor at position {donor_positions[0]} and acceptor at {acceptor_positions[0]}")
        
        # Extract coordinates
        distances = []
        for donor, acceptor in zip(donor_positions, acceptor_positions):
            # No need to convert to 0-based indexing anymore
            donor_mask = structure.res_id == donor
            acceptor_mask = structure.res_id == acceptor
            
            print(f"\nAtoms at donor position {donor}:")
            print(set(structure.atom_name[donor_mask]))
            print(f"Atoms at acceptor position {acceptor}:")
            print(set(structure.atom_name[acceptor_mask]))
            
            # Get S atom from Cys
            s_mask = donor_mask & (structure.atom_name == "SG")
            s_atoms = structure[s_mask]
            
            # Get alpha carbon from acceptor
            ca_mask = acceptor_mask & (structure.atom_name == "CA")
            ca_atoms = structure[ca_mask]
            
            if len(s_atoms) == 0 or len(ca_atoms) == 0:
                raise Exception(f"Could not find required atoms for positions {donor} and {acceptor}")
                
            # Calculate distance
            dist = struc.distance(s_atoms[0], ca_atoms[0])
            distances.append(dist)
            
        return distances

    def evaluate_sequence(self, sequence: str, 
                         sactibond_pairs: List[Tuple[int, int]], 
                         models: List[str] = None) -> Dict[str, List[float]]:
        """
        Evaluate a sequence using specified models
        sactibond_pairs: List of (donor_position, acceptor_position) pairs
        """
        if models is None:
            models = list(self.models.keys())
            
        results = {}
        
        for model_name in models:
            if model_name not in self.models:
                raise ValueError(f"Unknown model: {model_name}")
                
            # Get structure prediction
            pdb_path = self.models[model_name](sequence)
            
            # Calculate distances
            donor_positions = [pair[0] for pair in sactibond_pairs]
            acceptor_positions = [pair[1] for pair in sactibond_pairs]
            distances = self.calculate_distances(pdb_path, donor_positions, acceptor_positions)
            
            results[model_name] = distances
            
            # Cleanup
            os.unlink(pdb_path)
            
        return results

# Example usage
if __name__ == "__main__":
    # Example sequence and known sactibond pairs: subtilosin A
    sequence = "NKGCATCSIGIACLVDGPIPDFECAGATGLGLWG" # Subtilosin A sequence. Source: https://www.uniprot.org/uniprotkb/C0HLK6/entry#sequences
    sactibond_pairs = [(4, 31), (7, 28), (13, 22)]  # Known sactibond pairs. Source: https://pubmed.ncbi.nlm.nih.gov/15035610/
    
    evaluator = SactibondEvaluator()
    results = evaluator.evaluate_sequence(sequence, sactibond_pairs)
    
    # Print results
    for model, distances in results.items():
        print(f"\n{model} results:")
        for (donor, acceptor), distance in zip(sactibond_pairs, distances):
            print(f"Distance between Cys{donor} and position {acceptor}: {distance:.2f} Å") 