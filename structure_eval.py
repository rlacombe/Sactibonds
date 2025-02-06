import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import numpy as np
import requests
import tempfile
import os
from typing import List, Tuple, Dict
import json

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

    def calculate_sacti_rmsd(self, distances: List[float], ideal_distance: float = 1.81) -> float:
        """
        Calculate RMSD between predicted sactibond distances and ideal distance.
        
        Args:
            distances: List of predicted distances in Angstroms
            ideal_distance: Theoretical ideal distance (default 1.81 Å)
            
        Returns:
            RMSD value in Angstroms
        """
        # Calculate squared differences from ideal distance
        squared_diffs = [(d - ideal_distance)**2 for d in distances]
        # Calculate mean and take square root
        rmsd = np.sqrt(np.mean(squared_diffs))
        return rmsd

    def evaluate_sequence(self, sequence: str, 
                         sactibond_pairs: List[Tuple[int, int]], 
                         models: List[str] = None) -> Dict[str, Dict[str, float]]:
        """
        Evaluate a sequence using specified models
        sactibond_pairs: List of (donor_position, acceptor_position) pairs
        Returns: Dict with distances and RMSD for each model
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
            
            # Calculate RMSD
            rmsd = self.calculate_sacti_rmsd(distances)
            
            results[model_name] = {
                'distances': distances,
                'sacti_rmsd': rmsd
            }
            
            # Cleanup
            os.unlink(pdb_path)
            
        return results

    def evaluate_multiple_sequences(self, sequences: Dict[str, Tuple[str, List[Tuple[int, int]]]], 
                                  models: List[str] = None) -> Dict[str, Dict[str, float]]:
        """
        Evaluate multiple sequences and calculate average metrics
        
        Args:
            sequences: Dictionary of form {
                'protein_name': (sequence, [(donor1, acceptor1), (donor2, acceptor2), ...])
            }
            models: List of model names to use
            
        Returns:
            Dictionary containing per-model metrics including:
            - average_rmsd
            - per_protein_results
        """
        if models is None:
            models = list(self.models.keys())
            
        results = {model: {
            'average_rmsd': 0.0,
            'per_protein_results': {}
        } for model in models}
        
        for protein_name, (sequence, sactibond_pairs) in sequences.items():
            print(f"\nEvaluating {protein_name}...")
            
            try:
                protein_results = self.evaluate_sequence(sequence, sactibond_pairs, models)
                
                # Store results for each model
                for model in models:
                    results[model]['per_protein_results'][protein_name] = protein_results[model]
                    
            except Exception as e:
                print(f"Error processing {protein_name}: {str(e)}")
                continue
        
        # Calculate average RMSD for each model
        for model in models:
            rmsds = [data['sacti_rmsd'] 
                    for data in results[model]['per_protein_results'].values()]
            if rmsds:
                results[model]['average_rmsd'] = np.mean(rmsds)
                results[model]['std_rmsd'] = np.std(rmsds)
            
        return results

# Update example usage
if __name__ == "__main__":
    # Load protein data from JSON
    try:
        with open('sactipeptides.json', 'r') as f:
            protein_data = json.load(f)
    except FileNotFoundError:
        print("Error: sactipeptides.json not found!")
        exit(1)
        
    # Convert JSON data to format needed by evaluate_multiple_sequences
    test_sequences = {
        name: (data['sequence'], data['sactibonds']) 
        for name, data in protein_data.items()
    }
    
    evaluator = SactibondEvaluator()
    results = evaluator.evaluate_multiple_sequences(test_sequences)
    
    # Print results
    for model, data in results.items():
        print(f"\n{model} Results:")
        print(f"Average Sacti-RMSD: {data['average_rmsd']:.2f} ± {data['std_rmsd']:.2f} Å")
        print("\nPer-protein results:")
        for protein, protein_data in data['per_protein_results'].items():
            print(f"\n{protein}:")
            print(f"  Sacti-RMSD: {protein_data['sacti_rmsd']:.2f} Å")
            print("  Individual distances:")
            for distance in protein_data['distances']:
                print(f"    {distance:.2f} Å") 