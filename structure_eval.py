import biotite.structure.io as strucio  # Changed import to use general IO
import biotite.structure as struc
import numpy as np
import os
from typing import List, Tuple, Dict
import json

class SactibondEvaluator:
    def load_structure(self, file_path: str):
        """
        Load structure from either PDB or mmCIF file.
        """
        try:
            # Use general structure loader
            structure = strucio.load_structure(file_path)
            
            # If we have multiple models, take the first one
            if isinstance(structure, struc.AtomArrayStack):
                structure = structure[0]
                
            return structure
            
        except Exception as e:
            raise ValueError(f"Failed to load structure from {file_path}: {str(e)}")

    def calculate_distances(self, file_path: str, donor_positions: List[int], 
                          acceptor_positions: List[int]) -> List[float]:
        """
        Calculate distances between Cys-S and acceptor C-alpha atoms.
        Supports both PDB and mmCIF formats.
        """
        structure = self.load_structure(file_path)
        
        # Extract coordinates
        distances = []
        for donor, acceptor in zip(donor_positions, acceptor_positions):
            donor_mask = structure.res_id == donor
            acceptor_mask = structure.res_id == acceptor
            
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
        """
        
        squared_diffs = [(d - ideal_distance)**2 for d in distances]
        rmsd = np.sqrt(np.mean(squared_diffs))
        return rmsd

    def evaluate_structure(self, file_path: str, 
                         sactibond_pairs: List[Tuple[int, int]]) -> Dict[str, float]:
        """
        Evaluate a single structure given the path to the PDB file and the sactibond pairs.
        """

        donor_positions = [pair[0] for pair in sactibond_pairs]
        acceptor_positions = [pair[1] for pair in sactibond_pairs]
        
        distances = self.calculate_distances(file_path, donor_positions, acceptor_positions)
        rmsd = self.calculate_sacti_rmsd(distances)
        
        return {
            'distances': distances,
            'sacti_rmsd': rmsd
        }

    def evaluate_structures(self, protein_data: Dict, structures_dir: str = "structures") -> Dict:
        """
        Evaluate all structures in the structures directory.
        """
        results = {}
        
        for filename in os.listdir(structures_dir):
            if not filename.endswith(('.pdb', '.cif')):
                continue
                
            # Parse filename to get model and protein name (remove either .pdb or .cif extension)
            model_name, protein_name = filename.rsplit('.', 1)[0].split('_', 1)
            
            # Try to find the matching protein name in protein_data
            protein_key = None
            for key in protein_data.keys():
                if key.lower().replace('i', '') == protein_name.lower().replace('i', ''):
                    protein_key = key
                    break
                
            if protein_key is None:
                print(f"Warning: No matching protein data found for {protein_name}")
                continue
            
            if model_name not in results:
                results[model_name] = {
                    'average_rmsd': 0.0,
                    'std_rmsd': 0.0,
                    'per_protein_results': {}
                }
            
            try:
                file_path = os.path.join(structures_dir, filename)
                protein_results = self.evaluate_structure(
                    file_path, 
                    protein_data[protein_key]['sactibonds']
                )
                results[model_name]['per_protein_results'][protein_key] = protein_results
                
            except Exception as e:
                print(f"Error evaluating {filename}: {str(e)}")
                continue
        
        # Calculate statistics for each model
        for model_name in results:
            rmsds = [data['sacti_rmsd'] 
                    for data in results[model_name]['per_protein_results'].values()]
            if rmsds:
                results[model_name]['average_rmsd'] = np.mean(rmsds)
                results[model_name]['std_rmsd'] = np.std(rmsds)
        
        return results

if __name__ == "__main__":
    # Load protein data
    try:
        with open('sactipeptides.json', 'r') as f:
            protein_data = json.load(f)
    except FileNotFoundError:
        print("Error: sactipeptides.json not found!")
        exit(1)
        
    # Evaluate structures
    evaluator = SactibondEvaluator()
    results = evaluator.evaluate_structures(protein_data)
    
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