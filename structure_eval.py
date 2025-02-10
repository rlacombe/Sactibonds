import biotite.structure.io as strucio  # Changed import to use general IO
import biotite.structure as struc
import numpy as np
import os
from typing import List, Tuple, Dict
import json
import argparse

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

    def calculate_sacti_rmsd(self, distances: List[float], ideal_distance: float = 1.8) -> float:
        """
        Calculate RMSD between predicted sactibond distances and ideal distance.
        """
        
        squared_diffs = [(d - ideal_distance)**2 for d in distances]
        rmsd = np.sqrt(np.mean(squared_diffs))
        return rmsd

    def calculate_gdt_ts(self, distances: List[float], 
                        ideal_distance: float = 1.8,
                        offsets: List[float] = [1.0, 2.0, 4.0, 8.0]) -> Tuple[float, List[float]]:
        """
        Calculate GDT_TS score for the sactibond distances.
        GDT_TS = (P1 + P2 + P4 + P8)/4
        where Px is the percentage of distances within ideal_distance + x Angstroms.
        
        Args:
            distances: List of S-CA distances for known sactibonds
            ideal_distance: The ideal sactibond distance (1.8 Å)
            offsets: List of distance offsets for each P score
        
        Returns:
            Tuple of (gdt_ts score, list of individual P scores)
        """
        scores = []
        total_bonds = len(distances)
        
        # For each offset, count bonds within ideal_distance + offset
        for offset in offsets:
            cutoff = ideal_distance + offset
            aligned = sum(1 for d in distances if d <= cutoff)
            percent = (aligned / total_bonds) * 100 if total_bonds > 0 else 0
            scores.append(percent)
        
        gdt_ts = sum(scores) / len(scores)
        return gdt_ts, scores

    def evaluate_structure(self, file_path: str, 
                         sactibond_pairs: List[Tuple[int, int]]) -> Dict[str, float]:
        """
        Evaluate a single structure given the path to the file and the sactibond pairs.
        For multi-model files (like AF2 predictions), evaluate all models and return best GDT_TS.
        """
        donor_positions = [pair[0] for pair in sactibond_pairs]
        acceptor_positions = [pair[1] for pair in sactibond_pairs]
        
        # Load structure - might be multiple models
        structure = strucio.load_structure(file_path)
        
        # If single model, wrap in list for consistent processing
        if isinstance(structure, struc.AtomArray):
            structures = [structure]
        else:
            structures = [structure[i] for i in range(len(structure))]
        
        # Evaluate each model
        model_results = []
        for model_structure in structures:
            try:
                distances = self.calculate_distances_from_structure(
                    model_structure, 
                    donor_positions, 
                    acceptor_positions
                )
                rmsd = self.calculate_sacti_rmsd(distances)
                gdt_ts, gdt_scores = self.calculate_gdt_ts(distances)
                
                model_results.append({
                    'distances': distances,
                    'sacti_rmsd': rmsd,
                    'gdt_ts': gdt_ts,
                    'gdt_scores': {
                        'P1': gdt_scores[0],  # Within 2.8 Å (1.8 + 1.0)
                        'P2': gdt_scores[1],  # Within 3.8 Å (1.8 + 2.0)
                        'P4': gdt_scores[2],  # Within 5.8 Å (1.8 + 4.0)
                        'P8': gdt_scores[3]   # Within 9.8 Å (1.8 + 8.0)
                    }
                })
            except Exception as e:
                print(f"Warning: Failed to evaluate model: {str(e)}")
                continue
        
        if not model_results:
            raise Exception("No valid models could be evaluated")
        
        # Return the best result (highest GDT_TS)
        best_result = max(model_results, key=lambda x: x['gdt_ts'])
        return best_result

    def calculate_distances_from_structure(self, structure: struc.AtomArray,
                                        donor_positions: List[int],
                                        acceptor_positions: List[int]) -> List[float]:
        """
        Calculate distances between Cys-S and acceptor C-alpha atoms for a single structure.
        """
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

    def evaluate_structures(self, protein_data: Dict, base_dir: str = "structures") -> Dict:
        """
        Evaluate all structures in the structures directory.
        Supports both:
        - Model-specific subdirectories (structures/modelname/peptide.pdb)
        - Flat directory with model prefix (structures/modelname_peptide.pdb)
        """
        results = {}
        
        # First check if we have model subdirectories
        has_subdirs = any(os.path.isdir(os.path.join(base_dir, d)) 
                         for d in os.listdir(base_dir))
        
        if has_subdirs:
            # Handle directory structure: structures/modelname/peptide.pdb
            for model_name in os.listdir(base_dir):
                model_dir = os.path.join(base_dir, model_name)
                if not os.path.isdir(model_dir):
                    continue
                    
                self._evaluate_model_structures(
                    model_name, model_dir, protein_data, results, 
                    get_protein_name=lambda f: os.path.splitext(f)[0]
                )
        else:
            # Handle flat structure: structures/modelname_peptide.pdb
            self._evaluate_model_structures(
                None, base_dir, protein_data, results,
                get_protein_name=lambda f: '_'.join(f.split('_')[1:]).rsplit('.', 1)[0]
            )
        
        return results

    def _evaluate_model_structures(self, model_name: str, directory: str, 
                                 protein_data: Dict, results: Dict,
                                 get_protein_name: callable):
        """
        Helper method to evaluate structures for a specific model.
        """
        current_model = model_name  # Initialize current_model with provided model_name
        
        for filename in os.listdir(directory):
            if not filename.endswith(('.pdb', '.cif')):
                continue
            
            # Get model name from filename if not provided
            if model_name is None:
                current_model = filename.split('_')[0]
                
            # Initialize model results if needed
            if current_model not in results:
                results[current_model] = {
                    'average_rmsd': 0.0,
                    'std_rmsd': 0.0,
                    'average_gdt_ts': 0.0,
                    'std_gdt_ts': 0.0,
                    'per_protein_results': {}
                }
            
            # Get protein name using provided function and remove model prefix if present
            protein_name = get_protein_name(filename)
            if '_' in protein_name:
                protein_name = '_'.join(protein_name.split('_')[1:])  # Remove model prefix
            
            # Try to find the matching protein name in protein_data
            protein_key = None
            for key in protein_data.keys():
                # Normalize both strings: lowercase, remove spaces and underscores
                norm_key = key.lower().replace('_', '').replace(' ', '')
                norm_name = protein_name.lower().replace('_', '').replace(' ', '')
                if norm_key == norm_name:
                    protein_key = key
                    break
            
            if protein_key is None:
                print(f"Warning: No matching protein data found for {protein_name}")
                continue
            
            try:
                file_path = os.path.join(directory, filename)
                protein_results = self.evaluate_structure(
                    file_path, 
                    protein_data[protein_key]['sactibonds']
                )
                results[current_model]['per_protein_results'][protein_key] = protein_results
                
            except Exception as e:
                print(f"Error evaluating {filename}: {str(e)}")
                continue
        
        # Calculate statistics if we have results
        if current_model in results:
            rmsds = [data['sacti_rmsd'] 
                    for data in results[current_model]['per_protein_results'].values()]
            gdt_scores = [data['gdt_ts']
                         for data in results[current_model]['per_protein_results'].values()]
            
            if rmsds:
                results[current_model]['average_rmsd'] = np.mean(rmsds)
                results[current_model]['std_rmsd'] = np.std(rmsds)
                results[current_model]['average_gdt_ts'] = np.mean(gdt_scores)
                results[current_model]['std_gdt_ts'] = np.std(gdt_scores)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Evaluate structure predictions for sactibonds.')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Print detailed results for each protein')
    args = parser.parse_args()
    
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
        print(f"Average GDT_TS: {data['average_gdt_ts']:.2f} %")
        
        if args.verbose:
            print("\nPer-protein results:")
            for protein, protein_data in data['per_protein_results'].items():
                print(f"\n{protein}:")
                print(f"  Sacti-RMSD: {protein_data['sacti_rmsd']:.2f} Å")
                print(f"  GDT_TS: {protein_data['gdt_ts']:.2f}%")
                print("  GDT scores:")
                for cutoff, score in protein_data['gdt_scores'].items():
                    print(f"    {cutoff}: {score:.2f}%")
                print("  Individual distances:")
                for distance in protein_data['distances']:
                    print(f"    {distance:.2f} Å") 