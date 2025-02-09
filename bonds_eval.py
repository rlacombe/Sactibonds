import biotite.structure.io as strucio
import biotite.structure as struc
import numpy as np
import os
from typing import List, Tuple, Dict, Set
import json

class BondPredictor:
    def __init__(self, max_distance: float = 5.43):
        self.max_distance = max_distance
    
    def find_potential_bonds(self, structure: struc.AtomArray) -> List[Tuple[int, int]]:
        """
        Find all potential sactibonds in a structure by looking for S atoms
        within max_distance of alpha carbons.
        """
        # Find all Cys-S atoms
        s_mask = structure.atom_name == "SG"
        s_atoms = structure[s_mask]
        
        # Find all alpha carbons
        ca_mask = structure.atom_name == "CA"
        ca_atoms = structure[ca_mask]
        
        predicted_bonds = []
        
        # For each S atom
        for s_atom in s_atoms:
            s_res_id = s_atom.res_id
            
            # Calculate distances to all alpha carbons
            distances = []
            for ca_atom in ca_atoms:
                if ca_atom.res_id != s_res_id:  # Don't bond to own residue
                    dist = struc.distance(s_atom, ca_atom)
                    distances.append((dist, ca_atom.res_id))
            
            # If we found any close enough alpha carbons
            close_carbons = [(d, rid) for d, rid in distances if d <= self.max_distance]
            if close_carbons:
                # Get the closest one
                closest = min(close_carbons)
                predicted_bonds.append((int(s_res_id), int(closest[1])))  # Convert to int for hashability
        
        return predicted_bonds
    
    def evaluate_predictions(self, predicted_bonds: List[Tuple[int, int]], 
                           reference_bonds: List[Tuple[int, int]]) -> Dict[str, float]:
        """
        Calculate precision and recall for predicted bonds.
        """
        # Convert to sets of tuples
        predicted_set = {tuple(map(int, bond)) for bond in predicted_bonds}
        reference_set = {tuple(map(int, bond)) for bond in reference_bonds}
        
        true_positives = len(predicted_set & reference_set)
        false_positives = len(predicted_set - reference_set)
        false_negatives = len(reference_set - predicted_set)
        
        precision = true_positives / (true_positives + false_positives) if predicted_set else 0
        recall = true_positives / (true_positives + false_negatives) if reference_set else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        return {
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'predicted_bonds': list(predicted_set),
            'reference_bonds': list(reference_set)
        }
    
    def evaluate_structure(self, file_path: str, reference_bonds: List[Tuple[int, int]]) -> Dict:
        """
        Evaluate a structure file (potentially with multiple models).
        Returns metrics for best model by F1 score.
        """
        # Load structure
        structure = strucio.load_structure(file_path)
        
        # Convert to list of structures if multiple models
        if isinstance(structure, struc.AtomArray):
            structures = [structure]
        else:
            structures = [structure[i] for i in range(len(structure))]
        
        best_result = None
        best_f1 = -1
        
        # Evaluate each model
        for i, model in enumerate(structures):
            try:
                predicted_bonds = self.find_potential_bonds(model)
                metrics = self.evaluate_predictions(predicted_bonds, reference_bonds)
                
                if best_result is None or metrics['f1'] > best_f1:
                    best_result = {
                        'model_index': i,
                        'predicted_bonds': predicted_bonds,
                        'metrics': metrics
                    }
                    best_f1 = metrics['f1']
            except Exception as e:
                print(f"Warning: Failed to evaluate model {i}: {str(e)}")
                continue
        
        if best_result is None:
            raise Exception("No valid models could be evaluated")
            
        return best_result

if __name__ == "__main__":
    # Load protein data
    with open('sactipeptides.json', 'r') as f:
        protein_data = json.load(f)
    
    predictor = BondPredictor()
    results = {}
    
    # Evaluate all structure files
    for filename in os.listdir('structures'):
        if not filename.endswith(('.pdb', '.cif')):
            continue
            
        # Parse filename
        model_name, protein_name = filename.rsplit('.', 1)[0].split('_', 1)
        
        # Find matching protein data
        protein_key = None
        for key in protein_data.keys():
            if key.lower().replace('i', '') == protein_name.lower().replace('i', ''):
                protein_key = key
                break
        
        if protein_key is None:
            print(f"Warning: No reference data found for {protein_name}")
            continue
        
        if model_name not in results:
            results[model_name] = {
                'per_protein_results': {},
                'average_metrics': {}
            }
        
        try:
            file_path = os.path.join('structures', filename)
            evaluation = predictor.evaluate_structure(
                file_path,
                protein_data[protein_key]['sactibonds']
            )
            results[model_name]['per_protein_results'][protein_key] = evaluation
            
        except Exception as e:
            print(f"Error evaluating {filename}: {str(e)}")
            continue
    
    # Calculate average metrics for each model
    for model_name, model_results in results.items():
        metrics = ['precision', 'recall', 'f1']
        for metric in metrics:
            values = [r['metrics'][metric] 
                     for r in model_results['per_protein_results'].values()]
            model_results['average_metrics'][metric] = {
                'mean': np.mean(values),
                'std': np.std(values)
            }
    
    # Print results
    for model_name, model_results in results.items():
        print(f"\n{model_name} Results:")
        print("Average Metrics:")
        for metric, stats in model_results['average_metrics'].items():
            print(f"  {metric}: {stats['mean']:.3f} ± {stats['std']:.3f}")
        
        print("\nPer-protein results:")
        for protein, evaluation in model_results['per_protein_results'].items():
            print(f"\n{protein}:")
            print(f"  Model: {evaluation['model_index']}")
            print(f"  Precision: {evaluation['metrics']['precision']:.3f}")
            print(f"  Recall: {evaluation['metrics']['recall']:.3f}")
            print(f"  F1 Score: {evaluation['metrics']['f1']:.3f}")
            print("  Predicted bonds:")
            for bond in evaluation['predicted_bonds']:
                print(f"    Cys{bond[0]} -> CA{bond[1]}") 