import biotite.structure.io as strucio
import biotite.structure as struc
import numpy as np
import os
from typing import List, Tuple, Dict, Set
import json
import argparse
from datetime import datetime
import sys

class Logger:
    def __init__(self, log_file=None):
        self.terminal = sys.stdout
        self.log_file = log_file

    def write(self, message):
        self.terminal.write(message)
        if self.log_file:
            self.log_file.write(message)

    def flush(self):
        # Needed for Python 3 compatibility
        pass

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
        
        # Calculate metrics, handling edge cases
        if len(predicted_set) > 0:
            precision = true_positives / len(predicted_set)
        else:
            precision = 0.0  # No bonds predicted
        
        if len(reference_set) > 0:
            recall = true_positives / len(reference_set)
        else:
            recall = 0.0  # No reference bonds
        
        # Calculate F1 score
        if precision + recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0  # Both precision and recall are 0
        
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
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Evaluate bond predictions from structure models.')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Print detailed results for each protein')
    parser.add_argument('-l', '--log', action='store_true',
                      help='Save output to a log file')
    args = parser.parse_args()
    
    # Setup logging if requested
    if args.log:
        # Create logs directory if it doesn't exist
        os.makedirs('logs', exist_ok=True)
        
        # Create log file with timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_path = os.path.join('logs', f'bonds_eval_{timestamp}.log')
        log_file = open(log_path, 'w')
        sys.stdout = Logger(log_file)
        print(f"Logging to: {log_path}")
    
    # Load protein data
    with open('sactipeptides.json', 'r') as f:
        protein_data = json.load(f)
    
    predictor = BondPredictor()
    results = {}
    base_dir = "structures"
    
    # Iterate through model directories
    for model_name in os.listdir(base_dir):
        model_dir = os.path.join(base_dir, model_name)
        if not os.path.isdir(model_dir):
            continue
            
        results[model_name] = {
            'per_protein_results': {},
            'average_metrics': {}
        }
        
        # Evaluate each structure file in the model directory
        for filename in os.listdir(model_dir):
            if not filename.endswith(('.pdb', '.cif')):
                continue
                
            # Get protein name (without extension) and remove model prefix if present
            protein_name = os.path.splitext(filename)[0]
            if '_' in protein_name:
                protein_name = '_'.join(protein_name.split('_')[1:])  # Remove model prefix
            
            # Find matching protein data
            protein_key = None
            for key in protein_data.keys():
                # Normalize both strings: lowercase, remove spaces and underscores
                norm_key = key.lower().replace('_', '').replace(' ', '')
                norm_name = protein_name.lower().replace('_', '').replace(' ', '')
                if norm_key == norm_name:
                    protein_key = key
                    break
            
            if protein_key is None:
                print(f"Warning: No reference data found for {protein_name}")
                continue
            
            try:
                file_path = os.path.join(model_dir, filename)
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
        
        if args.verbose:
            print("\nPer-protein results:")
            for protein, evaluation in model_results['per_protein_results'].items():
                print(f"\n{protein}:")
                print(f"  Model: {evaluation['model_index']}")
                print(f"  Metrics:")
                for metric, value in evaluation['metrics'].items():
                    if metric not in ['predicted_bonds', 'reference_bonds']:
                        print(f"    {metric}: {value:.3f}")
                print("  Predicted bonds:")
                for bond in evaluation['metrics']['predicted_bonds']:
                    print(f"    Cys{bond[0]} -> CA{bond[1]}")
                print("  Reference bonds:")
                for bond in evaluation['metrics']['reference_bonds']:
                    print(f"    Cys{bond[0]} -> CA{bond[1]}") 
    
    # Close log file if we were logging
    if args.log:
        log_file.close()
        sys.stdout = sys.stdout.terminal  # Restore normal stdout 