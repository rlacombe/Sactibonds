import requests
import os
from typing import List, Dict
import json

class PeptideFolding:
    def __init__(self):
        self.models = {
            'esmfold': self.predict_esmfold
        }
        
    def read_fasta(self, fasta_path: str) -> str:
        """Read sequence from FASTA file"""
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return sequence

    def predict_esmfold(self, sequence: str, output_path: str):
        """
        Predict structure using ESMFold API and save to output_path
        """
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        headers = {'Content-Type': 'application/json'}
        response = requests.post(url, headers=headers, data=sequence)
        
        if response.status_code != 200:
            raise Exception(f"ESMFold prediction failed: {response.status_code}")
            
        # Save PDB file
        with open(output_path, 'wb') as f:
            f.write(response.content)

    def predict_structures(self, protein_data: Dict, base_dir: str = "structures"):
        """
        Predict structures for all proteins using all models.
        Each model's predictions are saved in a subdirectory.
        """
        for model_name, predict_func in self.models.items():
            # Create model-specific directory
            model_dir = os.path.join(base_dir, model_name)
            os.makedirs(model_dir, exist_ok=True)
            
            for protein_name, data in protein_data.items():
                print(f"\nProcessing {protein_name} with {model_name}...")
                sequence = self.read_fasta(data['fasta'])
                
                output_path = os.path.join(model_dir, f"{protein_name}.pdb")
                
                try:
                    predict_func(sequence, output_path)
                    print(f"Generated structure for {protein_name}")
                except Exception as e:
                    print(f"Error predicting structure for {protein_name}: {str(e)}")

if __name__ == "__main__":
    # Load protein data
    try:
        with open('sactipeptides.json', 'r') as f:
            protein_data = json.load(f)
    except FileNotFoundError:
        print("Error: sactipeptides.json not found!")
        exit(1)
        
    # Generate structures
    folder = PeptideFolding()
    folder.predict_structures(protein_data) 