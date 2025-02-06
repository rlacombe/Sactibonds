import requests
import tempfile
import os
from typing import List, Dict
import json

class PeptideFolding:
    def __init__(self):
        self.models = {
            'esmfold': self.predict_esmfold
            # We'll add others back as we implement them properly
        }
        
    def read_fasta(self, fasta_path: str) -> str:
        """Read sequence from FASTA file"""
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return sequence

    """
    ESMFold model:
    https://github.com/facebookresearch/esm
    """

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

    """
    AlphaFold 3 model: https://github.com/deepmind/alphafold
    Prediction was done via the AlphaFold webserver interface: https://alphafoldserver.com/
    Output is in mmCIF format. Seed: 42.
    """       

    def predict_structures(self, protein_data: Dict, output_dir: str = "structures"):
        """
        Predict structures for all proteins using all models
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        for protein_name, data in protein_data.items():
            print(f"\nProcessing {protein_name}...")
            sequence = self.read_fasta(data['fasta'])
            
            for model_name, predict_func in self.models.items():
                output_path = os.path.join(output_dir, f"{model_name}_{protein_name}.pdb")
                
                try:
                    predict_func(sequence, output_path)
                    print(f"Generated structure for {protein_name} using {model_name}")
                except Exception as e:
                    print(f"Error predicting structure for {protein_name} with {model_name}: {str(e)}")

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