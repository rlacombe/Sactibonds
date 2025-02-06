# Sactibonds
Evaluating protein structure prediction models by how closely they match known sulfur-to-alpha-carbon post-translational modifications.

![Structures of some known sactipeptides](sactibonds.png)

### Usage

To compute the predicted structures for the peptides, run:
```bash
python peptide_folding.py
```

To evaluate the predicted structures, run:

```bash
python structure_eval.py
```

### Sources

Sactibond structures for the following peptides:

- Subtilosin A
- Thuringin H
- Thuringin CD alpha
- Thuringin CD beta

Taken from:

[Flühe L, Marahiel MA. Radical S-adenosylmethionine enzyme catalyzed thioether bond formation in sactipeptide biosynthesis. Curr Opin Chem Biol. 2013;17(4):605-612. doi:10.1016/j.cbpa.2013.06.031](https://www.sciencedirect.com/science/article/pii/S1367593113001269)


### Models

Models evaluated to date:

- ESMFold: https://github.com/facebookresearch/esm
- AlphaFold 3: https://github.com/deepmind/alphafold


