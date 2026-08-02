import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def sdf_into_smiles(comp):
    supplier = Chem.SDMolSupplier(str(comp))
    compounds = []

    # Save as a list of SMILES
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        smiles = Chem.MolToSmiles(mol, canonical=True)
        compounds.append(smiles)
    return compounds

# Input file
input_file = 'top_5_AS.sdf'
sdf_into_smiles(input_file)