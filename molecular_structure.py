# Visualize molecular structure of standardized molecules
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

### Disable RDKit informational messages ###
RDLogger.DisableLog('rdApp.*')

#######################
### Draw structures ###
#######################

def draw_structures(smiles):

    ### Get SMILES ###
    supply = Chem.SmilesMolSupplier(
    smiles,
    delimiter='\t',
    titleLine=False
    )

    mols = []
    
    for mol in supply:
        if mol is not None:
            mols.append(mol)

    ### Draw Structures ###
    img = Draw.MolsToGridImage(mols, molsPerRow=2)

    return img

### Input and Call of function ###
smiles_input = 'standardized_smiles.smi'
draw_structures(smiles_input)