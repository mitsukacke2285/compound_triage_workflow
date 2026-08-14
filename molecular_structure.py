#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Visualize molecular structure of standardized molecules
import config
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
smiles_input = config.INPUT_SDF
draw_structures(smiles_input)

