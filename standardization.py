import rdkit
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.MolStandardize import rdMolStandardize

### Disable RDKit informational messages ###
RDLogger.DisableLog('rdApp.info')

### Setup ###
input_smi = "top_5_AS.smi"
smiles_list = []
normalizer = rdMolStandardize.Normalizer()
uncharger = rdMolStandardize.Uncharger()
te = rdMolStandardize.TautomerEnumerator()

###################
### Standardize ###
###################

def standardize_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)
    
    if mol is not None:
    
        ### General cleanup ###
        mol = rdMolStandardize.Cleanup(mol)
        
        ### Keep largest fragment (remove salts/counterions) ###
        mol = rdMolStandardize.FragmentParent(mol)
        
        ### Normalize functional groups ###
        mol = normalizer.normalize(mol)
        
        ### Neutralize charges where appropriate ###
        mol = uncharger.uncharge(mol)
        
        ### Canonical tautomer ###
        mol = te.Canonicalize(mol)         
        
        print("Extracted ligand fixed and aligned!")
        return mol
    else:
        print("Error! Invalid SMILES!")


#####################################
### Get SMILES and convert to mol ###
#####################################

supplier = Chem.SmilesMolSupplier(
    input_smi,
    titleLine=False,
    nameColumn=-1
)

for mol in supplier:
    if mol is not None:
        smi = Chem.MolToSmiles(mol)
        smiles_list.append(smi)

for smi in smiles_list:
    if smi is not None:
        standardize_smiles(smi)

#########################################
### Save as "standardized_smiles.smi" ###
#########################################

with open("standardized_smiles.smi", "w") as f:
    for smi in smiles_list:
        f.write(smi + "\n")
        
print("Results have been successfully saved as 'standardized_smiles.smi'!")