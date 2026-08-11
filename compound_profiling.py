%pip install admet-ai
import csv
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import Lipinski
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP
from admet_ai import ADMETModel

### Disable RDKit informational messages ###
RDLogger.DisableLog('rdApp.*')

#############
### Setup ###
#############

supply = Chem.SmilesMolSupplier(
    'standardized_smiles.smi',
    delimiter='\t',
    titleLine=False
    )

model = ADMETModel()

#######################
### Get Descriptors ###
#######################

def calculate_descriptors(mol, admet_model=None):

    smiles = Chem.MolToSmiles(mol)
    
    descriptors = {
        "Compound": Chem.MolToSmiles(mol),
        "MW": round(Descriptors.MolWt(mol),2),
        "LogP": round(Crippen.MolLogP(mol),2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Aromatic_Ring_Number": rdMolDescriptors.CalcNumAromaticRings(mol),
        "TPSA": round(Descriptors.TPSA(mol),2),
        "RotB": Descriptors.NumRotatableBonds(mol),
        "QED": round(QED.qed(mol),2),
        "Formal_Charge": Chem.GetFormalCharge(mol), 
    }

    # Add hERG prediction
    if admet_model is not None:
        preds = admet_model.predict(smiles=smiles)
        descriptors["hERG"] = round(preds.get("hERG", None),2)
        descriptors["Caco2_Wang"] = round(preds.get("Caco2_Wang", None),2)
        descriptors["BBB_Martins"] = round(preds.get("BBB_Martins", None),2)
        descriptors["Clearance_Hepatocyte_AZ"] = round(preds.get("Clearance_Hepatocyte_AZ", None),2)
        descriptors["PAINS_alert"] = preds.get("PAINS_alert", None)
        descriptors["BRENK_alert"] = preds.get("BRENK_alert", None)

    else:
        descriptors["hERG"] = None
        descriptors["Caco2_Wang"] = None
        descriptors["BBB_Martins"] = None
        descriptors["Clearance_Hepatocyte_AZ"] = None
        descriptors["PAINS_alert"] = None
        descriptors["BRENK_alert"] = None
    
    return descriptors # Returns a dictionary containing the above defined ADMET parameters

descriptor_table = []

for i, mol in enumerate(supply, start=0):
    if supply is not None:
        descriptor_table.append(calculate_descriptors(mol, admet_model=model))
        
##############
### Output ###
##############

### Save results as csv file ###
df = pd.DataFrame(descriptor_table)
saved_path = 'standardized_smiles_admet.csv'
df.to_csv(saved_path, index=False)

print(f'The results have been successfully saved as {saved_path}!') 