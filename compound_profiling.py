import os
import sys
import config
import pandas as pd

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import Lipinski
from rdkit.Chem import RDConfig
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from admet_ai import ADMETModel

from file_io import load_smiles, save_csv


# ============================================================
# RDKit setup
# ============================================================

RDLogger.DisableLog("rdApp.*")


# ============================================================
# ADMET
# ============================================================

def create_admet_model():

    return ADMETModel(
        include_physchem=False
    )


# ============================================================
# PAINS / BRENK / NIH
# ============================================================

params_pains = FilterCatalogParams()
params_pains.AddCatalog(
    FilterCatalogParams.FilterCatalogs.PAINS_A
)
catalog_pains = FilterCatalog(params_pains)


params_brenk = FilterCatalogParams()
params_brenk.AddCatalog(
    FilterCatalogParams.FilterCatalogs.BRENK
)
catalog_brenk = FilterCatalog(params_brenk)


params_nih = FilterCatalogParams()
params_nih.AddCatalog(
    FilterCatalogParams.FilterCatalogs.NIH
)
catalog_nih = FilterCatalog(params_nih)


# ============================================================
# SA Score
# ============================================================

sys.path.append(
    os.path.join(RDConfig.RDContribDir, "SA_Score")
)

import sascorer


# ============================================================
# Structural alerts
# ============================================================

def calculate_pains_brenk(mol):

    pains_match = catalog_pains.HasMatch(mol)

    brenk_match = catalog_brenk.HasMatch(mol)

    nih_match = catalog_nih.HasMatch(mol)

    return pains_match, brenk_match, nih_match


# ============================================================
# Profile one compound
# ============================================================

def profile_compound(mol, admet_model=None):

    smiles = Chem.MolToSmiles(mol)

    pains_match, brenk_match, nih_match = \
        calculate_pains_brenk(mol)

    profile = {

        "Compound": smiles,

        "MW": round(
            Descriptors.MolWt(mol), 2
        ),

        "LogP": round(
            Crippen.MolLogP(mol), 2
        ),

        "HBD":
            Lipinski.NumHDonors(mol),

        "HBA":
            Lipinski.NumHAcceptors(mol),

        "Aromatic_Ring_Number":
            rdMolDescriptors.CalcNumAromaticRings(mol),

        "TPSA": round(
            Descriptors.TPSA(mol), 2
        ),

        "RotB":
            Descriptors.NumRotatableBonds(mol),

        "QED": round(
            QED.qed(mol), 2
        ),

        "Formal_Charge":
            Chem.GetFormalCharge(mol),

        "SA_score": round(
            sascorer.calculateScore(mol), 2
        ),

        "PAINS":
            pains_match,

        "BRENK":
            brenk_match,

        "NIH":
            nih_match
    }


    # ========================================================
    # ADMET
    # ========================================================

    if admet_model is not None:

        preds = admet_model.predict(
            smiles=smiles
        )

        for endpoint in config.ADMET_ENDPOINTS:

            value = preds.get(endpoint)

            profile[endpoint] = (
                round(value, 2)
                if value is not None
                else None
            )

    else:

        for endpoint in config.ADMET_ENDPOINTS:

            profile[endpoint] = None


    return profile


# ============================================================
# Profile complete file
# ============================================================

def profile_file(
    input_file,
    admet_model,
    output_file
):

    supply = load_smiles(input_file)

    profile_table = []


    for mol in supply:

        if mol is not None:

            profile = profile_compound(
                mol,
                admet_model=admet_model
            )

            profile_table.append(profile)


    # Convert to DataFrame
    df = pd.DataFrame(profile_table)


    # Save
    save_csv(
        df,
        output_file
    )


    print(
        f"The profiling results have been "
        f"successfully saved as {output_file}!"
    )


    return df