import config
import pandas as pd
from rdkit import Chem


#########################################
### Screen for structural liabilities ###
#########################################

def structural_filter(df):

    conditions = (
        ((df["PAINS"] == True) & config.REJECT_PAINS) |
        ((df["BRENK"] == True) & config.REJECT_BRENK) |
        ((df["NIH"] == True) & config.REJECT_NIH) |
        (df["SA_score"] > config.SA_SCORE_MAX) |
        (df["MW"] > config.MW_MAX) |
        (df["TPSA"] > config.TPSA_MAX) |
        (df["HBD"] > config.HBD_MAX) |
        (df["HBA"] > config.HBA_MAX) |
        (df["Aromatic_Ring_Number"] > config.AROMATIC_RING_NUMBER_MAX) |
        (df["RotB"] > config.ROTB_MAX) |
        (df["LogP"] < config.LOGP_MIN) |
        (df["LogP"] > config.LOGP_MAX)
    )

    df["Structural_flag"] = conditions

    return df


####################################
### Screen for ADMET liabilities ###
####################################

def admet_filter(df):

    conditions = (
        (df["hERG"] > config.HERG_MAX) |
        (df["Caco2_Wang"] > config.CACO2_WANG_MAX) |
        (df["Clearance_Hepatocyte_AZ"] >
         config.CLEARANCE_HEPATOCYTE_AZ_MAX) |
        (df["Clearance_Microsome_AZ"] >
         config.CLEARANCE_MICROSOME_AZ_MAX) |
        (df["PAMPA_NCATS"] < config.PAMPA_NCATS_MIN) |
        (df["PPBR_AZ"] > config.PPBR_AZ_MAX) |
        (df["Solubility_AqSolDB"] <
         config.SOLUBILITY_AQSOLDB_MIN)
    )

    df["ADMET_flag"] = conditions

    return df


#############################
### Combined filter       ###
#############################

def combined_filter(df):

    df = structural_filter(df)
    df = admet_filter(df)

    return df


####################################
### Select compounds for docking ###
####################################

def select_compounds(df):

    return df[
        (df["Structural_flag"] == False) &
        (df["ADMET_flag"] == False)
    ].copy()


####################################
### Save filtered compounds      ###
####################################

def save_filtered_csv(df, filename):

    df.to_csv(filename, index=False)


####################################
### Save compounds as SDF        ###
####################################

def save_filtered_sdf(df, filename):

    writer = Chem.SDWriter(filename)

    for smiles in df["Compound"]:

        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            mol.SetProp("SMILES", smiles)
            writer.write(mol)

    writer.close()