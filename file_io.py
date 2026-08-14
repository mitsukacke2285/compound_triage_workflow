# file_io.py

import pandas as pd
from rdkit import Chem


# ============================================================
# SMILES
# ============================================================

def load_smiles(filename):
    """
    Load molecules from a SMILES file.

    Parameters
    ----------
    filename : str
        Path to the SMILES file.

    Returns
    -------
    supplier : Chem.SmilesMolSupplier
        RDKit molecule supplier.
    """

    supplier = Chem.SmilesMolSupplier(
        filename,
        delimiter="\t",
        titleLine=False
    )

    return supplier

def save_smiles(smiles_list, filename):

    with open(filename, "w") as f:

        for smiles in smiles_list:
            f.write(smiles + "\n")



# ============================================================
# CSV
# ============================================================

def load_csv(filename):
    """
    Load a CSV file into a pandas DataFrame.

    Parameters
    ----------
    filename : str
        Path to the CSV file.

    Returns
    -------
    df : pandas.DataFrame
    """

    df = pd.read_csv(filename)

    return df


def save_csv(df, filename):
    """
    Save a pandas DataFrame as a CSV file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to save.

    filename : str
        Output CSV path.
    """

    df.to_csv(filename, index=False)


# ============================================================
# SDF
# ============================================================

def load_sdf(filename):
    """
    Load molecules from an SDF file.

    Parameters
    ----------
    filename : str
        Path to the SDF file.

    Returns
    -------
    supplier : Chem.SDMolSupplier
        RDKit molecule supplier.
    """

    supplier = Chem.SDMolSupplier(filename)

    return supplier


def save_sdf(df, filename, smiles_column="Compound"):
    """
    Save molecules from a DataFrame to an SDF file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing SMILES.

    filename : str
        Output SDF path.

    smiles_column : str
        Name of the column containing SMILES.
    """

    writer = Chem.SDWriter(filename)

    for smiles in df[smiles_column]:

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            continue

        # Store SMILES as an SDF property
        mol.SetProp("SMILES", smiles)

        writer.write(mol)

    writer.close()