import config

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from file_io import load_smiles, save_smiles


# ============================================================
# RDKit setup
# ============================================================

RDLogger.DisableLog("rdApp.*")

normalizer = rdMolStandardize.Normalizer()
uncharger = rdMolStandardize.Uncharger()
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()


# ============================================================
# Standardize individual molecule
# ============================================================

def standardize_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # --------------------------------------------------------
    # General cleanup
    # --------------------------------------------------------

    mol = rdMolStandardize.Cleanup(mol)

    # --------------------------------------------------------
    # Keep largest fragment
    # --------------------------------------------------------

    if config.REMOVE_SALTS:
        mol = rdMolStandardize.FragmentParent(mol)

    # --------------------------------------------------------
    # Normalize functional groups
    # --------------------------------------------------------

    mol = normalizer.normalize(mol)

    # --------------------------------------------------------
    # Neutralize charges where appropriate
    # --------------------------------------------------------

    if config.NEUTRALIZE:
        mol = uncharger.uncharge(mol)

    # --------------------------------------------------------
    # Canonical tautomer
    # --------------------------------------------------------

    mol = tautomer_enumerator.Canonicalize(mol)

    # --------------------------------------------------------
    # Convert to canonical SMILES
    # --------------------------------------------------------

    standardized_smiles = Chem.MolToSmiles(mol)

    return standardized_smiles


# ============================================================
# Standardize complete SMILES file
# ============================================================

def standardize_file(input_file, output_file):

    supplier = load_smiles(input_file)

    standardized_smiles = []

    total = 0
    invalid = 0

    # --------------------------------------------------------
    # Process compounds
    # --------------------------------------------------------

    for mol in supplier:

        total += 1

        if mol is None:
            invalid += 1
            continue

        smiles = Chem.MolToSmiles(mol)

        standardized = standardize_smiles(smiles)

        if standardized is not None:
            standardized_smiles.append(standardized)

    # --------------------------------------------------------
    # Remove duplicates
    # --------------------------------------------------------

    before_duplicates = len(standardized_smiles)

    if config.REMOVE_DUPLICATES:

        standardized_smiles = list(
            dict.fromkeys(standardized_smiles)
        )

    duplicates_removed = (
        before_duplicates - len(standardized_smiles)
    )

    # --------------------------------------------------------
    # Save standardized SMILES
    # --------------------------------------------------------

    save_smiles(
        standardized_smiles,
        output_file
    )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    print("\n========================================")
    print("      STANDARDIZATION COMPLETE")
    print("========================================")
    print(f"Input compounds:       {total}")
    print(f"Invalid compounds:     {invalid}")
    print(f"Duplicates removed:    {duplicates_removed}")
    print(f"Output compounds:      {len(standardized_smiles)}")
    print(f"Output file:            {output_file}")
    print("========================================\n")

    return standardized_smiles