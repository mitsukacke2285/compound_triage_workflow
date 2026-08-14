#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import config

from standardization import standardize_smiles, standardize_file
from compound_profiling import create_admet_model, profile_file
from filter import combined_filter, select_compounds
from file_io import load_csv, save_csv, save_sdf

def main():

    # ========================================================
    # 1. Standardization
    # ========================================================

    print("Starting compound standardization...")
    
    standardize_file(
        config.INPUT_SMILES,
        config.STANDARDIZED_SMILES
    )

    print(
        f"Standardized compounds saved to: "
        f"{config.STANDARDIZED_SMILES}"
    )

    # ========================================================
    # 2. Creation of ADMET Model
    # ========================================================

    print("Loading ADMET model...")
    
    admet_model = create_admet_model()

    # ========================================================
    # 3. Compound profiling
    # ========================================================

    print("Starting compound profiling...")

    profile_file(
        config.STANDARDIZED_SMILES,
        admet_model,
        config.PROFILE_OUTPUT
    )

    # ========================================================
    # 4. Load profile table
    # ========================================================

    df = load_csv(config.PROFILE_OUTPUT)

    # ========================================================
    # 5. Filter
    # ========================================================

    print("Applying structural and ADMET filters...")

    df_filtered = combined_filter(df)

    df_filtered_select = select_compounds(df_filtered)


    # ========================================================
    # 6. Save filtering results
    # ========================================================

    save_csv(
        df_filtered,
        config.FILTERED_OUTPUT
    )

    save_csv(
        df_filtered_select,
        config.FILTERED_SELECT_OUTPUT
    )

    save_sdf(
        df_filtered_select,
        config.FILTERED_COMPOUNDS_SDF
    )

    # ========================================================
    # 7. Summary
    # ========================================================

    print("\nCompound triage completed.")

    print(
        f"Total profiled: {len(df)}"
    )

    print(
        f"Passed filters: {len(df_filtered_select)}"
    )

    print(
        f"Filtered table: {config.FILTERED_OUTPUT}"
    )


if __name__ == "__main__":
    main()


# In[ ]:




