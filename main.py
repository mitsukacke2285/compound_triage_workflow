from io import load_compounds
from standardization import standardize
from descriptors import calculate_descriptors
from scoring import score_compounds
from reports import export_report

def main():

    compounds = load_compounds("top_5_AS.sdf")

    compounds = standardize(compounds)

    compounds = calculate_descriptors(compounds)

    compounds = score_compounds(compounds)

    export_report(compounds)

if __name__ == "__main__":
    main()