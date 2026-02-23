# Streamlit code
import pandas as pd
import app.chem
import app.interactions
import app.utils

URL_base = "https://pubchem.ncbi.nlm.nih.gov"

# Asking for compound SMILES code
smiles_code = input("Enter the SMILES code: ")
compound = app.chem.retrieve_compound(smiles_code)
if compound:
    print("Compound found!")
    compound_info = app.chem.compound_information(compound)
    print("\n")
    print(compound_info)
else:
    print("No compound found, check the SMILES code entered again")

# 1. Creates the URL
index_url = f"{app.interactions.URL_base}/rest/pug_view/index/compound/{compound.cid}/JSON"
# 2. Retrieves the index JSON
index_json = app.interactions.get_json(index_url)
# 3. Saves the index JSON to a file
app.interactions.save_json(index_json, f"indexes/compound_{compound.cid}_{compound.synonyms[0]}__index.json")
# 4. Gets all sections in the index JSON
sections = []
sections = app.interactions.get_all_sections(compound)
print(sections)


# Retrieving Compound information from PubChem
#data = app.chem.retrieve_compoundURL(compound)
#print("Extracted tables: ", data)

# : CC(C[N+](C)(C)C)OC(=O)N  
# 5831 C[N+](C)(C)CCOC(=O)N.[Cl-]