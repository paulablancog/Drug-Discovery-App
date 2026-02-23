# Streamlit code
import pandas as pd
import app.chem
app.interactions

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

index_url = f"{app.interactions.URL_base}/rest/pug_view/index/compound/{compound.cid}/JSON/"
index_data = app.interactions.get_json(index_url)
print(index_data)
       

# Retrieving Compound information from PubChem
#data = app.chem.retrieve_compoundURL(compound)
#print("Extracted tables: ", data)

# : CC(C[N+](C)(C)C)OC(=O)N  
# 5831 C[N+](C)(C)CCOC(=O)N.[Cl-]