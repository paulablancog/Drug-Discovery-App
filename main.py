# Streamlit code
import pandas as pd
import app.chem
import app.interactions
import app.utils

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
index_url = f"{app.utils.URL_base}/rest/pug_view/index/compound/{compound.cid}/JSON"
# 2. Retrieves the index JSON
index_json = app.utils.get_json(index_url)
# 3. Saves the index JSON to a file
app.utils.save_json(index_json, f"indexes/compound_{compound.cid}_{compound.synonyms[0]}__index.json")
# 4. Gets all sections in the index JSON
sections = []
# 5. Finds the Interactions and Pathways section and retrieves the data of interactions
out, data = app.interactions.get_all_sections(compound)

# 6. Get the table externally
rows = app.interactions.get_interactions_table(compound)
print(len(rows))
json_path = app.utils.save_rows_json(rows, compound, f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json")
print("\nSuccessfully saved interactions table as JSON for compound: "+compound.synonyms[0])
csv_path = app.utils.save_rows_csv(rows, compound, f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.csv")
print("\nSuccessfully saved interactions table as CSV for compound: "+compound.synonyms[0])

# CC(C[N+](C)(C)C)OC(=O)N  
# 5831 C[N+](C)(C)CCOC(=O)N.[Cl-]