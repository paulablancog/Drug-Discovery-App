import requests
import json
import app.utils

compound = "Ginsenoside rb1"

# Appear in main so that the user can fetch the proteins in the interactions of the compound?
def retrieve_proteins(compound):
    # Retrieve the Interactions json
    rows = app.utils.load_json(f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json", "1.Chemical-Target Interactions")
    if rows is None:
        print("No Interactions JSON retrieved")
        return None

    # Find the desired protein
    proteins_list = list() # a set is useful for unique protacxn values but it does not conserve order
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            protacxn_id = row.get("protacxn")
            proteins_list.append(protacxn_id)
    print("Successfully retrieved Proteins IDs")

    # Save protein set in folder
    filename = f"protein_data"
    folder = f"{compound.synonyms[0]}"
    app.utils.create_file(filename, folder)
    app.utils.write_file(proteins_list, filename, folder)
    
    print(f"Saved proteins of {compound} to {filename}")
    return filename

retrieve_proteins(compound)