import requests
import json
import app.utils

# Appear in main so that the user can fetch the proteins in the interactions of the compound?
def retrieve_targets_1(compound):
    # Retrieve the Interactions json
    rows = app.utils.load_json(f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json", "1.Chemical-Target Interactions")
    if rows is None:
        print("No Interactions JSON retrieved")
        return None

    # Find the desired protein
    target_list = list() # a set is useful for unique protacxn values but it does not conserve order
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            protacxn_id = row.get("srctargetname")
            target_list.append(protacxn_id)
    target_clean_list = list(filter(None, target_list))
    print("Successfully retrieved Proteins IDs")

    # Save protein set in folder
    filename = f"protein_data"
    folder = f"{compound.synonyms[0]}"
    app.utils.create_file(filename, folder)
    app.utils.write_file(filename, folder, target_clean_list)
    
    print(f"Saved proteins of {compound} to {filename}")
    return filename