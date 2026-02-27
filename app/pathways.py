import pandas as pd
import pubchempy as pcp
import requests
from pathlib import Path
import json
import app.utils

compound = "Ginsenoside rb1"

# Retrieve Pathways IDs of each compound
def retrieve_pathways(compound):
    # Retrieve the Interactions json
    #rows = utils.load_json(f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json", "1.Pathways")
    rows = app.utils.load_json(f"compound_9898279_{compound}_interactionstable.json", "1.Pathways")
    if rows is None:
        print("No Pathways JSON retrieved")
        return None

    # Find the desired protein
    pathways_list = list() # a set is useful for unique protacxn values but it does not conserve order
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            protacxn_id = row.get("srctargetname")
            pathways_list.append(protacxn_id)
    pathways_clean_list = list(filter(None, pathways_list))
    print("Successfully retrieved Pathway IDs")

    # Save protein set in folder
    filename = f"pathways_data"
    #folder = f"{compound.synonyms[0]}"
    folder = "Ginsenoside rb1"
    app.utils.create_file(filename, folder)
    app.utils.write_file(filename, folder, pathways_clean_list)
    
    print(f"Saved proteins of {compound} to {filename}")
    return filename

retrieve_pathways(compound)