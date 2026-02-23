import pandas as pd
import pubchempy as pcp
import requests 
from urllib.parse import quote
from pathlib import Path
import json
import os
import csv
import app.chem

# PIPELINE to retrieve interactions and pathways of a compound given its cid

# 1) Download PUG View index of the compound (UTILS)
# 2) Checks if the compound has an "Interactions and Pathways" section
# 3) If it does, download that section as a JSON file
# 4) Scans JSON that points to proteins, genes and pathways (for later use PUG-REST information geneID, proteinID and pathway info + interactions)
# 5) Returns those IDs as sets
 
# Inside a Record JSON, there are Sections with TOCHeading and Section (subsections)
def get_all_sections(compound):
    # TODO -> do safely the synonym thing
    index_json = app.utils.load_index_json(f"compound_{compound.cid}_{compound.synonyms[0]}__index.json")
    
    record = index_json.get("Record", {}) 
    sections = record.get("Section", []) or []
    
    out = [] #just in case there are no sections

    # For each Section, get all the subsections
    for section in sections:
        # for each Section, get all the subsections recursively
        out.extend(get_sections(section)) # what does extend do? it adds the elements of the list returned by get_sections to the out list, instead of adding the list itself as a single element
        
    # Find Interactions and Pathways TOCHeading
    found = None
    for section in out:
        if section.get("TOCHeading", "").lower() == "interactions and pathways":
            found = section
            break

    if found:
        print("\nFound Interactions and Pathways section in the index JSON")
        data = load_interactions_and_pathways_data(compound)
        save_data(data, compound)
        return out, data # both the sections and the data of interactions
    else:
        print("\nNo Interactions and Pathways section found in this compound")
        return out, None # just the sections, no interactions data

   
# Inside a Section, there are subsections with TOCHeading and Section 
def get_sections(section):
    out = [section] 
    for subsection in section.get("Section", []):
        out.extend(get_sections(subsection)) # We call recursively get_sections to get each subsection possible until there are no more subsections
        print("\n")
        print(subsection)
    return out # returns a list of Sections inside a Record (key Sections) i guess..
   

def load_interactions_and_pathways_data(compound):
    # 1. Creates the URL
    url = f"{app.utils.URL_base}/rest/pug_view/data/compound/{compound.cid}/JSON?heading=Interactions%20and%20Pathways"
    # 2. Retrieves the index JSON data of interactions and pathways
    data = app.utils.get_json(url) # This data is going to be saved as a JSON file in the interactions_and_pathways folder 

    if data is None:
        print("\nFailed to retrieve Interactions and Pathways data JSON for compound: "+compound.synonyms[0])
        return None
    return data
    

def save_data(data, compound):
    if data is not None:
        app.utils.save_json(data, f"interactions_and_pathways/compound_{compound.cid}_{compound.synonyms[0]}__interactions_and_pathways.json")
        print("\nInteractions and Pathways data saved successfully for compound: "+compound.synonyms[0])
    else:
        print("\nNo Interactions and Pathways data to save for compound: "+compound.synonyms[0])

        
