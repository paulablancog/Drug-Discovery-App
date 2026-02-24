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
        #print("\n")
        #print(subsection)
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

        
# Build a SDQ query for PubChem for External Tables
def sdq_query(collection, where, select = "*", start = 1, limit = 1000, order = "cid, asc"):
    query = {
        "select": select, # which columnds you want back
        "collection": collection, # # which table to query
        "order": [order], # sort order ("cid, asc" means sort by cid ascending)
        "start": start, # pagination start row (1=start row)
        "limit": limit, # how many rows to return (1000 is a safe choice)
        "where": where, # the filter condition (example: cid must equal 5831)
        "width": 10000000  # allow wide fields (for long text citations for example)
    }
    # endpoint for SDQ queries
    url = f"{app.utils.URL_base}/sdq/sphinxql.cgi"
    params = {
        "infmt": "json",
        "outfmt": "json",
        "query": json.dumps(query)
    }

    for attempt in range(3): # Try 3 times to retrieve the data
        print(f"\nAttempting to retrieve External Table data from URL: {url} (Attempt {attempt + 1}/3)")
        try:
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                print(f"Successfully retrieved data from URL in attempt {attempt + 1}")
                return response.json()
            else:
                print(f"Request failed with status code {response.status_code}. Attempt {attempt + 1}/3")
        # Catch any request exceptions
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}. Attempt {attempt + 1}/3")   
    print("SDQ request failed after 3 attempts.")
    return None


def get_interactions_table(compound, page_size = 1000):
    collection = "consolidatedcompoundtarget"
    where = {"ands": [{"cid": str(compound.cid)}]} # the WHERE condition for the query

    request = sdq_query(collection, where, start = 1, limit = min(page_size,1000))
    out_set = request.get("SDQOutputSet", [])
    if not out_set:
        print("No data found in the SDQ query")
        return []
    
    block = out_set[0]
    rows = block.get("rows", []) or []
    total = int(block.get("totalCount", len(rows)))

    if len(rows) >= total:
        return rows
    
    all_rows = list(rows)
    start = 1+len(rows)

    while len(all_rows) < total:
        data = sdq_query(collection, where, start = start, limit = min(page_size,1000))
        rows_page = data.get("SDQOutputSet", [{}])[0].get("rows", []) or []
        if not rows_page:
            print("No data found in the SDQ query")
            break
        all_rows.extend(rows_page)
        start += len(rows_page)

    return all_rows
