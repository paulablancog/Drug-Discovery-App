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

# 1) Download PUG View index of the compound
# 2) Checks if the compound has an "Interactions and Pathways" section
# 3) If it does, download that section as a JSON file
# 4) Scans JSON that points to proteins, genes and pathways (for later use PUG-REST information geneID, proteinID and pathway info + interactions)
# 5) Returns those IDs as sets

URL_base = "https://pubchem.ncbi.nml.nih.gov"

# Util Methods

# Check problems and timeout
def get_json(url):
    for attempt in range(3):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"Request failed with status code {response.status_code}. Attempt {attempt + 1}/3")
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}. Attempt {attempt + 1}/3")   
    print("Failed to retrieve data after 3 attempts.")
    return None

# Save JSON method
def save_json(data, filename):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=4))
    return str(path)

def get_sections(section):
    out = [section] # where do I get the sections?
    for sub in section.get("Section", []):
        out.extend(get_sections(sub))
    return out
    
def get_all_sections(record_json):
    record = record_json.get("Record", {}) #why like this?
    sections = record.get("Section", []) or []
    out = []
    for section in sections:
        out.extend(get_sections(section))
    return out

def get_headings(index_json):
    headings = []
    record = index_json.get("Record", {})
    sections = record.get("Section", []) or []

    def walk(sec):
        h = sec.get("TOCHeading")
        if isinstance(h, str):
            headings.append(h)
        for child in sec.get("Section", []) or []:
            walk(child)

    for s in sections:
        walk(s)

    # unique while preserving order
    seen = set()
    out = []
    for h in headings:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return out
