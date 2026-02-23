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
    index_json = app.utils.load_index_json(f"compound_{compound.cid}_{compound.synonyms[0]}__index.json")
    record = index_json.get("Record", {}) 
    sections = record.get("Section", []) or []
    out = [] #just in case there are no sections

    # For each Section, get all the subsections
    for section in sections:
        # for each Section, get all the subsections recursively
        out.extend(get_sections(section)) # what does extend do? it adds the elements of the list returned by get_sections to the out list, instead of adding the list itself as a single element
    return out # out is a flat list of every section in the entire JSON index

# Inside a Section, there are subsections with TOCHeading and Section 
def get_sections(section):
    out = [section] 
    for sub in section.get("Section", []):
        out.extend(get_sections(sub))
    return out # returns a list of Sections inside a Record (key Sections) i guess..
   
   
def get_headings(index_json):
    headings = []
    record = index_json.get("Record", {})
    sections = record.get("Section", []) or [] # if it does not have any Sections, empty

    def walk(section):
        header = section.get("TOCHeading")
        if isinstance(header, str):
            headings.append(header)
        for subsection in section.get("Section", []) or []:
            walk(subsection)

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
