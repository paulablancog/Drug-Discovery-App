import pandas as pd
import pubchempy as pcp
import requests 
from urllib.parse import quote
import os
import csv

# Methods 
# Step 1: Recognize the drug by SMILES code with PubChem database

# CHANGE TO SAFE MODE
def retrieve_compound(smiles_code):
    compounds = pcp.get_compounds(smiles_code, namespace='smiles') #returns a list of matches
    if not compounds: 
        return None
    compound = compounds[0]
    return compound

def compound_information(compound):
    return{
        'cid': compound.cid,
        'synonyms': compound.synonyms[0],
        'name': compound.iupac_name,
        'molecular_formula': compound.molecular_formula,
        'molecular_weight': compound.molecular_weight,
    }


# Step 2: Get Interactions and Pathways CSV from compound -> get it from the PUG REST API 
def retrieve_compoundURL(compound):
    cid = compound.cid
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"

    response = requests.get(url) #in a Request form
    if response.status_code != 200:
        print("Invalid URL. Failed to retrieve Interactions & Pathway of: "+compound.synonyms[0])
        return None
    
    print("Valid URL")

    data = response.json() #parse to JSON 
    index = data.get("Record", {}).get("Section", []) #get Sections in the json order

    section = None
    for header in index :
        if header.get("TOCHeading", "").lower() == "interactions and pathways":
            section = header
            break
    
    if section == None:
        print("Compound does not interact with anything... Check JSON")
        return []
    
    print(section)

    subdata = section.get("Section", [])
    subsection = None
    downloaded = []

    for subheader in subdata:
        subsection = subheader.get("TOCHeading", "")
        
        if subsection.lower() == "chemical-target interactions":
            print(f"\n Interaction subheader entered: {subsection}") 
            file = download_chemical_target(cid) # PUG REST
            downloaded.append(file)
        elif subsection.lower() == "drug-drug interactions":
            print(f"\n Interaction subheader entered: {subsection}")
            compound_name = compound.synonyms[0]
            file = download_drug_drug(compound_name) # DrugBank collection
            downloaded.append(file)
        elif subsection.lower() == "drug-food interactions":
            print(f"\n Interaction subheader entered: {subsection}")
            #download_drug_food(cid)
    return downloaded

        #SUBHEADER ENTERED (drug-drug, drug-protein, drug-food...)
        #for info in subheader.get("Information", []):
         #   interaction_table = info.get("Value", {}).get("ExternalTableName", "")
          #  print("\n This is the final table found:")
           # print(interaction_table)
#
 #           if interaction_table and interaction_table not in interaction_tables:
  #              download_tables(interaction_table, cid)
   #             interaction_tables.add(interaction_table) #only stored the table names in strings
    #print("Extracted tables final for loop: ", interaction_tables)
    #return interaction_tables


def download_chemical_target (cid, save_folder = "interaction_tables"):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/consolidatedcompoundtarget/CSV"
    print("Downloading Chemical-Target table in CSV format")
    
    os.makedirs(save_folder, exist_ok=True)
    
    # ERROR HANDLING
    response = requests.get(url, timeout=60)

    if response.status_code != 200:
        print("Status: ", response.status_code)
        print("PubChem says (first 1200 chars): ", response.text[:1200])
        print("Download failed for chemical-target interactions")
        return None

    file_csv = os.path.join(save_folder, f"chemical_target_{cid}_consolidatedcompoundtarget.csv")


    with open(file_csv, "wb") as f:
        f.write(response.content)

    #csv_path = os.path.join(save_folder, f"chemical_target_{cid}.csv")
    #filename_csv = txt_to_csv(file_txt, csv_path)

    print("Download completed for chemical-target interactions")
    
    return file_csv


def download_drug_drug (compound_name, save_folder = "interaction_tables"):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/collection/drugbankddi/query/name/{compound_name}/CSV"
    print("Downloading Chemical-Target table")

    response = requests.get(url)
    if response.status_code != 200:
        print("Status: ", response.status_code)
        print("Download failed for chemical-target interactions")
    
    os.makedirs(save_folder, exist_ok=True)
    filename=os.path.join(save_folder, f"drugbankddi_{compound_name}.csv")

    with open(filename, "wb") as f:
        f.write(response.content)

    print("Download completed for chemical-target interactions")
    return filename


def txt_to_csv(filename_txt, filename_csv):
    
    with open(filename_txt, "r") as txt_file:
        geneids = [line.strip() for line in txt_file if line.strip()]

    with open(filename_csv, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["GeneID"])
        for geneid in geneids:
            writer.writerow([geneid])

    print("Conversion from TXT to CSV completed for chemical-target interactions")
    return filename_csv
