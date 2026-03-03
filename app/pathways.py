import pandas as pd
import pubchempy as pcp
import requests
from pathlib import Path
import json
import app.utils
import app.chem
import app.interactions

# TODO: mirar como conseguir en pathways repetidos solamente 1 y con el que tenga mayor numero de proteinas -> si tienen el mismo nombre


# Retrieve Pathways IDs of each compound
def retrieve_pathways(compound):
    # Retrieve the Interactions JSON
    rows = app.utils.load_json(f"compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json", "1.Pathways")
    if rows is None:
        print("No Pathways JSON retrieved")
        return None

    # Send rows that contain Pathway table JSON information 
    retrieve_proteins_from_pathway(compound, rows)

    # Find the desired pathway
    pathways_list = list() 
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            pathway_name = row.get("name")
            pathways_list.append(pathway_name)
    pathways_clean_list = list(filter(None, pathways_list))
    print("Successfully retrieved Pathway IDs")

    # Save protein set in folder
    filename = f"pathways_data"
    folder = f"{compound.synonyms[0]}"
    app.utils.create_file(filename, folder)
    app.utils.write_file(filename, folder, pathways_clean_list)
    
    print(f"Saved proteins of {compound} to {filename}")
    return filename


def retrieve_proteins_from_pathway(compound, rows):
    # Cada Pathway ID obtenido en el .txt se busca en PubChem
    # Se saca JSON del Pathway y se va a interactions 
    # Se descargan las tablas de interactions = Proteins
    # Se vuelve a hacer target count data por cada pathway
    compound_folder = f"{compound.synonyms[0]}"
    pathway_folder = app.utils.make_subfolder(compound_folder, "Pathways")

    for row in rows:
        if isinstance (row,dict):
            pathway_id = row.get("pathwayid") or ""
            pwacc = row.get("pwacc") or ""

            if not pwacc:
                continue
            safe_pwacc = safe_filename(pwacc)

            # Save pathway JSON (it has the taxonomy)
            #url_pathway = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/pathway/{pwacc}/JSON"
            #pathway_data = app.utils.get_json(url_pathway)
            #if pathway_data is None:
             #   print("No Pathway JSON to save")
              #  continue
            #else:
             #   app.utils.save_json(pathway_data, f"{safe_pwacc}_Information.json", pathway_folder)
            

            # Get the Protein subsection JSON to get the table name
            url_protein_subsection = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/pathway/{pwacc}/JSON?heading=Proteins"
            proteins_data = app.utils.get_json(url_protein_subsection)
            if proteins_data is None:
                print("No Proteins in pathway: "+row.get("name"))
                continue
            #else:
                #app.utils.save_json(proteins_data, f"{safe_pwacc}proteinSection_data.json", pathway_folder)

            # Get from the Protein subsection the External table:pcget.cgi...
            url_proteins = pcget_pathway_protein_url(pathway_id)
            proteins_table_json = app.utils.get_json(url_proteins)

           
            if proteins_table_json is None:
                print("Failed to get proteins from pathway")
            else:
                print(f"Creation of protein table for: {safe_pwacc}...")
                #app.utils.save_json(proteins_table, f"{safe_pwacc}pcgetprotein.json", proteins_folder)
                retrieve_targets2(safe_pwacc, compound, proteins_table_json)

            print("pwacc: "+pwacc)    


def safe_filename(text):
    text = str(text)
    replace= '<>:"//|?*'
    for ch in replace:
        text = text.replace(ch,"_")
    return text.strip()

#compound = app.chem.retrieve_compound("CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C")
#retrieve_pathways(compound)

def pcget_pathway_protein_url(pathwayid, start=1, limit = 10000000):
    return(f"{app.utils.URL_base}/assay/pcget.cgi?task=pathway_protein&pathwayid={pathwayid}&start={start}&limit={limit}&infmt=json&outfmt=json")


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

def retrieve_targets2(safe_pwacc, compound, json):
    # Retrieve the Proteins in pathway json
    folder = f"{compound.synonyms[0]}"
    pathway_folder = app.utils.make_subfolder(folder, "Pathways")
    proteins_folder = app.utils.make_subfolder(pathway_folder, "Proteins")
    #json = app.utils.load_json(f"{safe_pwacc}_Information.json", proteins_folder)
    if json is None:
        print("No Interactions JSON retrieved")
        return None
    
    # Find protein names & id
    target_list = list() 
    rows = json.get("SDQOutputSet", [{}])[0].get("rows", []) or []
    status = json.get("SDQOutputSet", [{}])[0].get("status", {}) or []
    
    if status.get("code") != 0:
        print("No proteins via pcget for: "+safe_pwacc)
    else:
        for row in rows:
            if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
                acc_id = row.get("acc") or ""
                protname = row.get("protname")  or ""
                if acc_id:
                    target_list.append(f"{acc_id}\t{protname}")
                    print("Successfully retrieved Proteins IDs")
                else:
                    print("No proteins extracted from the json file")

        # Save protein set in folder
        filename = f"pathway_protein_data_{safe_pwacc}"
        app.utils.create_file(filename, proteins_folder)
        app.utils.write_file(filename, proteins_folder, target_list)
        
        print(f"Saved proteins of {compound} to {filename}")
        return filename


def read_all_pathways(pathwaystxt, compound):
    folder = compound
    path = Path(folder) / pathwaystxt

    if not path.exists():
        print("Skipping compound, no paths")
        return pd.DataFrame(columns=["compound", "pathway_id"])

    lines = app.utils.read_file_lines(pathwaystxt, folder)
    pathway_list = []

    for p in lines:
        p = p.strip()
        if not p:
            continue
        pathway_list.append({"compound": compound, "pathway_id": p})
    
    return pd.DataFrame(pathway_list)

def map_pathways(df_pathways):
    summary = (
        df_pathways.groupby("pathway_id")
            .agg(
                total_count=("pathway_id", "size"),
                n_compounds=("compound", "nunique"),
                compounds=("compound", lambda x: ";".join(sorted(set(x)))),
            )
            .reset_index()
            .sort_values(["total_count", "n_compounds"], ascending=[False, False])
    )

    return summary