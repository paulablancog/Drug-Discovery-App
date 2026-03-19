import pandas as pd
import app.utils

# TODO: mirar como conseguir en pathways repetidos solamente 1 y con el que tenga mayor numero de proteinas -> si tienen el mismo nombre


def retrieve_pathways(compound):
    """Retrieves all the present Pathways in a compound's Pathway section. Returns the DataFrame containing all the Pathways"""
    compound_name = app.utils.safe_filename(compound.synonyms[0] if compound.synonyms else compound.cid)
    # Retrieve the Interactions JSON
    rows = app.utils.load_json(
        f"compound_{compound.cid}_{compound_name}_interactionstable.json", "1.Pathways")
    if rows is None:
        return None

    # Send rows that contain Pathway table JSON information 
    df_proteinspathway = retrieve_proteins_from_pathway(compound, rows)

    # Find the desired pathway
    pathways_list = []
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            pathway_name = row.get("name")
            if pathway_name:
                pathways_list.append(pathway_name)

    # Save protein set in folder
    app.utils.create_text_file("pathways_data", compound_name)
    app.utils.write_text_file("pathways_data", compound_name, pathways_list)
   
    return df_proteinspathway


def retrieve_proteins_from_pathway(compound, rows):
    """Extracts the Protein subsection from a compound's Pathway and retrieves all the proteins from the Pathway"""
    # Cada Pathway ID obtenido en el .txt se busca en PubChem
    # Se saca JSON del Pathway y se va a interactions 
    # Se descargan las tablas de interactions = Proteins
    # Se vuelve a hacer target count data por cada pathway
    dfs = []
    compound_name = app.utils.safe_filename(compound.synonyms[0] if compound.synonyms else compound.cid)
    
    for row in rows:
        if isinstance (row,dict):
            pathway_id = row.get("pathwayid") or ""
            pwacc = row.get("pwacc") or ""

            if not pwacc:
                continue

            safe_pwacc = app.utils.safe_filename(pwacc)

            # Get the Protein subsection JSON to get the table name
            # For this pwacc, does the pathway have Protein section?
            url_protein_subsection = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/pathway/{pwacc}/JSON?heading=Proteins"
            proteins_data = app.utils.get_json(url_protein_subsection)
            if proteins_data is None:
                continue
            # The table we really want, give me the protein table for this pathway ID
            url_proteins = pcget_pathway_protein_url(pathway_id)
            proteins_table_json = app.utils.get_json(url_proteins)
            if proteins_table_json is None:
                continue
    
            # If protein table is downloaded, convert it into a clean DataFrame
            df_targetlist = retrieve_pathway_proteins(safe_pwacc, compound_name, proteins_table_json)

            if df_targetlist is not None and not df_targetlist.empty:
                dfs.append(df_targetlist)
                  
    if not dfs:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound"])
    
    return pd.concat(dfs, ignore_index=True)



def pcget_pathway_protein_url(pathwayid, start=1, limit = 10000000):
    """Creates the URL to get the Protein Section from a specific Pathway"""
    return(f"{app.utils.URL_BASE}/assay/pcget.cgi?task=pathway_protein&pathwayid={pathwayid}&start={start}&limit={limit}&infmt=json&outfmt=json")


def retrieve_pathway_proteins(safe_pwacc, compound_name, pathway_json):
    """Retrieves proteins from each pathway and returns the list of proteins"""
    # Retrieve the Proteins in pathway json
    if pathway_json is None:
        return None
    
    # Find protein names & id
    rows = pathway_json.get("SDQOutputSet", [{}])[0].get("rows", []) or []
    status = pathway_json.get("SDQOutputSet", [{}])[0].get("status", {}) or []
    
    if status.get("code") != 0:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound"])
    
    # Extracts protein information one by one
    target_list = []
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            acc_id = row.get("acc") or ""
            protname = row.get("protname")  or ""
            if acc_id:
                target_list.append({"uniprot_accession": acc_id, 
                                       "protein_name": protname,
                                       "pathway": safe_pwacc,
                                       "compound": compound_name})

    return pd.DataFrame(target_list) 


"""def read_all_pathways(pathwaystxt, compound):
    Reads pathways extracted from a single compound and returns a DataFrame with pathway's ids
    compound_name = app.utils.safe_filename(compound.synonyms[0] if compound.synonyms else compound.cid)
    path = app.utils.create_folder(compound_name) / pathwaystxt
    if not path.exists():
        return pd.DataFrame(columns=["compound", "pathway_id"])

    lines = app.utils.read_file_lines(pathwaystxt, compound_name)
    pathway_list = []

    for p in lines:
        p = p.strip()
        if p:
            pathway_list.append({"compound": compound_name, "pathway_id": p})
    
    return pd.DataFrame(pathway_list)
"""

"""def map_pathways(df_pathways):
    Summarizes how many times a pathwayid is present in the list of compounds
    return  (
        df_pathways.groupby("pathway_id")
            .agg(
                total_count=("pathway_id", "size"),
                n_compounds=("compound", "nunique"),
                compounds=("compound", lambda x: ";".join(sorted(set(x)))),
            )
            .reset_index()
            .sort_values(["total_count", "n_compounds"], ascending=[False, False])
    )"""