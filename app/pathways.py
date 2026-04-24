import pandas as pd
import app.utils

# TODO: mirar como conseguir en pathways repetidos solamente 1 y con el que tenga mayor numero de proteinas -> si tienen el mismo nombre


def retrieve_pathways(compound, rows, compound_name):
    """Retrieves all the present Pathways in a compound's Pathway section. Returns the DataFrame containing all the Pathways"""
    if not rows:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound", "cid"])
    
    # Send rows that contain Pathway table JSON information 
    df_proteinspathway = retrieve_proteins_from_pathway(compound, rows, compound_name)

    return df_proteinspathway


def retrieve_proteins_from_pathway(compound, rows, compound_name):
    """Extracts the Protein subsection from a compound's Pathway and retrieves all the proteins from the Pathway"""
    # Cada Pathway ID obtenido en el .txt se busca en PubChem
    # Se saca JSON del Pathway y se va a interactions 
    # Se descargan las tablas de interactions = Proteins
    # Se vuelve a hacer target count data por cada pathway
    dfs = []
    for row in rows:
        if not isinstance(row, dict):
            continue

        pathway_id = row.get("pathwayid") or ""
        pwacc = row.get("pwacc") or ""

        if not pathway_id or not pwacc:
            continue

        # The table we really want, give me the protein table for this pathway ID
        url_proteins = pcget_pathway_protein_url(pathway_id)
        proteins_table_json = app.utils.get_json(url_proteins)
        if proteins_table_json is None:
            continue
    
        # If protein table is downloaded, convert it into a clean DataFrame
        compound_cid = getattr(compound, "cid", None)
        df_targetlist = retrieve_pathway_proteins(pwacc, compound_name, compound_cid, proteins_table_json)

        if df_targetlist is not None and not df_targetlist.empty:
            dfs.append(df_targetlist)
                  
    if not dfs:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound", "cid"])
    
    # TODO revisar esta linea
    df = pd.concat(dfs, ignore_index=True)
    df = df.drop_duplicates(subset=["uniprot_accession", "protein_name", "pathway", "compound", "cid"], keep="first").reset_index(drop=True)
    return df



def pcget_pathway_protein_url(pathwayid, start=1, limit = 10000000):
    """Creates the URL to get the Protein Section from a specific Pathway"""
    return(f"{app.utils.URL_BASE}/assay/pcget.cgi?task=pathway_protein&pathwayid={pathwayid}&start={start}&limit={limit}&infmt=json&outfmt=json")


def retrieve_pathway_proteins(pwacc, compound_name, compound_cid, pathway_json):
    """Retrieves proteins from each pathway and returns the list of proteins"""
    # Retrieve the Proteins in pathway json
    if pathway_json is None:
        return None
    
    # Find protein names & id
    rows = pathway_json.get("SDQOutputSet", [{}])[0].get("rows", []) or []
    status = pathway_json.get("SDQOutputSet", [{}])[0].get("status", {}) or []
    
    if status.get("code") != 0:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound", "cid"])
    
    # Extracts protein information one by one
    target_list = []
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            acc_id = row.get("acc") or ""
            protname = row.get("protname")  or ""
            if acc_id:
                target_list.append({"uniprot_accession": acc_id, 
                                       "protein_name": protname,
                                       "pathway": pwacc,
                                       "compound": compound_name,
                                       "cid": compound_cid})

    return pd.DataFrame(target_list) 

def group_pathways(df_pathways):
    """Groups pathway proteins so each pathway is only present once"""
    if df_pathways is None or df_pathways.empty:
        return pd.DataFrame(columns=["pathway", "n_proteins", "n_compounds", "compounds", "compound_cid_pairs", "proteins", "uniprot_accessions"])
    df = df_pathways.copy()

    for col in ["pathway", "protein_name", "uniprot_accession", "compound", "cid"]:
        if col in df.columns:
            df[col] = df[col].astype("string").fillna("").str.strip()

    df = df[df["pathway"] != ""].copy()
    
    grouped = (
        df.groupby("pathway", as_index=False)
        .agg(
            n_proteins=("uniprot_accession", lambda x: x.replace("", pd.NA).dropna().nunique()),
            n_compounds=("compound", lambda x: x.replace("", pd.NA).dropna().nunique()),
            proteins =("protein_name", lambda x: ";".join(sorted(set(v for v in x if v)))),
            compounds=("compound", lambda x: ";".join(sorted(set(v for v in x if v)))),
            uniprot_accessions=("uniprot_accession", lambda x: ";".join(sorted(set(v for v in x if v)))),
        ).sort_values(["n_compounds", "pathway"], ascending=[False, True]).reset_index(drop=True)
    )

    return grouped


def group_compounds(df_pathways, selected_pathway = None):
    if df_pathways is None or df_pathways.empty:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "count", "compounds"])
    
    df_pathwaysproteins = df_pathways.copy()

    for col in ["pathway", "uniprot_accession", "protein_name", "compound"]:
        if col in df_pathwaysproteins.columns:
            df_pathwaysproteins[col] = df_pathwaysproteins[col].fillna("").astype(str).str.strip()

    if selected_pathway is not None:
        selected_pathway = str(selected_pathway).strip()
        df_selected = df_pathwaysproteins[df_pathwaysproteins["pathway"] ==selected_pathway].copy()
    else:
        df_selected = df_pathwaysproteins.copy()

    df_uniprot = df_selected[df_selected["uniprot_accession"] != ""].copy()
                                          
    grouped = (
        df_uniprot.groupby(["uniprot_accession", "protein_name"], as_index=False)
        .agg(
            count = ("compound", "size"),
            compounds = ("compound", lambda x: ";".join(sorted(set(v for v in x if v))
                                                        ))
                                                        ).sort_values(["count","uniprot_accession", "protein_name"], ascending=[False,True, True])
                                                                                                     .reset_index(drop=True)
        )  
    return grouped                
