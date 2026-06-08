import pandas as pd
import app.proteins
import app.interactions
import app.utils

# Step 1: Retrieve pathways rows already retrieved from PubChem for the identified compound
# Step 2: For each pathway ID, request the pathway-protein table (ExternalTableName)
# Step 3: Extract UniProt accessions and protein names for each pathway
# Step 4: Enrich protein information with UniProt metadata: protein name, gene symbol, taxid, taxname
# Step 5: Return a DataFrame of pathway-protein-compound relationships

PATHWAY_COLUMNS = [
    "uniprot_accession",
    "protein_name",
    "symbol",
    "pathway",
    "pathway_name",
    "compound",
    "cid",
    "taxid",
    "taxname",
]

def empty_pathway_df():
    """Returns an empty DataFrame with the columns of the Pathway DataFrame"""
    return pd.DataFrame(columns=PATHWAY_COLUMNS)

def retrieve_pathways(compound, rows, compound_name, selected_tax_ids=None):
    """Checks whether pathway rows exist, and if they do, retrieves the proteins from each pathway."""
    if not rows:
        return empty_pathway_df()
    
    # Send rows that contain Pathway table JSON information 
    df_proteinspathway = retrieve_proteins_from_pathway(compound, rows, compound_name, selected_tax_ids=selected_tax_ids)

    if df_proteinspathway is None or df_proteinspathway.empty:
        return empty_pathway_df()
    
    return df_proteinspathway


def retrieve_proteins_from_pathway(compound, rows, compound_name, selected_tax_ids=None):
    """Given the rows of Pathway section, retrieves the proteins in each pathway and returns a df with the information of the proteins, pathways and compound"""
    
    dfs = []

    for row in rows:
        if not isinstance(row, dict):
            continue

        pathway_id = row.get("pathwayid") or ""
        pwacc = row.get("pwacc") or ""
        pathway_name = row.get("name") or ""

        if not pathway_id or not pwacc:
            continue

        # The table we really want, give me the protein table for this pathway ID
        url_proteins = pcget_pathway_protein_url(pathway_id)
        proteins_table_json = app.utils.get_json(url_proteins)
        
        if proteins_table_json is None:
            continue
    
        # If protein table is downloaded, convert it into a clean DataFrame
        compound_cid = getattr(compound, "cid", None)

        df_targetlist = retrieve_pathway_proteins(pwacc, pathway_name, compound_name, compound_cid, proteins_table_json, selected_tax_ids=selected_tax_ids)

        if df_targetlist is not None and not df_targetlist.empty:
            dfs.append(df_targetlist)
                  
    if not dfs:
        return empty_pathway_df()
    
    df = pd.concat(dfs, ignore_index=True)
    df = df.drop_duplicates(subset=["uniprot_accession", "protein_name", "symbol", "pathway", "compound", "cid", "taxid", "taxname"], keep="first").reset_index(drop=True)
    return df



def pcget_pathway_protein_url(pathwayid, start=1, limit = 10000000):
    """Creates the URL to get the Protein Section from a specific Pathway"""
    return(f"{app.utils.URL_BASE}/assay/pcget.cgi?task=pathway_protein&pathwayid={pathwayid}&start={start}&limit={limit}&infmt=json&outfmt=json")


def retrieve_pathway_proteins(pwacc, pathway_name, compound_name, compound_cid, pathway_json, selected_tax_ids=None):
    """Retrieves proteins from a pathway JSON response"""
    # Retrieve the Proteins in pathway json
    if pathway_json is None:
        return empty_pathway_df()
  
    # Find protein names & id
    rows = pathway_json.get("SDQOutputSet", [{}])[0].get("rows", []) or {}
    status = pathway_json.get("SDQOutputSet", [{}])[0].get("status", {}) or {}
    status_code = status.get("code", 0)

    if str(status_code) != "0":
        return empty_pathway_df()
    
    if not isinstance(rows, list):
        rows = []

    # Extracts protein information one by one
    target_list = []

    for row in rows:
        if isinstance(row,dict): 
            acc_id = row.get("acc") or ""
            protname = row.get("protname")  or ""
            if acc_id:
                target_list.append({"uniprot_accession": acc_id, 
                                       "protein_name": protname,
                                       "pathway": pwacc,
                                       "pathway_name": pathway_name,
                                       "compound": compound_name,
                                       "cid": compound_cid})
    df = pd.DataFrame(target_list)

    if df.empty:
        return empty_pathway_df()
    
    accessions = (df["uniprot_accession"].dropna().astype(str).str.strip().unique().tolist())
    accessions = [acc for acc in accessions if acc]

    df_uniprot_info = pd.DataFrame(
        columns=[
            "uniprot_accession",
            "protein_name",
            "symbol",
            "taxid",
            "taxname",
        ]
    )

    if accessions:
        df_uniprot_info = app.proteins.map_uniprot_to_info(accessions)
        df_uniprot_info = df_uniprot_info.rename(
            columns={
                "protein_name": "uniprot_protein_name",
                "mapped_symbol": "symbol",
            }
        )
        df = df.merge(df_uniprot_info, on= "uniprot_accession", how = "left")

    else:
        df["uniprot_protein_name"] = ""
        df["symbol"] = ""
        df["taxid"] = ""
        df["taxname"] = ""

    df["protein_name"] = df.apply(
        lambda row: str(row.get("protein_name", "")).strip()
        if pd.notna(row.get("protein_name")) and str(row.get("protein_name")).strip()
        else str(row.get("uniprot_protein_name", "")).strip(),
        axis = 1
    )

    for col in ["symbol", "taxid", "taxname"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str).str.strip()
    
    df = df.drop(columns=["uniprot_protein_name"], errors = "ignore")


    selected_tax_ids = app.interactions.normalize_taxonomy_ids(selected_tax_ids)
    if selected_tax_ids:
        df = df[df["taxid"].astype(str).str.strip().isin(selected_tax_ids)].copy()
    

    return df[
        [
            "uniprot_accession", "protein_name", "symbol", "pathway", "pathway_name", "compound", "cid", "taxid", "taxname"]].reset_index(drop=True)

def group_pathways(df_pathways):
    """Groups pathway proteins so each pathway is only present once"""
    
    if df_pathways is None or df_pathways.empty:
        return empty_pathway_df()
    df = df_pathways.copy()

    for col in ["pathway", "pathway_name","protein_name", "uniprot_accession", "compound", "cid", "taxid", "taxname"]:
        if col in df.columns:
            df[col] = df[col].astype("string").fillna("").str.strip()

    df = df[df["pathway"] != ""].copy()
    
    grouped = (
        df.groupby(["pathway", "pathway_name"], as_index=False)
        .agg(
            n_proteins=("uniprot_accession", lambda x: x.replace("", pd.NA).dropna().nunique()),
            n_compounds=("compound", lambda x: x.replace("", pd.NA).dropna().nunique()),
            proteins =("protein_name", lambda x: ";".join(sorted(set(v for v in x if v)))),
            compounds=("compound", lambda x: ";".join(sorted(set(v for v in x if v)))),
            uniprot_accessions=("uniprot_accession", lambda x: ";".join(sorted(set(v for v in x if v)))),
            taxid=("taxid", lambda x: ";".join(sorted(set(v for v in x if v)))),
            taxname=("taxname", lambda x: ";".join(sorted(set(v for v in x if v)))),
        ).sort_values(["n_compounds", "pathway", "pathway_name"], ascending=[False, True, True]).reset_index(drop=True)
    )

    return grouped[["pathway", "pathway_name", "n_proteins", "n_compounds", "proteins", "compounds", "uniprot_accessions", "taxid", "taxname"]]


def group_compounds(df_pathways, selected_pathway = None):
    """Groups pathway proteins so each pathway is only present once"""

    if df_pathways is None or df_pathways.empty:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "count", "compounds", "taxid", "taxname"])
    
    df_pathwaysproteins = df_pathways.copy()

    for col in ["pathway", "pathway_name", "uniprot_accession", "protein_name", "symbol", "compound", "taxid", "taxname"]:
        if col in df_pathwaysproteins.columns:
            df_pathwaysproteins[col] = df_pathwaysproteins[col].fillna("").astype(str).str.strip()

    if selected_pathway is not None:
        selected_pathway = str(selected_pathway).strip()
        df_selected = df_pathwaysproteins[df_pathwaysproteins["pathway"] ==selected_pathway].copy()
    else:
        df_selected = df_pathwaysproteins.copy()

    df_uniprot = df_selected[df_selected["uniprot_accession"] != ""].copy()

    if df_uniprot.empty:
        return pd.DataFrame(columns=["uniprot_accession", "protein_name", "symbol", "count", "compounds", "taxid", "taxname"])

    grouped = (
        df_uniprot.groupby(["uniprot_accession", "protein_name", "symbol"], as_index=False)
        .agg(
            count=("compound", "size"),
            compounds=("compound", lambda x: ";".join(sorted(set(v for v in x if v)))),
            taxid=("taxid", lambda x: ";".join(sorted(set(v for v in x if v)))),
            taxname=("taxname", lambda x: ";".join(sorted(set(v for v in x if v)))),
        )
        .sort_values(["count", "uniprot_accession", "protein_name"], ascending=[False, True, True])
        .reset_index(drop=True)
    )
    return grouped[["uniprot_accession", "protein_name", "symbol", "count", "compounds", "taxid", "taxname"]]
