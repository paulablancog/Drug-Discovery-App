from pathlib import Path
import json
import app.utils
import re
import mygene
import pandas as pd
import pubchempy as pcp
from Bio import Entrez
import requests
import time


def chunk_list(items,size):
    for i in range(0,len(items), size):
        yield items[i:i + size]

# -- MY METHOD
def translate_geneid_to_protein(email, protein_data, compound, api_key=None, batch_size=200):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Retries
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 20

    folder = compound.synonyms[0]
    path = Path(folder) / protein_data

    if not path.exists():
        print("Skipping compound, no protein (Weird, check interactions)")
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])

    lines = app.utils.read_file_lines(protein_data, folder)
    protein_list = []

    if not lines:
        print(f"[warn] {compound}: no numeric GeneIDs found in {protein_data}")
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])
    
    for batch in chunk_list(lines, batch_size):
        ids= ",".join(batch)
        try:
            handle = Entrez.esummary(db="gene", id = ids, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print("Failed batching")
            continue

        documents = records["DocumentSummarySet"]["DocumentSummary"]

        for rec in documents:
            gid = str(rec.attributes.get("uid", ""))
            #print(gid)
            symbol = str(rec.get("NomenclatureSymbol") or rec.get("Name") or "").upper()
            description = str(rec.get("Description") or "")
            
            protein_list.append({
                "compound": compound.synonyms[0],
                "geneid": gid,
                "description": description,
                "symbol": symbol,
            })
            #small pause
            time.sleep(0.2)

    out = pd.DataFrame(protein_list)
    return out

UniProt_url = "https://rest.uniprot.org"

# -- MAP GENEID TO UNIPROTKB ACCESSION CODES (1128 -> P11229)
def get_idmapping_db(db_name):
# Look in Uniprot ID-mapping database code GeneID and UniProtKB
    req = requests.get(f"{UniProt_url}/configure/idmapping/fields", timeout=30).json()
    print("Printing idmapping dbs")
    
    for group in req.get("groups", []):
        for item in group.get("items", []):
            if item.get("displayName") == db_name:
                return item["name"]
    raise ValueError(f"Could not find mapping database with displayName:{db_name!r}") #TODO what does the r! mean?

def input_idmapping_dbs(from_db, to_db, gene_list):
# Tell UniProt to do the translation process with both GeneID db and UniProtKB db with the gene_ids list
    data = {"from": from_db, "to":to_db, "ids": ",".join(map(str,gene_list))}
    req = requests.post(f"{UniProt_url}/idmapping/run", data = data, timeout=60)
    req.raise_for_status()
    return req.json()["jobId"] # not the result but the response of the process of translation

def wait_for_job(jobId, repeats = 2):
# Loop with sleep + timeout until the job is FINISHED 
    t = time.time()
    timeout=450
    while True:
        req = requests.get(f"{UniProt_url}/idmapping/status/{jobId}", timeout=30, allow_redirects=False)
        req.raise_for_status()
        print(req)
        
        if req.status_code == 303:
            print("Job finished!")
            return

        json = req.json()
        status = json.get("jobStatus") or json.get("status") #TODO: check this field
        
        if status in ("FINISHED", "DONE"):
            return
        if status in ("FAILED", "ERROR"):
            raise RuntimeError(f"UniProt mapping job failed: {json}")
        if time.time() - t > timeout:
            raise TimeoutError(f"UniProt mapping job timed out after {timeout} seconds")
        
        time.sleep(repeats)


# TODO: que es el codigo del link
def download_results(jobId):
# Retrieve the JSON results and download them
    url = f"{UniProt_url}/idmapping/results/{jobId}"
    params = {"format": "json"}
    out = {}

    while url:
        req = requests.get(url, params, timeout=60)
        req.raise_for_status()
        json = req.json()
    # Now let's download the resulting table with the accession code
        for row in json.get("results", []):
            geneid=str(row.get("from")).strip()
            accession = row.get("to")
            if not geneid or not accession:
                continue

            out.setdefault(geneid,[])
            if accession not in out[geneid]:
                out[geneid].append(accession)
        params = None

        link = req.headers.get("Link")
        next_url = None
        if link:
            m = re.search(r'<([^>]+)>;\s*rel="next"', link)
            if m:
                next_url = m.group(1)
        url = next_url
    return out

def P_accession(accessions):
    if not accessions:
        return None
    for a in accessions:
        if str(a).startswith("P"):
            return a
    return accessions[0]

def map_genes_to_uniprot(filename, compound):
    folder = compound.synonyms[0]
    path = Path(folder) / filename

    if not path.exists():
        print("Skipping compound, no GeneID file found")
        return pd.DataFrame(columns=["geneid", "uniprot_accession", "uniprot_accessions"])

    lines = app.utils.read_file_lines(filename, folder)

    if not lines:
        print(f"[warn] {compound}: no numeric GeneIDs found in {filename}")
        return pd.DataFrame(columns=["compound", "uniprot_accession", "uniprot_accessions"])
    
    #TODO mirar esta linea interesante
    geneids = [str(x).strip() for x in lines if str(x).strip()]

    from_db = get_idmapping_db("GeneID")
    to_db = get_idmapping_db("UniProtKB")

    jobId = input_idmapping_dbs(from_db, to_db, geneids)
    wait_for_job(jobId)

    results = download_results(jobId)

    rows = []
    for gene in geneids:
        accessions = results.get(str(gene), [])
        reviewed_acc = P_accession(accessions)
        rows.append({
            "geneid": str(gene),
            "uniprot_accession": reviewed_acc,
            "uniprot_accessions": ";".join(accessions) if accessions else None,
        })
    return pd.DataFrame(rows)


# Appear in main so that the user can fetch the proteins in the interactions of the compound?
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
            gene_id = row.get("geneid")
            target_list.append(gene_id)
    target_clean_list = list(filter(None, target_list))
    print("Successfully retrieved Gene IDs for Chemical-Target interactions")

    # Save protein set in folder
    filename = f"protein_data"
    folder = f"{compound.synonyms[0]}"
    app.utils.create_file(filename, folder)
    app.utils.write_file(filename, folder, target_clean_list)
    
    print(f"Saved proteins of {compound} to {filename}")
    return filename


def normalize_protein(proteintxt, compound):
    # in the main this has to be in a loop so it can reach every protein.txt from each compound
    folder = compound
    gene_symbol = re.compile(r"^[A-Z0-9]{2,15}$")
    lines = app.utils.read_file_lines(proteintxt, folder)
    protein_list = []

    for p in lines:
        p = p.strip()
        if not p:
            continue

        is_gene = bool(gene_symbol.match(p)) and p.upper() == p

        if is_gene:
            normalized = p.upper() # keep the gene symbol annotation (HGNC)

        else: 
            normalized = p
            normalized = re.sub(r"\([^)]*\)", " ", normalized) # replace parenthesis with space
            normalized = normalized.lower()
            normalized = re.sub(r"[^a-z0-9\s]", " ", normalized) #replace punctionation and symbols with spaces
            normalized = re.sub(r"\s+", " ", normalized).strip() # remove whitespaces
        
        protein_list.append({"compound": compound,
                             "protein": p,
                             "normalized_target": normalized,
                             "is_gene": is_gene})
    
    df = pd.DataFrame(protein_list)
    print(df)
    return df

# Now I need to call normalized_protein for each compound, store the returned DataFrames in a list and then concat
# ask what this does -> TODO: ESTE METODO YA NO ME SIRVE (UTILIZO LOS DE UNIPROT_ACCESSION)

def read_target_proteins (targets_pathway):
    df = pd.read_csv(targets_pathway)
    if "protein" not in df.columns or "normalized_target" not in df.columns or "compound" not in df.columns:
        print("Incorrect CSV.")

    df["protein"] = df["protein"].astype(str).str.strip()
    df["normalized_target"] = df["normalized_target"].astype(str).str.strip()

    if "is_gene" not in df.columns:
        print("There is no identifier in the CSV. Check information CSV")

    if df["is_gene"].dtype == object:
        # Converts string into booleans True/False
        df["is_gene"] = df["is_gene"].map({"True":True, "False": False}).fillna(df["is_gene"])
    
    df["is_gene"] = df["is_gene"].astype(bool)

    query = (df.loc[~df["is_gene"], "normalized_target"].dropna().astype(str).unique().tolist())

    mg = mygene.MyGeneInfo()
    results = mg.querymany(query, scopes=["symbol", "name", "alias", "uniprot"], fields=["symbol", "name", "entrezgene", "score"], species="human", returnall=False, as_dataframe = False, verbose = False,)
    
    map_rows = []
    for r in results:
        response = r.get("query")
        if not response:
            continue
        map_rows.append({
                "normalized_target": response,
                "gene_symbol_mygene": (str(r.get("symbol")).upper()
                                       if r.get("symbol") else None),
                "score": r.get("score"),
                "notfound": bool(r.get("notfound", False)),                      
        })

        df_map = pd.DataFrame(map_rows)

        if not df_map.empty and "score" in df_map.columns:
            df_map = df_map.sort_values("score", ascending=False).drop_duplicates("normalized_target", keep="first")

        df_mapped = df.merge(df_map, on="normalized_target", how="left")

        df_mapped["canonical_target"] = df_mapped["normalized_target"]
        df_mapped["canonical_target"] = df_mapped["canonical_target"]
        mask = (~df_mapped["is_gene"] & df_mapped["gene_symbol_mygene"].notna())
        df_mapped.loc[mask, "canonical_target"] = df_mapped.loc[mask, "gene_symbol_mygene"]
    return df_mapped, df_map


#omg what does this do
def results_summary_count(csv_files):
    dfs = []

    for file in csv_files:
        df = pd.read_csv(file)

        df["compound"] = df["compound"].astype(str).str.strip()
        df["geneid"] = df["geneid"].astype(str).str.strip()
        df["symbol"] = df["symbol"].astype(str).str.strip()
        df["uniprot_accession"] = df["uniprot_accession"].astype(str).str.strip()

        df = df[df["uniprot_accession"].notna()]
        df = df[df["uniprot_accession"] != ""]
        df = df[df["uniprot_accession"].str.lower() != "nan"]

        dfs.append(df)
    if not dfs:
        return pd.DataFrame(columns=["uniprot_accession", "total_count", "n_compounds", "compounds", "symbol", "geneid"
        ])
    
    all_df = pd.concat(dfs, ignore_index=True)

    summary =(
        all_df.groupby("uniprot_accession").agg(
            total_count = ("uniprot_accession", "size"),
            n_compounds = ("compound", "nunique"),
            compounds = ("compound", lambda x: ";".join(sorted(set(x)))),
            symbol = ("symbol", lambda x: ";".join(sorted(set(x)))),
            geneid = ("geneid", lambda x: ";".join(sorted(set(x)))),
        ).reset_index().sort_values(["total_count"], ascending = [False])
    )

    return summary
