from pathlib import Path
import re
import pandas as pd
from Bio import Entrez
import requests
import time


import app.utils

UNIPROT_URL = "https://rest.uniprot.org"

def chunk_list(items,size):
    for i in range(0,len(items), size):
        yield items[i:i + size]

# Appear in main so that the user can fetch the proteins in the interactions of the compound?
def retrieve_targets_1(compound):
    compound_name = app.utils.safe_filename(compound.synonyms[0] if compound.synonyms else compound.cid)
    
    # Retrieve the Interactions json
    rows = app.utils.load_json(f"compound_{compound.cid}_{compound_name}_interactionstable.json", "1.Chemical-Target Interactions")
    if rows is None:
        return None
    
    target_list = [] 
    for row in rows:
        if isinstance(row,dict): # is row a dictionary object? -> rows should be a list of dictionaries (.get() only exists in dictionaries)
            gene_id = row.get("geneid")
            if gene_id:
                target_list.append(str(gene_id).strip())

    
    #app.utils.create_text_file("protein_data", compound_name)
    #app.utils.write_text_file("protein_data", compound_name, target_list)
   
    df_geneids = pd.DataFrame(
        {
            "compound": [compound_name] * len(target_list),
            "geneid": target_list,
        },
        columns = ["compound", "geneid"]
    )

    return df_geneids

# -- MY METHOD
def translate_geneid_to_protein(email, df_geneids, compound, api_key=None, batch_size=200):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Retries
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 20

    compound_name = app.utils.safe_filename(compound.synonyms[0] if compound.synonyms else compound.cid)

    if df_geneids is None or df_geneids.empty or "geneid" not in df_geneids.columns:
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])
    
    lines = (df_geneids["geneid"].dropna().astype(str).str.strip())
    lines = [x for x in lines.unique().tolist() if x] #removes blank strings and keeps only unique geneids to a Python list

    if not lines:
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])
    
    protein_list = []

    for batch in chunk_list(lines, batch_size):
        ids= ",".join(batch)
        try:
            handle = Entrez.esummary(db="gene", id = ids, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
        except Exception as e:
            continue

        # In Entrez JSON this is how you access the databases??? TODO
        documents = records["DocumentSummarySet"]["DocumentSummary"]

        for rec in documents:
            gid = str(rec.attributes.get("uid", ""))
            symbol = str(rec.get("NomenclatureSymbol") or rec.get("Name") or "").upper()
            description = str(rec.get("Description") or "")
            
            protein_list.append({
                "compound": compound_name,
                "geneid": gid,
                "description": description,
                "symbol": symbol,
            })
            #small pause
            time.sleep(0.2)

    return pd.DataFrame(protein_list,
                        columns=["compound", "geneid","symbol","description"],
                        )
 

# -- MAP GENEID TO UNIPROTKB ACCESSION CODES (1128 -> P11229)
def get_idmapping_db(db_name):
# Look in Uniprot ID-mapping database code GeneID and UniProtKB
    req = requests.get(f"{UNIPROT_URL}/configure/idmapping/fields", timeout=30).json()
    
    for group in req.get("groups", []):
        for item in group.get("items", []):
            if item.get("displayName") == db_name:
                return item["name"]
            
    raise ValueError(f"Could not find mapping database with displayName:{db_name!r}") #TODO what does the r! mean?


def input_idmapping_dbs(from_db, to_db, gene_list):
# Tell UniProt to do the translation process with both GeneID db and UniProtKB db with the gene_ids list
    data = {"from": from_db, "to":to_db, "ids": ",".join(map(str,gene_list))}
    req = requests.post(f"{UNIPROT_URL}/idmapping/run", data = data, timeout=60)
    req.raise_for_status()
    return req.json()["jobId"] # not the result but the response of the process of translation


def wait_for_job(jobId, repeats = 2):
# Loop with sleep + timeout until the job is FINISHED 
    start_time = time.time()
    timeout=450
    while True:
        req = requests.get(f"{UNIPROT_URL}/idmapping/status/{jobId}", timeout=30, allow_redirects=False)
        req.raise_for_status()
        
        if req.status_code == 303:
            print("Job finished!")
            return

        data = req.json()
        status = data.get("jobStatus") or data.get("status") #TODO: check this field
        
        if status in ("FINISHED", "DONE"):
            return
        if status in ("FAILED", "ERROR"):
            raise RuntimeError(f"UniProt mapping job failed: {data}")
        if time.time() - start_time > timeout:
            raise TimeoutError(f"UniProt mapping job timed out after {timeout} seconds")
        
        time.sleep(repeats)


# TODO: que es el codigo del link
def download_results(jobId):
# Retrieve the JSON results and download them
    url = f"{UNIPROT_URL}/idmapping/results/{jobId}"
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


def map_genes_to_uniprot(df_geneids):
    if df_geneids is None or df_geneids.empty or "geneid" not in df_geneids.columns:
        return pd.DataFrame(columns=["geneid", "uniprot_accession", "uniprot_accessions"])
    
    geneids = (df_geneids["geneid"].dropna().astype(str).str.strip())
    geneids = [x for x in geneids.unique().tolist() if x] #removes blank strings and keeps only unique geneids to a Python list
    
    if not geneids:
        return pd.DataFrame(columns=["geneid", "uniprot_accession", "uniprot_accessions"])

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

    return pd.DataFrame(rows,
                        columns=["geneid", "uniprot_accession", "uniprot_accessions"],
                        )

def map_uniprot_to_symbol(accessions):
    accessions = [str(x).strip() for x in accessions if pd.notna(x) and str(x).strip()]

    if not accessions:
        return pd.DataFrame(columns=["uniprot_accession", "mapped_symbol"])
    
    rows = []

    with requests.Session() as session:
        for chunk in chunk_list(accessions, 100):
            query = " OR ".join(f"accession:{acc}" for acc in chunk)

            response = session.get(f"{UNIPROT_URL}/uniprotkb/search",
                                   params = {"query": query,
                                             "fields": "accession,gene_primary",
                                             "format": "json",
                                             "size": len(chunk),
                                             },
                                             timeout = 60,
                                        )
            response.raise_for_status()
            data = response.json()

            for item in data.get("results", []):
                accession = item.get("primaryAccession", "") or ""

                mapped_symbol = ""
                genes = item.get("genes", []) or []
                if genes:
                    gene_name = genes[0].get("geneName") or {}
                    mapped_symbol = gene_name.get("value", "") or ""

                rows.append({
                    "uniprot_accession": accession,
                    "mapped_symbol": mapped_symbol,
                })
    return pd.DataFrame(rows, columns=["uniprot_accession", "mapped_symbol"])
            

# ESTE SE UTILIZA??
def results_proteins(df_interactions, df_pathways):
    df = df_interactions.copy()
    df["compound"] = df["compound"].astype(str).str.strip()
    df["geneid"] = df["geneid"].astype(str).str.strip()
    df["symbol"] = df["symbol"].astype(str).str.strip()
    df["uniprot_accession"] = df["uniprot_accession"].astype(str).str.strip()

    df = df[df["uniprot_accession"].notna()]
    df = df[df["uniprot_accession"] != ""]
    df = df[df["uniprot_accession"].str.lower() != "nan"]
    
   
    summary =(
        df.groupby("uniprot_accession").agg(
            total_count = ("uniprot_accession", "size"),
            n_compounds = ("compound", "nunique"),
            compounds = ("compound", lambda x: ";".join(sorted(set(x)))),
            symbol = ("symbol", lambda x: ";".join(sorted(set(x)))),
            geneid = ("geneid", lambda x: ";".join(sorted(set(x)))),
        ).reset_index().sort_values(["total_count"], ascending = [False])
    )

    return summary

# Se podria elegir que aspecto quieres buscar ({biological_process, molecular_function, cellular_component})
def fetch_goterms(df, aspects = None):
    if df is None or df.empty or "uniprot_accession" not in df.columns:
        return pd.DataFrame(columns=["uniprot_accession", "go_id", "aspect"])
    
    accessions = df["uniprot_accession"].dropna().astype(str).str.strip()
    accessions = accessions[accessions != ""].unique().tolist() #removes blank strings and keeps only unique accessions to a Python list
    
    if aspects is None:
        aspects = [
            "biological_process",
            "molecular_function",
            "cellular_component",
        ]

    check_aspects = {
        "biological_process",
        "molecular_function",
        "cellular_component",
    }

    aspects = [a.strip() for a in aspects if str(a).strip()]
    invalid = [a for a in aspects if a not in check_aspects]
    if invalid:
        raise ValueError(
            "Invalid aspects"
            f"Allowed aspects are: {check_aspects}"
        )


    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    headers = {"Accept": "application/json"}

    
    rows = []
    limit = 200
    
    # I have to access more than one page (pagination)
    with requests.Session() as session:
        for acc in accessions:
            page = 1

            # For pages to add each loop
            while True:
                parameters = {
                    "geneProductId": f"UniProtKB:{acc}",
                    "aspect": ",".join(aspects),
                    "limit":200,
                    "page": page,
                }
                request = session.get(url, params=parameters, headers=headers, timeout=60)
                request.raise_for_status()
                result_json = request.json()

                results = result_json.get("results", [])
                if not results:
                    break

                for row in results:
                    rows.append({
                        "uniprot_accession": acc,
                        "go_id": row.get("goId"),
                        "aspect": row.get("goAspect"),
                    })
                if len(results) < limit:
                    break
                page +=1

    return pd.DataFrame(rows,
                        columns=["uniprot_accession", "go_id", "aspect"],
                        )


def fetch_gonames(go_ids):
    go_ids = sorted({str(x).strip() for x in go_ids if pd.notna(x) and str(x).strip()})
    if not go_ids:
        return pd.DataFrame(columns=["go_id", "go_name"])
    
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{ids}"
    headers = {"Accept": "application/json"}

    rows = []
    # there are thousands of go_ids, so send them 100 by 100
    chunk_size = 100

    with requests.Session() as session:
        for i in range(0, len(go_ids), chunk_size):
            chunk = go_ids[i:i+chunk_size]
            ids_str = ",".join(chunk)

            request = session.get(url.format(ids=ids_str), headers=headers, timeout=60)
            request.raise_for_status()
            result_json = request.json()

            for row in result_json.get("results", []):
                rows.append({
                    "go_id": row.get("id"),
                    "go_name": row.get("name"),
                })
    
    return pd.DataFrame(rows,
                        columns=["go_id", "go_name"],
                        )

def summarize_goaspect(df_go, aspect, prefix):
    df_aspect = df_go[df_go["aspect"] == aspect].copy()

    if df_aspect.empty:
        return pd.DataFrame(columns=["uniprot_accession", f"go_{prefix}_ids", f"go_{prefix}_names",])
    
    return (
        df_aspect.groupby("uniprot_accession", as_index=False).agg(
            **{
                f"go_{prefix}_ids": ("go_id", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
                f"go_{prefix}_names": ("go_name", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            }
        )
    )