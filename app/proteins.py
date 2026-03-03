from pathlib import Path
import json
import app.utils
import re
import mygene
import pandas as pd
import pubchempy as pcp
from Bio import Entrez

# -- MANUAL DICTIONARY FOR PROTEIN-GENES
"""OVERRIDES = {
    "m1 receptor": "CHRM1",
    "m2 receptor": "CHRM2",
    "m3 receptor": "CHRM3",
    "m4 receptor": "CHRM4",
    "m5 receptor": "CHRM5",
    "muscarinic acetylcholine receptor m1": "CHRM1",
    "muscarinic acetylcholine receptor m2": "CHRM2",
    "muscarinic acetylcholine receptor m3": "CHRM3",
    "muscarinic acetylcholine receptor m4": "CHRM4",
    "muscarinic acetylcholine receptor m5": "CHRM5",
    "adenosine receptor a1": "ADORA1",
    "a1 receptor": "ADORA1",
    "a2a receptor": "ADORA2A",
    "a2b receptor": "ADORA2B",
    "a3 receptor": "ADORA3",
    "net": "SLC6A2",  # norepinephrine transporter
    "norepinephrine transporter": "SLC6A2",
}
"""


def get_entrez_gene_summary(
    gene_name, email, organism="human", max_gene_ids=100
):
    """Returns the 'Summary' contents for provided input
    gene from the Entrez Gene database. All gene IDs 
    returned for input gene_name will have their docsum
    summaries 'fetched'.
    
    Args:
        gene_name (string): Official (HGNC) gene name 
           (e.g., 'KAT2A')
        email (string): Required email for making requests
        organism (string, optional): defaults to human. 
           Filters results only to match organism. Set to None
           to return all organism unfiltered.
        max_gene_ids (int, optional): Sets the number of Gene
           ID results to return (absolute max allowed is 10K).
        
    Returns:
        dict: Summaries for all gene IDs associated with 
           gene_name (where: keys → [orgn][gene name],
                      values → gene summary)
    """
    Entrez.email = email

    query = (
        f"{gene_name}[Gene Name]"
        if not organism
        else f"({gene_name}[Gene Name]) AND {organism}[Organism]"
    )
    handle = Entrez.esearch(db="gene", term=query, retmax=max_gene_ids)
    record = Entrez.read(handle)
    handle.close()

    gene_summaries = defaultdict(dict)
    gene_ids = record["IdList"]

    print(
        f"{len(gene_ids)} gene IDs returned associated with gene {gene_name}."
    )
    for gene_id in gene_ids:
        print(f"\tRetrieving summary for {gene_id}...")
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="docsum")
        gene_dict = xmltodict.parse(
            "".join([x.decode(encoding="utf-8") for x in handle.readlines()]),
            dict_constructor=dict,
        )
        gene_docsum = gene_dict["eSummaryResult"]["DocumentSummarySet"][
            "DocumentSummary"
        ]
        name = gene_docsum.get("Name")
        summary = gene_docsum.get("Summary")
        gene_organism = gene_docsum.get("Organism")["CommonName"]
        gene_summaries[gene_organism][name] = summary
        handle.close()
        time.sleep(0.34)  # Requests to NCBI are rate limited to 3 per second

    return gene_summaries



def translate_geneid_to_protein(email, protein_data, compound):
    Entrez.email = email
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
    

    for geneid in lines:
        handle = Entrez.esummary(db="gene", id = geneid, retmode="xml")
        records = Entrez.read(handle)
        print(handle)
        print(records)
        handle.close()
        for rec in records:
            gid = str(rec["Id"])
            print(gid)
            symbol = str(rec.get("NomenclatureSymbol") or rec.get("Name") or ""),
            description = str(rec.get("Description" or ""))
            
            protein_list.append({
                "compound": compound,
                "geneid": gid,
                "description": description,
                "symbol": symbol,
            })

    out = pd.DataFrame(protein_list)
    print(out)
    return out



def translate_geneid_to_proteins(email: str, protein_data: str, compound: str, sleep_s: float = 0.34):
    """
    Reads GeneIDs from {compound}/{protein_data} and maps them via NCBI Gene esummary.
    Returns a DataFrame: compound, geneid, symbol, description
    """
    Entrez.email = email

    folder = compound
    path = Path(folder) / protein_data

    if not path.exists():
        print(f"[skip] {compound}: no {protein_data}")
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])

    # read file lines
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    gene_ids = [str(x).strip() for x in lines if str(x).strip().isdigit()]

    if not gene_ids:
        print(f"[warn] {compound}: no numeric GeneIDs found in {protein_data}")
        return pd.DataFrame(columns=["compound", "geneid", "symbol", "description"])

    # chunk requests
    chunk_size = 200
    rows = []

    for i in range(0, len(gene_ids), chunk_size):
        chunk = gene_ids[i:i+chunk_size]

        handle = Entrez.esummary(db="gene", id=",".join(chunk), retmode="xml")
        data = Entrez.read(handle)
        handle.close()

        docs = data["DocumentSummarySet"]["DocumentSummary"]
        for rec in docs:
            gid = str(rec["Id"])
            symbol = str(rec.get("NomenclatureSymbol") or rec.get("Name") or "")
            desc = str(rec.get("Description") or "")
            rows.append({
                "compound": compound,
                "geneid": gid,
                "symbol": symbol,
                "description": desc,
            })

        time.sleep(sleep_s)

    return pd.DataFrame(rows)

def translate_geneid_to_protein(protein_data: str, compound: str):
    folder = compound
    path = Path(folder) / protein_data

    if not path.exists():
        print(f"[skip] {compound}: no {protein_data}")
        return None

    gene_ids = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    mapping = geneids_to_symbols(gene_ids)

    # build rows
    rows = []
    for gid in gene_ids:
        gid = str(gid).strip()
        if not gid:
            continue
        rec = mapping.get(gid, {})
        rows.append({
            "compound": compound,
            "geneid": gid,
            "symbol": rec.get("symbol"),
            "description": rec.get("description"),
            "name": rec.get("name"),
        })

    return rows

Entrez.email = "your_email@example.com"

def geneid_to_protein_ids(gene_id: str) -> list[str]:
    """
    Returns linked protein record IDs from NCBI protein database.
    """
    h = Entrez.elink(dbfrom="gene", db="protein", id=str(gene_id))
    linkset = Entrez.read(h)
    h.close()

    protein_ids = []
    for ls in linkset:
        for linkdb in ls.get("LinkSetDb", []):
            if linkdb.get("LinkName") == "gene_protein":
                protein_ids.extend([l["Id"] for l in linkdb.get("Link", [])])
    return protein_ids

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
# ask what this does

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
        df_mapped["canonical_target"] = df_mapped["canonical_target"].replace(OVERRIDES)
        mask = (~df_mapped["is_gene"] & df_mapped["gene_symbol_mygene"].notna())
        df_mapped.loc[mask, "canonical_target"] = df_mapped.loc[mask, "gene_symbol_mygene"]
    return df_mapped, df_map


#omg what does this do
def results_summary_count(df_mapped):
    summary = (
        df_mapped.groupby("canonical_target")
            .agg(
                total_count=("canonical_target", "size"),
                n_compounds=("compound", "nunique"),
                compounds=("compound", lambda x: ";".join(sorted(set(x)))),
            )
            .reset_index()
            .sort_values(["total_count", "n_compounds"], ascending=[False, False])
    )

    return summary
