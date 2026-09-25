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


def retrieve_proteins_from_pathway(
    compound,
    rows,
    compound_name,
    selected_tax_ids=None,
):
    """Given pathway rows, retrieves the proteins associated with each pathway."""

    dfs = []

    for row in rows:
        if not isinstance(row, dict):
            continue

        raw_pwacc = row.get("pwacc") or ""
        pwacc = str(raw_pwacc).replace("\\:", ":").strip()
        pathway_name = row.get("name") or ""
        geneids = row.get("geneids") or ""

        print("[Pathway DEBUG] pwacc:", repr(pwacc))
        print("[Pathway DEBUG] pathway:", pathway_name)
        print("[Pathway DEBUG] geneids:", repr(geneids))

        if not pwacc:
            continue

        # ---------------------------------------------------------
        # Reactome / PharmGKB:
        # retrieve protein accessions directly from the pathway
        # ---------------------------------------------------------
        if pwacc.startswith(("Reactome:", "PharmGKB:")):

            url_proteins = pcget_pathway_protein_url(pwacc)

            print("[Pathway DEBUG] protein URL:", url_proteins)

            pathway_json = app.utils.get_json(url_proteins)

            print(
                "[Pathway DEBUG] protein response type:",
                type(pathway_json),
            )

            if pathway_json is None:
                print("[Pathway DEBUG] protein response is None")
                continue

            df_targetlist = retrieve_pathway_proteins(
                pwacc,
                pathway_name,
                compound_name,
                getattr(compound, "cid", None),
                pathway_json,
                selected_tax_ids=selected_tax_ids,
            )

        # ---------------------------------------------------------
        # WikiPathways:
        # retrieve Gene IDs from the pathway
        # ---------------------------------------------------------
        elif pwacc.startswith("WikiPathways:"):

            url_genes = pcget_pathway_gene_url(pwacc)

            print("[Pathway DEBUG] gene URL:", url_genes)

            pathway_json = app.utils.get_json(url_genes)

            print(
                "[Pathway DEBUG] gene response type:",
                type(pathway_json),
            )

            if pathway_json is None:
                print("[Pathway DEBUG] gene response is None")
                continue

            df_targetlist = retrieve_pathway_proteins_from_genes(
                pwacc,
                pathway_name,
                compound_name,
                getattr(compound, "cid", None),
                pathway_json,
                selected_tax_ids=selected_tax_ids,
            )

        else:
            print(
                "[Pathway DEBUG] Unsupported pathway source:",
                pwacc.split(":", 1)[0],
            )
            continue

        if df_targetlist is not None and not df_targetlist.empty:
            dfs.append(df_targetlist)

    if not dfs:
        return empty_pathway_df()

    df = pd.concat(dfs, ignore_index=True)

    df = df.drop_duplicates(
        subset=[
            "uniprot_accession",
            "protein_name",
            "symbol",
            "pathway",
            "compound",
            "cid",
            "taxid",
            "taxname",
        ],
        keep="first",
    ).reset_index(drop=True)

    return df

def pcget_pathway_gene_url(pwacc):
    """Creates the PubChem PUG REST URL to retrieve Gene IDs for a pathway."""

    return (
        f"{app.utils.URL_BASE}/rest/pug/pathway/pwacc/"
        f"{pwacc}/geneids/JSON"
    )


def pcget_pathway_protein_url(pwacc):
    """Creates the PubChem PUG REST URL to retrieve protein accessions for a pathway."""

    return (
        f"{app.utils.URL_BASE}/rest/pug/pathway/pwacc/"
        f"{pwacc}/accessions/JSON"
    )


def build_pathway_protein_dataframe(
    pwacc,
    pathway_name,
    compound_name,
    compound_cid,
    accessions,
    selected_tax_ids=None,
):
    """Build pathway-protein DataFrame and enrich proteins with UniProt metadata."""

    if not accessions:
        return empty_pathway_df()

    # Create pathway-protein rows
    target_list = [
        {
            "uniprot_accession": acc,
            "protein_name": "",
            "pathway": pwacc,
            "pathway_name": pathway_name,
            "compound": compound_name,
            "cid": compound_cid,
        }
        for acc in accessions
    ]

    df = pd.DataFrame(target_list)

    # Retrieve UniProt metadata
    df_uniprot_info = app.proteins.map_uniprot_to_info(accessions)

    if df_uniprot_info is not None and not df_uniprot_info.empty:

        df_uniprot_info = df_uniprot_info.rename(
            columns={
                "protein_name": "uniprot_protein_name",
                "mapped_symbol": "symbol",
            }
        )

        df = df.merge(
            df_uniprot_info,
            on="uniprot_accession",
            how="left",
        )

    else:

        df["uniprot_protein_name"] = ""
        df["symbol"] = ""
        df["taxid"] = ""
        df["taxname"] = ""

    # Prefer the original protein name if available,
    # otherwise use the UniProt protein name.
    df["protein_name"] = df.apply(
        lambda row:
            str(row.get("protein_name", "")).strip()
            if (
                pd.notna(row.get("protein_name"))
                and str(row.get("protein_name")).strip()
            )
            else str(row.get("uniprot_protein_name", "")).strip(),
        axis=1,
    )

    for col in ["symbol", "taxid", "taxname"]:

        if col not in df.columns:
            df[col] = ""

        df[col] = (
            df[col]
            .fillna("")
            .astype(str)
            .str.strip()
        )

    df = df.drop(
        columns=["uniprot_protein_name"],
        errors="ignore",
    )

    # Apply taxonomy filter
    selected_tax_ids = app.interactions.normalize_taxonomy_ids(
        selected_tax_ids
    )

    if selected_tax_ids:
        df = df[
            df["taxid"]
            .astype(str)
            .str.strip()
            .isin(selected_tax_ids)
        ].copy()

    return df[
        [
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
    ].reset_index(drop=True)


def retrieve_pathway_proteins_from_genes(
    pwacc,
    pathway_name,
    compound_name,
    compound_cid,
    pathway_json,
    selected_tax_ids=None,
):
    """Retrieve UniProt proteins associated with WikiPathways Gene IDs."""

    if pathway_json is None:
        return empty_pathway_df()

    geneids = []

    if isinstance(pathway_json, dict):

        information_list = pathway_json.get("InformationList", {})

        if isinstance(information_list, dict):

            information = information_list.get("Information", [])

            if isinstance(information, list):

                for item in information:

                    if not isinstance(item, dict):
                        continue

                    ids = item.get("GeneID", [])

                    if isinstance(ids, list):
                        geneids.extend(ids)

                    elif ids is not None:
                        geneids.append(ids)

    geneids = (
        pd.Series(geneids, dtype="string")
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda x: x != ""]
        .unique()
        .tolist()
    )

    print(
        "[Pathway DEBUG] Gene IDs found:",
        len(geneids),
    )

    print(
        "[Pathway DEBUG] Gene IDs:",
        geneids[:20],
    )

    if not geneids:
        print("[Pathway DEBUG] NO GENE IDS FOUND")
        return empty_pathway_df()

    # -------------------------------------------------------------
    # Map NCBI Gene IDs to UniProt accessions
    # -------------------------------------------------------------
    df_geneids = pd.DataFrame({"geneid": geneids})

    df_mapped = app.proteins.map_genes_to_uniprot(df_geneids)

    print(
        "[Pathway DEBUG] GeneID → UniProt mapping rows:",
        len(df_mapped),
    )

    if df_mapped is None or df_mapped.empty:
        print("[Pathway DEBUG] NO UNIPROT ACCESSIONS FOUND FOR GENE IDS")
        return empty_pathway_df()

    accessions = (
        df_mapped["uniprot_accession"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda x: x != ""]
        .unique()
        .tolist()
    )

    print(
        "[Pathway DEBUG] UniProt accessions from Gene IDs:",
        len(accessions),
    )
    print("[Pathway DEBUG] accessions:", accessions[:20])

    if not accessions:
        print("[Pathway DEBUG] NO UNIPROT ACCESSIONS FOUND FOR GENE IDS")
        return empty_pathway_df()

    print(
        "[Pathway DEBUG] UniProt accessions from Gene IDs:",
        len(accessions),
    )

    print(
        "[Pathway DEBUG] accessions:",
        accessions[:20],
    )

    if not accessions:
        print(
            "[Pathway DEBUG] NO UNIPROT ACCESSIONS FOUND FOR GENE IDS"
        )
        return empty_pathway_df()

    return build_pathway_protein_dataframe(
        pwacc,
        pathway_name,
        compound_name,
        compound_cid,
        accessions,
        selected_tax_ids,
    )


def retrieve_pathway_proteins(
    pwacc,
    pathway_name,
    compound_name,
    compound_cid,
    pathway_json,
    selected_tax_ids=None,
):
    """Extract UniProt protein accessions from a PubChem pathway response."""

    if pathway_json is None:
        return empty_pathway_df()

    accessions = []

    if isinstance(pathway_json, dict):

        information_list = pathway_json.get("InformationList", {})

        if isinstance(information_list, dict):

            information = information_list.get("Information", [])

            if isinstance(information, list):

                for item in information:

                    if not isinstance(item, dict):
                        continue

                    accs = item.get("ProteinAccession", [])

                    if isinstance(accs, list):
                        accessions.extend(accs)

                    elif isinstance(accs, str):
                        accessions.append(accs)

    accessions = (
        pd.Series(accessions, dtype="string")
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda x: x != ""]
        .unique()
        .tolist()
    )

    print(
        "[Pathway DEBUG] UniProt accessions found:",
        len(accessions),
    )

    print(
        "[Pathway DEBUG] accessions:",
        accessions[:20],
    )

    if not accessions:
        print("[Pathway DEBUG] NO PROTEIN ACCESSIONS FOUND")
        return empty_pathway_df()

    return build_pathway_protein_dataframe(
        pwacc,
        pathway_name,
        compound_name,
        compound_cid,
        accessions,
        selected_tax_ids,
    )


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
