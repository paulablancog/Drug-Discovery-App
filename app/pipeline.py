from pathlib import Path
import pandas as pd

import app.chem
import app.interactions
import app.proteins
import app.pathways
import app.utils

def fetch_pubchem_compound(smiles_code, email):
    compound = app.chem.compound_retrieval(smiles_code)
    if compound is None:
        raise ValueError("No compound found. Check the SMILES code")
    
    compound_info = app.chem.compound_information(compound)

    load_save_interactions(compound)

    df_geneids = app.proteins.retrieve_targets_1(compound)
    proteins_data = app.proteins.translate_geneid_to_protein(email, df_geneids, compound)
    df_map = app.proteins.map_genes_to_uniprot(df_geneids)

    protein_data = proteins_data.drop_duplicates(subset=["geneid"]) if "geneid" in proteins_data.columns else proteins_data

    if "geneid" not in df_map.columns:
        df_map = pd.DataFrame(columns=["geneid", "uniprot_accession", "uniprot_accessions"])
    
    if protein_data.empty:
        df_proteins = pd.DataFrame(columns=["compound", "geneid", "symbol", "description", "uniprot_accessions"])
    else:
        df_proteins = protein_data.merge(df_map, on="geneid", how="left")

    compound_name = app.utils.safe_filename(app.chem.compound_display_name(compound))
    df_proteins.to_csv(f"{compound_name}_UniProt_proteins.csv", index=False)

    df_pathways = app.pathways.retrieve_pathways(compound)
    if df_pathways is None:
        df_pathways = pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound"])

    df_pathways.to_csv(f"{compound_name}_ProteinsAllPathways.csv", index=False)
    
    return compound, compound_info, df_proteins


def load_save_interactions(compound):
    compound_name = app.utils.safe_filename(app.chem.compound_display_name(compound))

    # 1. Creates the URL
    index_url = f"{app.utils.URL_BASE}/rest/pug_view/index/compound/{compound.cid}/JSON"
    # 2. Retrieves the index JSON
    index_json = app.utils.get_json(index_url)

    if index_json is None:
        raise ValueError("Failed to retrieve index JSON")

    # 3. Saves the index JSON to a file inside the indexes folder
    app.utils.save_json(index_json, f"compound_{compound.cid}_{compound_name}_index.json", "indexes")
    if not app.interactions.has_interactions_and_pathways(index_json):
        raise ValueError("No Interactions and Pathways section found")
    
    # 4. Gets all sections in the index JSON
    # 5. Finds the Interactions and Pathways section and retrieves the data of interactions
    data = app.interactions.load_interactions_and_pathways_data(compound)
    if data is None:
        raise ValueError("No Interactions data for this compound")

    # 6. Get external tables
    tables = app.interactions.retrieve_externaltable(data)
    
    # 7. Download and save each table
    where_pathways = {"ands": [{"cid": str(compound.cid)}, {"core": "1"}]}

    for subsection, table_list in tables:
        subsection_1 = (subsection or "").lower()

        clean_tables = []
        seen = set()
        for table_name in table_list:
            table_name = str(table_name).strip()
            if table_name not in seen:
                seen.add(table_name)
                clean_tables.append(table_name)
            
            # -- PATHWAY INTERACTIONS --
            if subsection_1 == "pathways":
                rows = app.interactions.get_interactions_table(compound, "pathway", where = where_pathways, order="pathwayid,asc")   
                app.utils.save_rows_json(rows, f"1.{subsection}/compound_{compound.cid}_{compound_name}_interactionstable.json")
                app.utils.save_rows_csv(rows, f"1.{subsection}/compound_{compound.cid}_{compound_name}_interactionstable.csv")
    

            # -- CHEMICAL-TARGET INTERACTIONS -- 
            elif subsection_1 == "chemical-target interactions": 
                for table_name in clean_tables:
                    if table_name.lower().startswith("collection="):
                        continue
                
                    rows = app.interactions.get_interactions_table(compound, table_name, order="geneid,asc")
                    app.utils.save_rows_json(rows, f"1.{subsection}/compound_{compound.cid}_{compound_name}_interactionstable.json")
                    app.utils.save_rows_csv(rows, f"1.{subsection}/compound_{compound.cid}_{compound_name}_interactionstable.csv")
    
            else:
                continue

            

def fetch_interactions_summary(compound_names):
    csv_files = [f"{app.utils.safe_filename(name)}_UniProt_proteins.csv" for name in compound_names]

    print(compound_names)

    dfs_int = []
    required_columns = [
        "compound",
        "geneid", 
        "symbol", 
        "description",
        "uniprot_accession",
        "uniprot_accessions",
    ]


    for csv_file in csv_files:
        print("Checking: "+csv_file)
        if not Path(csv_file).exists():
            continue
    
    # So even if it is empty, the application does not crash
        df = pd.read_csv(csv_file)
        for col in required_columns:
            if col not in df.columns:
                df[col] = ""

        dfs_int.append(pd.read_csv(csv_file))

    df_interactions = (
        pd.concat(dfs_int, ignore_index=True) 
        if dfs_int 
        else pd.DataFrame(columns = ["uniprot_accession", "compound", "symbol", "geneid"])
    )

    df_interactions.to_csv("ProteinInteractions.csv", index=False)
    return df_interactions


def fetch_pathway_summary(compound_names):
    dfs_proteins = []

    for name in compound_names:
        filename = f"{app.utils.safe_filename(name)}_ProteinsAllPathways.csv"
        compound_file = Path(filename)

        if compound_file.exists():
            dfs_proteins.append(pd.read_csv(compound_file))
            continue

    df_pathways = (
        pd.concat(dfs_proteins, ignore_index=True) 
        if dfs_proteins 
        else pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound"])
    )

    df_pathways.to_csv("AllProteinsPathways.csv", index=False)
    return df_pathways

def fill_missing_symbols(final_summary):
    final_summary = final_summary.copy()

    missing = (final_summary["symbol"].isna() | (final_summary["symbol"].astype(str).str.strip() == ""))
    missing_accessions = (final_summary.loc[missing, "uniprot_accession"].dropna().astype(str).str.strip())
    missing_accessions = [x for x in missing_accessions.unique() if x]

    if not missing_accessions:
        return final_summary
    df_symbols = app.proteins.map_uniprot_to_symbol(missing_accessions)

    final_summary = final_summary.merge(df_symbols, on="uniprot_accession", how="left")

    final_summary["symbol"] = final_summary.apply(
        lambda row: str(row["symbol"]).strip()
        if pd.notna(row["symbol"]) and str(row["symbol"]).strip()
        else (
            str(row["mapped_symbol"]).strip()
            if pd.notna(row["mapped_symbol"]) and str(row["mapped_symbol"]).strip()
            else ""
            ),
        axis=1
    )

    final_summary = final_summary.drop(columns=["mapped_symbol"], errors = "ignore")
    return final_summary


# I NEED TO CHECK IF THOSE FILES ALREADY EXIST OR NOT
def build_final_summary(df_interactions, df_pathways):
    df_interactions = df_interactions.copy()
    df_pathways = df_pathways.copy()

    df_interactions.columns = df_interactions.columns.str.strip()
    df_pathways.columns = df_pathways.columns.str.strip()

    df_interactions["uniprot_accession"] = (df_interactions["uniprot_accession"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NaN": ""})
    )
    df_interactions = df_interactions[df_interactions["uniprot_accession"] != ""].copy()

    df_pathways["uniprot_accession"] = (df_pathways["uniprot_accession"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NaN": ""})
    )
    df_pathways = df_pathways[df_pathways["uniprot_accession"] != ""].copy()

    interactions_summary = (
        df_interactions.groupby("uniprot_accession", as_index=False).agg(
            total_count = ("uniprot_accession", "size"),
            n_compounds = ("compound", "nunique"),
            compounds=("compound", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            symbol=("symbol", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            geneid=("geneid", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    pathway_summary = (
        df_pathways.groupby("uniprot_accession", as_index=False).agg(
            n_pathways = ("pathway", "nunique"),
            pathways = ("pathway", lambda x: ";".join(sorted(set(x)))),
            pathway_compounds = ("compound", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    final_summary = interactions_summary.merge(pathway_summary, on="uniprot_accession", how="outer")

    final_summary = final_summary.fillna({
        "total_count":0,
        "n_compounds":0,
        "compounds": "",
        "symbol": "",
        "geneid": "",
        "n_pathways": 0,
        "pathways": "",
        "pathway_compounds": "",
    })

    final_summary["total_count"] = final_summary["total_count"].astype(int)
    final_summary["n_pathways"] = final_summary["n_pathways"].astype(int)


    final_summary["source"] = final_summary.apply(
        lambda row:
        "interaction_and_pathway" if row["total_count"] >0 and row["n_pathways"] >0
        else "interaction_predominant" if row["total_count"] >0
        else "pathway_predominant" if row["n_pathways"] >0
        else "",
        axis = 1
    )

    final_summary = fill_missing_symbols(final_summary)


    def merge_compound_strings(*values):
        items = []
        seen = set()

        for value in values:
            if pd.isna(value) or str(value).strip() == "":
                continue

            for item in str(value).split(";"):
                item = item.strip()
                if item and item not in seen:
                    seen.add(item)
                    items.append(item)

        return ";".join(sorted(items))


    # Pathway-only proteins now inherit the compound(s) from pathways
    final_summary["compounds"] = final_summary.apply(
        lambda row: merge_compound_strings(
            row.get("compounds", ""), 
            row.get("pathway_compounds", "")
            ),
            axis=1
        )

    final_summary["n_compounds"] = final_summary["compounds"].apply(lambda x: len([v for v in str(x).split(";") if v.strip()]) if str(x).strip() else 0)
    
    final_summary = final_summary.sort_values(["total_count", "n_compounds"], ascending=[False,False])
    
    final_summary.to_csv("Protein_final_summary.csv", index=False)
    return final_summary


def build_go_enrichment(final_summary):
    df_go = app.proteins.fetch_goterms(final_summary, 
                                       aspects=["biological_process", "molecular_function", "cellular_component"],
                                       )
    
    df_go_empty =pd.DataFrame(columns=["uniprot_accession", "go_id", "go_name", "symbol", "aspect"])
    df_go_empty_aspect = pd.DataFrame(columns=["uniprot_accession", "go_id", "go_name", "symbol"])
    
    
    if df_go.empty:
        final_summaryGO = final_summary.copy()
        final_summaryGO["n_go_terms"] = 0
        final_summaryGO["go_ids"] = ""
        final_summaryGO["go_names"] = ""

        final_summaryGO["go_bp_ids"] = ""
        final_summaryGO["go_bp_names"] = ""
            
        final_summaryGO["go_mf_ids"] = ""
        final_summaryGO["go_mf_names"] = ""
            
        final_summaryGO["go_cc_ids"] = ""
        final_summaryGO["go_cc_names"] = ""
        
        return {
            "df_go": df_go_empty,
            "df_go_bp": df_go_empty_aspect.copy(),
            "df_go_mf": df_go_empty_aspect.copy(),
            "df_go_cc": df_go_empty_aspect.copy(),
            "final_summaryGO": final_summaryGO
        }

    # Add GO names to the GO ids
    go_name = app.proteins.fetch_gonames(df_go["go_id"].dropna().unique())
    df_go = df_go.merge(go_name, on="go_id", how="left")

    # Add symbols to the GO tables
    symbol_map = final_summary[["uniprot_accession", "symbol"]].drop_duplicates()
    df_go = df_go.merge(symbol_map, on="uniprot_accession", how="left")

    # Final DataFrame standardized and removing duplicates
    df_go = df_go[["uniprot_accession", "go_id","go_name", "symbol", "aspect"]].copy()
    df_go = df_go.drop_duplicates(subset=["uniprot_accession", "go_id", "aspect"]).copy()

    # Build 3 tables per aspect of GO aspects
    df_go_bp = (df_go[df_go["aspect"] == "biological_process"]
                [["uniprot_accession","go_id","go_name","symbol"]].drop_duplicates().reset_index(drop=True))
    df_go_mf = (df_go[df_go["aspect"] == "molecular_function"]
                [["uniprot_accession","go_id","go_name","symbol"]].drop_duplicates().reset_index(drop=True))
    df_go_cc = (df_go[df_go["aspect"] == "cellular_component"]
                [["uniprot_accession","go_id","go_name","symbol"]].drop_duplicates().reset_index(drop=True))

    # Build a summary table with all the aggregated GO information per protein
    go_summary = (
        df_go.groupby("uniprot_accession", as_index=False).agg(
        n_go_terms = ("go_id", lambda x: x.dropna().nunique()),
        go_ids = ("go_id", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        go_names = ("go_name", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    bp_summary = app.proteins.summarize_goaspect(df_go, "biological_process", "bp")
    mf_summary = app.proteins.summarize_goaspect(df_go, "molecular_function", "mf")
    cc_summary = app.proteins.summarize_goaspect(df_go, "cellular_component", "cc")


    final_summaryGO = final_summary.merge(go_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(bp_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(mf_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(cc_summary, on="uniprot_accession", how="left")
            
    final_summaryGO["n_go_terms"] = final_summaryGO["n_go_terms"].fillna(0).astype(int)

    for col in ["go_ids", "go_names", 
                "go_bp_ids", "go_bp_names", 
                "go_mf_ids", "go_mf_names", 
                "go_cc_ids", "go_cc_names"]:
        final_summaryGO[col] = final_summaryGO[col].fillna("")

    final_summaryGO.to_csv("Protein_final_summaryGO.csv", index=False)
    return {
        "df_go": df_go,
        "df_go_bp": df_go_bp,
        "df_go_mf": df_go_mf,
        "df_go_cc": df_go_cc,
        "final_summaryGO": final_summaryGO
    }


def run_full_pipeline(smiles_codes, email):
    compound_names = []
    print("Running Pipeline")

    for smiles in smiles_codes:
        try:
            compound, compound_info, df_proteins = fetch_pubchem_compound(smiles, email)
            compound_name = app.utils.safe_filename(app.chem.compound_display_name(compound))
            compound_names.append(compound_name)
            print(compound_name)
            print(compound_info)
        except Exception as e:
            print(f"Skipping {smiles}: {e}")
    
    df_interactions = fetch_interactions_summary(compound_names)
    df_pathways = fetch_pathway_summary(compound_names)

    final_summary = build_final_summary(df_interactions, df_pathways)

    go_results = build_go_enrichment(final_summary)

    go_results["final_summaryGO"].to_csv("Protein_final_summaryGO.csv", index=False)
    excelsummary = pd.read_csv("Protein_final_summaryGO.csv")
    excelsummary.to_excel("Protein_final_summaryGO.xlsx", index=False)

    return {
        "compound_names": compound_names,
        "df_interactions": df_interactions,
        "df_pathways": df_pathways,
        "final_summary": final_summary,
        "df_go": go_results["df_go"],
        "df_go_bp": go_results["df_go_bp"],
        "df_go_mf": go_results["df_go_mf"],
        "df_go_cc": go_results["df_go_cc"],
        "final_summaryGO": go_results["final_summaryGO"],
    }
