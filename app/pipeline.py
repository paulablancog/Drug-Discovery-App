import pandas as pd

import app.chem
import app.interactions
import app.proteins
import app.pathways
import app.utils

def fetch_pubchem_compound(smiles_code, email, selected_tax_ids=None):
    """Fetched the PubChem compound information, chemical-target interactions and pathways DataFrames for a given SMILES code. It includes a taxonomic filtering."""
    compound = app.chem.compound_retrieval(smiles_code)
    if compound is None:
        raise ValueError("No compound found. Check the SMILES code")
    
    compound_info = app.chem.compound_information(compound)
    compound_name = app.chem.compound_display_name(compound)

    interaction_data = load_interactions(compound, selected_tax_ids=selected_tax_ids) 
     
    df_geneids = app.proteins.retrieve_targets_1(compound_name, interaction_data["chemical_target_rows"], selected_tax_ids=selected_tax_ids)
    proteins_data = app.proteins.translate_geneid_to_protein(email, df_geneids, compound_name)
    tax_map = pd.DataFrame(columns=["geneid", "taxid", "taxname"])

    if df_geneids is not None and not df_geneids.empty:
        for col in ["geneid", "taxid", "taxname"]:
            if col not in df_geneids.columns:
                df_geneids[col] = ""
        tax_map = df_geneids[["geneid", "taxid", "taxname"]].drop_duplicates()
    
    df_map = app.proteins.map_genes_to_uniprot(df_geneids)
    df_uniprot_info = pd.DataFrame(columns=["uniprot_accession", "protein_name", "mapped_symbol", "taxid", "taxname"])

    if df_map is not None and not df_map.empty and "uniprot_accession" in df_map.columns:
        accessions = (
            df_map["uniprot_accession"].dropna().astype(str).str.strip()
        )
        accessions = [x for x in accessions.unique() if x]

        if accessions:
            df_uniprot_info = app.proteins.map_uniprot_to_info(accessions)

    protein_data = proteins_data.drop_duplicates(subset=["geneid"]) if "geneid" in proteins_data.columns else proteins_data

    if "geneid" not in df_map.columns:
        df_map = pd.DataFrame(columns=["geneid", "uniprot_accession"])
    
    if protein_data.empty:
        df_proteins = pd.DataFrame(columns=["compound", "cid", "geneid", "symbol", "description", "uniprot_accession", "protein_name", "taxid", "taxname"])
    else:
        df_proteins = protein_data.merge(df_map, on="geneid", how="left")
        df_proteins = df_proteins.merge(tax_map, on="geneid", how="left")

        # separate taxid and taxname from UniProt and PubChem (UniProt taxonomic information here)
        df_uniprot_info = df_uniprot_info.rename(columns={
            "taxid": "taxid_uniprot",
            "taxname": "taxname_uniprot"
        })

        df_proteins = df_proteins.merge(
            df_uniprot_info,
            on="uniprot_accession",
            how="left",
        )

        # PubChem taxonomy information separated from UniProt's
        df_proteins["pubchem_taxid"] = df_proteins["taxid"].fillna("").astype(str).str.strip()
        df_proteins["pubchem_taxname"] = df_proteins["taxname"].fillna("").astype(str).str.strip()

        # Prioritize the taxonomic information of UniProt instead of PubChem's
        df_proteins["taxid"] = df_proteins.apply(
            lambda row: str(row["taxid_uniprot"]).strip()
            if pd.notna(row.get("taxid_uniprot")) and str(row.get("taxid_uniprot")).strip()
            else str(row.get("pubchem_taxid", "")).strip(),
            axis = 1
        )

        df_proteins["taxname"] = df_proteins.apply(
            lambda row: str(row["taxname_uniprot"]).strip()
            if pd.notna(row.get("taxname_uniprot")) and str(row.get("taxname_uniprot")).strip()
            else str(row.get("pubchem_taxname", "")).strip(),
            axis = 1
        )

        selected_tax_ids = app.interactions.normalize_taxonomy_ids(selected_tax_ids)
        if selected_tax_ids:
            df_proteins = df_proteins[df_proteins["taxid"].astype(str).str.strip().isin(selected_tax_ids)].copy()

        df_proteins["symbol"] = df_proteins.apply(
            lambda row: str(row["symbol"]).strip()
            if pd.notna(row.get("symbol")) and str(row.get("symbol")).strip()
            else str(row.get("mapped_symbol", "")).strip(),
            axis = 1
        )

        df_proteins["compound"] = compound_name
        df_proteins["cid"] = compound_info.get("cid")
        df_proteins = df_proteins.drop(columns=["taxid_uniprot", "taxname_uniprot", "pubchem_taxid", "pubchem_taxname", "mapped_symbol"],
        errors = "ignore")

    df_pathways = app.pathways.retrieve_pathways(compound, interaction_data["pathway_rows"] ,compound_name, selected_tax_ids=selected_tax_ids)
    if df_pathways is None:
        df_pathways = pd.DataFrame(columns=["uniprot_accession", "protein_name", "symbol", "pathway", "pathway_name", "compound", "cid", "taxid", "taxname"])

    return compound, compound_info, compound_name, df_proteins, df_pathways


def load_interactions(compound, selected_tax_ids=None):
    """Load PubChem chemical-target interactions and pathways for a given compound, with optional taxonomy filtering."""
    # 1. Creates the URL
    index_url = f"{app.utils.URL_BASE}/rest/pug_view/index/compound/{compound.cid}/JSON"
    # 2. Retrieves the index JSON
    index_json = app.utils.get_json(index_url)

    if index_json is None:
        raise ValueError("Failed to retrieve index JSON")

    if not app.interactions.has_interactions_and_pathways(index_json):
        return{
            "chemical_target_rows": [],
            "pathway_rows": []
        }
    
    # 4. Gets all sections in the index JSON
    # 5. Finds the Interactions and Pathways section and retrieves the data of interactions
    data = app.interactions.load_interactions_and_pathways_data(compound)
    
    if data is None:
        return{
            "chemical_target_rows": [],
            "pathway_rows": []
        }

    # 6. Get external tables
    tables = app.interactions.retrieve_externaltable(data)

    chemical_target_rows = []
    pathway_rows = []

    # 7. Download and save each table
    where_pathways = {"ands": [{"cid": str(compound.cid)}, {"core": "1"}]}

    for subsection, table_list in tables:
        subsection_1 = (subsection or "").strip().lower()

        clean_tables = []
        seen = set()

        for table_name in table_list:
            table_name = str(table_name).strip()
            if table_name not in seen:
                seen.add(table_name)
                clean_tables.append(table_name)
            
        # -- PATHWAY INTERACTIONS --
        if subsection_1 == "pathways":
            rows = app.interactions.get_interactions_table(compound, "pathway", where = where_pathways, order="pathwayid,asc", selected_tax_ids=selected_tax_ids)   
            pathway_rows.extend(rows)
    

        # -- CHEMICAL-TARGET INTERACTIONS -- 
        elif subsection_1 == "chemical-target interactions": 
            for table_name in clean_tables:
                if table_name.lower().startswith("collection="):
                    continue
                
                rows = app.interactions.get_interactions_table(compound, table_name, order="geneid,asc", selected_tax_ids=selected_tax_ids)
                chemical_target_rows.extend(rows)
    
    return {
        "chemical_target_rows": chemical_target_rows,
        "pathway_rows": pathway_rows,
    }

            

def fetch_interactions_summary(proteins):
    """Given a list of DataFrames with chemical-target interactions, 
    it returns a DataFrame with the relevant information of the interaction proteins involved"""
    required_columns = [
        "compound",
        "cid",
        "geneid", 
        "symbol", 
        "description",
        "uniprot_accession",
        "protein_name",
        "taxid",
        "taxname",
    ]

    dfs_int = []
    for df in proteins:
        if df is None:
            continue
        df = df.copy()
        for col in required_columns:
            if col not in df.columns:
                df[col] = ""
        dfs_int.append(df[required_columns])

    return (pd.concat(dfs_int, ignore_index=True)
            if dfs_int
            else pd.DataFrame(columns=required_columns)
    )


def fetch_pathway_summary(pathways):
    """Given a list of DataFrames with pathway interactions, 
    it returns a DataFrame with the relevant information of the pathway proteins involved"""
    required_columns = [
        "uniprot_accession",
        "protein_name",
        "symbol",
        "pathway",
        "pathway_name", 
        "compound", 
        "cid",
        "taxid",
        "taxname"
    ]

    dfs_proteins = []

    for df in pathways:
        if df is None:
            continue
        df = df.copy()
        for col in required_columns:
            if col not in df.columns:
                df[col] = ""
        dfs_proteins.append(df[required_columns])

    df_pathways = (
        pd.concat(dfs_proteins, ignore_index=True) 
        if dfs_proteins 
        else pd.DataFrame(columns=required_columns)
    )

    df_groupedpathways = app.pathways.group_pathways(df_pathways)

    return df_pathways, df_groupedpathways

def fill_missing_symbols(final_summary):
    """Fills missing gene symbols in the final summary by mapping the UniProt accessions to gene symbols using the UniProt API"""
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


def build_final_summary(df_interactions, df_pathways):
    """Given the chemical-target interactions and pathway interactions DataFrames, it builds a final summary 
    DataFrame by merging information and aggregating it by protein (UniProt accession)"""
    df_interactions = df_interactions.copy()
    df_pathways = df_pathways.copy()

    df_interactions.columns = df_interactions.columns.str.strip()
    df_pathways.columns = df_pathways.columns.str.strip()

    for col in ["uniprot_accession", "protein_name", "compound", "geneid", "symbol", "description","taxid", "taxname"]:
        if col in df_interactions.columns:
            df_interactions[col] = (df_interactions[col]
                .astype(str)
                .str.strip()
                .replace({"nan": "", "None": "", "NaN": ""})
            )
        else:
            df_interactions[col] = ""

   
    df_interactions = df_interactions[df_interactions["uniprot_accession"] != ""].copy()
    
    for col in ["uniprot_accession", "protein_name", "pathway", "pathway_name", "compound", "taxid", "taxname"]:
        if col in df_pathways.columns:
            df_pathways[col] = (df_pathways[col]
                .astype(str)
                .str.strip()
                .replace({"nan": "", "None": "", "NaN": ""})
            )
        else:
            df_pathways[col] = ""

   
    df_pathways = df_pathways[df_pathways["uniprot_accession"] != ""].copy()

    interactions_summary = (
        df_interactions.groupby("uniprot_accession", as_index=False).agg(
            interaction_count = ("uniprot_accession", "size"),
            n_compounds = ("compound", "nunique"),
            compounds=("compound", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            symbol=("symbol", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            geneid=("geneid", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            interaction_protein_name=("protein_name", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            interaction_taxid = ("taxid", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            interaction_taxname = ("taxname", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    pathway_summary = (
        df_pathways.groupby("uniprot_accession", as_index=False).agg(
            pathway_count = ("uniprot_accession", "size"),
            n_pathways = ("pathway", "nunique"),
            pathways = ("pathway", lambda x: ";".join(sorted(set(x)))),
            pathway_names = ("pathway_name", lambda x: ";".join(sorted(set(x)))),
            pathway_compounds = ("compound", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            pathway_protein_name = ("protein_name", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            pathway_taxid = ("taxid", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            pathway_taxname = ("taxname", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    final_summary = interactions_summary.merge(pathway_summary, on="uniprot_accession", how="outer")

    final_summary = final_summary.fillna({
        "interaction_count":0,
        "pathway_count": 0,
        "n_compounds":0,
        "compounds": "",
        "symbol": "",
        "geneid": "",
        "n_pathways": 0,
        "pathways": "",
        "pathway_compounds": "",
        "pathway_names": "",
        "interaction_taxid": "",
        "interaction_taxname": "",
        "pathway_taxid": "",
        "pathway_taxname": "",
        "pathway_protein_name": "",
        "interaction_protein_name": "",
    })

    final_summary["interaction_count"] = final_summary["interaction_count"].astype(int)
    final_summary["pathway_count"] = final_summary["pathway_count"].astype(int)
    final_summary["n_pathways"] = final_summary["n_pathways"].astype(int)

    final_summary["total_count"] = (final_summary["interaction_count"] + final_summary["pathway_count"])

    final_summary["source"] = final_summary.apply(
        lambda row: "interaction_and_pathway" 
        if row["interaction_count"] >0 and row["pathway_count"] >0
        else "interaction_predominant" if row["interaction_count"] >0 and row["pathway_count"] == 0
        else "pathway_predominant" if row["pathway_count"] >0 and row["interaction_count"] == 0
        else "",
        axis = 1
    )

    final_summary = fill_missing_symbols(final_summary)

    def merge_compound_strings(*values):
        """Given multiple strings of compounds separated by semicolons, it merges them into a single string with unique compounds."""
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
    
    final_summary["taxid"] = final_summary.apply(
        lambda row: merge_compound_strings(
            row.get("interaction_taxid", ""),
            row.get("pathway_taxid", "")
        ),
        axis=1
    )

    final_summary["taxname"] = final_summary.apply(
        lambda row: merge_compound_strings( 
            row.get("interaction_taxname", ""),
            row.get("pathway_taxname", "")
        ),
        axis=1
    )
    
    final_summary["protein_name"] = final_summary.apply(
        lambda row: merge_compound_strings(
            row.get("pathway_protein_name", ""),
            row.get("interaction_protein_name", "")
        ),
        axis=1
    )

    final_summary["n_compounds"] = final_summary["compounds"].apply(
        lambda x: len([v for v in str(x).split(";") if v.strip()]) if str(x).strip() else 0
    )
    
    final_summary["source_order"] = final_summary["source"].map({
        "interaction_predominant":0,
        "interaction_and_pathway":1,
        "pathway_predominant":2,
        "": 3,
    })
    final_summary = final_summary.sort_values(
        ["source_order", "n_compounds", "total_count"],
        ascending=[True, False, False]
    ).drop(columns=["source_order"]).reset_index(drop=True)

    return final_summary


def build_go_enrichment(final_summary):
    """Given the final summary DataFrame, it builds the GO enrichment results by 
    fetching the GO terms of the proteins and grouping them by aspect."""
    df_go = app.proteins.fetch_goterms(final_summary, 
                                       aspects=["biological_process", "molecular_function", "cellular_component"],
                                       )
    
    df_go_empty =pd.DataFrame(columns=["uniprot_accession", "go_id", "go_name", "symbol", "aspect", "compounds"])
    df_go_empty_aspect = pd.DataFrame(columns=["uniprot_accession", "go_id", "go_name", "symbol", "compounds"])
    empty_grouped = pd.DataFrame(columns=["go_name", "go_id", "n_proteins", "n_compounds", "proteins", "compounds"])
    
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
            "df_go_bp_grouped": empty_grouped.copy(),
            "df_go_mf_grouped": empty_grouped.copy(),
            "df_go_cc_grouped": empty_grouped.copy(),
            "final_summaryGO": final_summaryGO
        }

    # Add GO names to the GO ids
    go_name = app.proteins.fetch_gonames(df_go["go_id"].dropna().unique())
    df_go = df_go.merge(go_name, on="go_id", how="left")

    # Add symbols to the GO tables
    symbol_map = final_summary[["uniprot_accession", "symbol"]].drop_duplicates()
    df_go = df_go.merge(symbol_map, on="uniprot_accession", how="left")

    # Add compounds to the GO tables
    compounds_map = final_summary[["uniprot_accession", "compounds"]].drop_duplicates()
    df_go = df_go.merge(compounds_map, on="uniprot_accession", how="left")

    # Final DataFrame standardized and removing duplicates
    df_go = df_go[["uniprot_accession", "go_id","go_name", "symbol", "aspect", "compounds"]].copy()
    df_go = df_go.drop_duplicates(subset=["uniprot_accession", "go_id", "aspect"]).copy()

    # Build 3 tables per aspect of GO aspects
    df_go_bp = (df_go[df_go["aspect"] == "biological_process"]
                [["uniprot_accession","go_id","go_name","symbol", "compounds"]].drop_duplicates().reset_index(drop=True))
    df_go_mf = (df_go[df_go["aspect"] == "molecular_function"]
                [["uniprot_accession","go_id","go_name","symbol", "compounds"]].drop_duplicates().reset_index(drop=True))
    df_go_cc = (df_go[df_go["aspect"] == "cellular_component"]
                [["uniprot_accession","go_id","go_name","symbol", "compounds"]].drop_duplicates().reset_index(drop=True))

    df_go_bp_grouped = app.proteins.group_goterms(df_go_bp)
    df_go_mf_grouped = app.proteins.group_goterms(df_go_mf)
    df_go_cc_grouped = app.proteins.group_goterms(df_go_cc)

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

    return {
        "df_go": df_go,
        "df_go_bp": df_go_bp,
        "df_go_mf": df_go_mf,
        "df_go_cc": df_go_cc,
        "df_go_bp_grouped": df_go_bp_grouped,
        "df_go_mf_grouped": df_go_mf_grouped,
        "df_go_cc_grouped": df_go_cc_grouped,
        "final_summaryGO": final_summaryGO
    }


def run_full_pipeline(smiles_codes, email, selected_tax_ids=None, ui = None):
    """Runs the full pipeline of fetching compound information, chemical-target and pathways information, 
    building summaries and GO enrichment for a list of SMILES codes with UI updates and taxonomic filtering."""
    compound_names = []
    all_compounds = []
    proteins = []
    pathways = []

    if ui:
        ui["status_box"].info("Running analysis... This may take a few minutes.")
        ui["progress_bar"].progress(5, text="Identifying compounds...")

    for i, smiles in enumerate(smiles_codes, start=1):
        try:
            compound, compound_info, compound_name, df_proteins, df_pathways = fetch_pubchem_compound(smiles, email, selected_tax_ids=selected_tax_ids)
            
            compound_names.append(compound_name)
            proteins.append(df_proteins)
            pathways.append(df_pathways)

            all_compounds.append({
                "smiles":smiles,
                "compound_name": compound_name,
                "cid": compound_info.get("cid"),
                "molecular_formula": compound_info.get("molecular_formula"),
                "molecular_weight": compound_info.get("molecular_weight"),
                "status": "Identified",
            })

            if ui:
                partial_df = pd.DataFrame(all_compounds)
                ui["compound_box"].markdown("### Compounds identified")
                ui["compound_box"].dataframe(partial_df, width="stretch")
                pct = 5 + int(25*i/max(len(smiles_codes),1))
                ui["progress_bar"].progress(pct, text=f"Compound identified {i} of {len(smiles_codes)}...")

        except Exception as e:
            all_compounds.append({
                "smiles":smiles,
                "compound_name": None,
                "cid": None,  
                "molecular_formula": None,
                "molecular_weight": None,
                "status": f"Error: {str(e)}",
            })
            print(f"Skipping {smiles}: {e}")
    
    compound_results = pd.DataFrame(all_compounds)
    
    df_interactions = fetch_interactions_summary(proteins)
    if ui:
        ui["status_box"].info("Compound-Protein interactions retrieved")
        ui["interactions_box"].markdown("### Compound-Protein interactions")
        ui["interactions_box"].dataframe(df_interactions, width="stretch")
        ui["progress_bar"].progress(55, text="Compound-Protein interactions completed.")

    df_pathways, df_groupedpathways = fetch_pathway_summary(pathways)
    if ui:
        ui["status_box"].info("Pathways retrieved")
        ui["pathway_box"].markdown("### Pathways")
        ui["pathway_box"].dataframe(df_groupedpathways, width="stretch")
        ui["progress_bar"].progress(70, text="Pathways completed.")

    final_summary = build_final_summary(df_interactions, df_pathways)
    if ui:
        ui["summary_box"].markdown("### Protein Summary")
        ui["summary_box"].dataframe(final_summary[[
            "uniprot_accession",
            "protein_name",
            "symbol",
            "taxid",
            "taxname",
            "interaction_count",
            "pathway_count",
            "total_count",
            "compounds",
            "n_compounds",
            "n_pathways",
            "pathways",
            "pathway_names",
            "pathway_compounds",
            "source",
        ]], width="stretch")
        ui["progress_bar"].progress(85, text="Protein summary completed.")

    go_results = build_go_enrichment(final_summary)

    if ui:
        ui["status_box"].success("Analysis completed!")
        ui["progress_bar"].progress(100, text="Done")

    return {
        "compound_names": compound_names,
        "compound_results": compound_results,
        "df_interactions": df_interactions,
        "df_pathways": df_pathways,
        "df_groupedpathways": df_groupedpathways,
        "final_summary": final_summary,
        "df_go": go_results["df_go"],
        "df_go_bp": go_results["df_go_bp"],
        "df_go_mf": go_results["df_go_mf"],
        "df_go_cc": go_results["df_go_cc"],
        "df_go_bp_grouped": go_results["df_go_bp_grouped"],
        "df_go_mf_grouped": go_results["df_go_mf_grouped"],
        "df_go_cc_grouped": go_results["df_go_cc_grouped"],
        "final_summaryGO": go_results["final_summaryGO"],
    }
