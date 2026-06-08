import io
import pandas as pd
import streamlit as st

page_selection = {
    "input_smiles": "Input SMILES",
    "compounds": "Compounds",
    "interactions": "Chemical-Target Interactions",
    "pathways_general": "Pathways - general table",
    "protein_pathways": "Pathways - proteins",
    "go_enrichment": "GO Enrichment - general table",
    "go_proteins": "GO Enrichment - proteins",
    "protein_summary": "Protein Summary",
}


def download_excel_analysis(selected_page = None):
    """Generates an Excel file with the analysis results, each section corresponds to a new page in the Excel file."""
    results = st.session_state.get("results") or {}
    submitted_smiles = st.session_state.get("submitted_smiles", [])

    if selected_page is None:
        selected_page = list(page_selection.keys())

    selected_page = set(selected_page)
    
    df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
    df_pathways = results.get("df_pathways", pd.DataFrame()).copy()
    df_groupedpathways = results.get("df_groupedpathways", pd.DataFrame()).copy()
    
    df_go_bp_grouped = results.get("df_go_bp_grouped", pd.DataFrame()).copy()
    df_go_mf_grouped = results.get("df_go_mf_grouped", pd.DataFrame()).copy()
    df_go_cc_grouped = results.get("df_go_cc_grouped", pd.DataFrame()).copy()
    
    df_go_bp = results.get("df_go_bp", pd.DataFrame()).copy()
    df_go_mf = results.get("df_go_mf", pd.DataFrame()).copy()
    df_go_cc = results.get("df_go_cc", pd.DataFrame()).copy()

    output = io.BytesIO()

    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # 1. Input SMILES
        if "input_smiles" in selected_page: 
            df_smiles = pd.DataFrame({"input_smiles": submitted_smiles})
            df_smiles.to_excel(writer, sheet_name = "1.Input_SMILES", header = True, index=False)

        # 2. Compounds
        if "compounds" in selected_page:
            compound_results = st.session_state.get("compound_results", pd.DataFrame()).copy()
            
            if compound_results.empty:
                compound_results = results.get("compound_results", pd.DataFrame()).copy()

            compound_table = [
                "smiles",
                "compound_name",
                "cid",
                "molecular_formula",
                "molecular_weight",
                "status",
            ]

            for col in compound_table:
                if col not in compound_results.columns:
                    compound_results[col] = ""

            compound_results = compound_results[compound_table]

            if not compound_results.empty:
                compound_results.to_excel(writer, sheet_name = "2.Compounds", header= True, index=False)
            else:
                pd.DataFrame(columns=compound_table).to_excel(writer, sheet_name = "2.Compounds", index=False)


            if not results:
                pd.DataFrame().to_excel(writer, sheet_name = "3.Interactions", index=False)
                pd.DataFrame().to_excel(writer, sheet_name = "4.Pathways", index=False)
                pd.DataFrame().to_excel(writer, sheet_name = "5.GO_Enrichment", index=False)
                pd.DataFrame().to_excel(writer, sheet_name = "6.Protein_Summary", index=False)
                output.seek(0)
                return output.getvalue()
            
        # 3. Interactions
        if "interactions" in selected_page:
            df_interactions.to_excel(writer, sheet_name = "3.Interactions", header = True, index=False)

        # 4.1 Pathways general table
        if "pathways_general" in selected_page:
            pathway_general_cols = [
                "pathway",
                "pathway_name",
                "n_compounds",
                "compounds",
                "n_proteins",
                "proteins",
                "uniprot_accessions",
                "taxid",
                "taxname",
            ]

            for col in pathway_general_cols:
                if col not in df_groupedpathways.columns:
                    df_groupedpathways[col] = ""

            df_groupedpathways = df_groupedpathways[pathway_general_cols]

            df_groupedpathways.to_excel(writer, sheet_name = "4.Pathways_General", header = True, index=False)

        # 4.2 Pathways proteins     
        if "protein_pathways" in selected_page:
            pathway_proteins_cols = [
                "pathway",
                "pathway_name",
                "uniprot_accession",
                "protein_name",
                "compound",
                "taxid",
                "taxname",
            ]

            for col in pathway_proteins_cols:
                if col not in df_pathways.columns:
                    df_pathways[col] = ""

            df_pathways = df_pathways[pathway_proteins_cols]

            df_pathways.to_excel(writer, sheet_name = "4.Pathways_Proteins", header = True, index=False)    

        # 5.1 GO Enrichment general table
        if "go_enrichment" in selected_page:
            go_general_cols = [
                "go_name",
                "go_id",
                "aspect",
                "n_compounds",
                "compounds",
                "n_proteins",
            ]
            go_general_aspects = []

            go_grouped_aspects = [
                ("Biological Process", df_go_bp_grouped),
                ("Molecular Function", df_go_mf_grouped),
                ("Cellular Component", df_go_cc_grouped),
            ]

            for aspect, grouped_df in go_grouped_aspects:
                grouped_df = grouped_df.copy()
                if grouped_df.empty:
                    continue
                grouped_df.insert(0, "aspect", aspect)

                for col in go_general_cols:
                    if col not in grouped_df.columns:
                        grouped_df[col] = ""

                go_general_aspects.append(grouped_df[go_general_cols])

            if go_general_aspects:
                df_go_general = pd.concat(go_general_aspects, ignore_index=True)
            else:
                df_go_general = pd.DataFrame(columns=go_general_cols)
            
            df_go_general.to_excel(writer, sheet_name = "5.GO_Enrichment_General", header = True, index=False)

        # 5.2 GO Enrichment proteins
        if "go_proteins" in selected_page:
            go_protein_cols = [
                "aspect",
                "go_name",
                "go_id",
                "uniprot_accession",
                "protein_name",
                "symbol",
                "taxid",
                "taxname",
            ]

            go_protein_parts = []

            go_detail_aspects = [
                ("Biological Process", df_go_bp),
                ("Molecular Function", df_go_mf),
                ("Cellular Component", df_go_cc),
            ]

            for aspect, detail_df in go_detail_aspects:
                detail_df = detail_df.copy()
                
                if detail_df.empty:
                    continue
                
                detail_df.insert(0, "aspect", aspect)
                go_protein_parts.append(detail_df)

            if go_protein_parts:
                df_go_proteins = pd.concat(go_protein_parts, ignore_index=True)
            else:
                df_go_proteins = pd.DataFrame(columns=go_protein_cols)

            # Add protein_name, taxid and taxname from summary final table
            df_summary = results.get("final_summaryGO", pd.DataFrame()).copy()
            if df_summary.empty:
                df_summary = results.get("final_summary", pd.DataFrame()).copy()
            
            summary_cols = ["uniprot_accession", "protein_name", "symbol", "taxid", "taxname"]
            
            if not df_summary.empty:
                for col in summary_cols:
                    if col not in df_summary.columns:
                        df_summary[col] = ""

                df_summary = df_summary[summary_cols].drop_duplicates(subset=["uniprot_accession"], keep="first")
                
                df_go_proteins = df_go_proteins.merge(df_summary, on=["uniprot_accession"], how="left", suffixes=("", "_summary"))

            for col in go_protein_cols:
                if col not in df_go_proteins.columns:
                    df_go_proteins[col] = ""
            df_go_proteins = df_go_proteins[go_protein_cols]

            df_go_proteins.to_excel(writer, sheet_name = "5.GO_Enrichment_Proteins", header = True, index=False)  
            
        # 6. Protein Summary
        if "protein_summary" in selected_page:
            df_proteinsummary = results.get("final_summaryGO", pd.DataFrame()).copy()
            df_proteinsummary.to_excel(writer, sheet_name = "6.Protein_Summary", header = True, index=False)

    output.seek(0)
    return output.getvalue()

