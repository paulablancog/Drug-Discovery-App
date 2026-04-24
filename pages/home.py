import io
import streamlit as st
import pandas as pd
from app.pipeline import run_full_pipeline
import app.chem

# This APP when entering new SMILES code different from the older ones, it DELETES the analysis results
# TODO -> maybe something to look at for the Add SMILES and the big text area (the big text area could do that
# and the Add SMILES could just add the new SMILES code without deleting the analysis and re running it again with
# that additional SMILES PREGUNTAR)

# For further load into the app
example_smiles = [
    "CC(C[N+](C)(C)C)OC(=O)N",
]

# Triple quotes allow multiple lines
# st.markdown() : for text you want formatted with Markdown
# st.write() is more general-purpose and can display text, dataFrames, lists and other Python objects
st.title("Drug Discovery Analysis")
st.markdown(
    """
    Enter one or more **SMILES codes** to prepare the analysis. 
    Use **one SMILES code per line**
    """)


def store_smiles(text):
    """Collects the SMILES code per line introduced by the user and stores it into a list of unique SMILES strings"""
    smiles = []
    for line in text.splitlines():
        value = line.strip()
        if value:
            smiles.append(value)
    unique_smiles = []
    seen = set()

    for item in smiles:
        if item not in seen:
            seen.add(item)
            unique_smiles.append(item)

    return unique_smiles


# For clearing the session:
def clear_session():
    for key in [
        "submitted_smiles", 
        "email_input",
        "submitted_email", 
        "results",
        "analysis_ready",
        "run_error",
        "new_smiles_input",
        "show_add_smiles",
        "smiles_to_delete",
        "show_delete_smiles",
        "compound_results",
        "compound_results_source"
        ]:
        if key in st.session_state:
            del st.session_state[key]

if "submitted_smiles" not in st.session_state:
    st.session_state["submitted_smiles"] = []

if "submitted_email" not in st.session_state:
    st.session_state["submitted_email"] = ""

if "new_smiles_input" not in st.session_state:
    st.session_state["new_smiles_input"] = ""

if "show_add_smiles" not in st.session_state:
    st.session_state["show_add_smiles"] = False

if "show_delete_smiles" not in st.session_state:
    st.session_state["show_delete_smiles"] = False

if "smiles_to_delete" not in st.session_state:
    st.session_state["smiles_to_delete"] = []

# This method is for clearing the results of the analysis only when the user modifies the SMILES input (should take into account the deletion of a SMILES code)
def clear_analysis_results():
    for key in ["results", "analysis_ready", "run_error"]:
        if key in st.session_state:
            del st.session_state[key]

# These methods modify the st.session_state of the widgets that were already initialized with a key
# It DELETES the analysis results previously made
def load_example():
    if "submitted_smiles" not in st.session_state:
        st.session_state["submitted_smiles"] = []

    st.session_state["smiles_input"] = "\n".join(example_smiles)
    clear_analysis_results()

def clear_inputs():
    st.session_state["smiles_input"] = ""
    clear_session()

def show_additional_input():
    st.session_state["show_add_smiles"] = True

def store_additional_smiles():
    value = st.session_state.get("new_smiles_input", "").strip()
    if not value:
        return
    
    # Here I initialize "submitted_smiles", but it is optimal if the 
    # user "initializes" it by firstly entering SMILES code into the 
    # big text area 
    if "submitted_smiles" not in st.session_state:
        st.session_state["submitted_smiles"] = []

    if value not in st.session_state["submitted_smiles"]:
        st.session_state["submitted_smiles"].append(value)

    st.session_state["new_smiles_input"] = ""
    st.session_state["show_add_smiles"] = False
    clear_analysis_results()

def show_delete_smiles():
    st.session_state["show_delete_smiles"] = True

def update_smiles():
    selected = st.session_state.get("smiles_to_delete", [])
    current_smiles = st.session_state.get("submitted_smiles", [])

    st.session_state["submitted_smiles"] = [
        smiles for smiles in current_smiles if smiles not in selected
    ]
    st.session_state["smiles_to_delete"] = []
    st.session_state["show_delete_smiles"] = False
    clear_analysis_results()


def save_input():
    cleaned_smiles = store_smiles(st.session_state.get("smiles_input", ""))
    
    if not cleaned_smiles and not st.session_state.get("submitted_smiles", []):
        st.session_state["run_error"] = "Please enter at least one SMILES code"
        return

    for smiles in cleaned_smiles:
        if smiles not in st.session_state.get("submitted_smiles", []):
            st.session_state["submitted_smiles"].append(smiles)
    

    email_input = st.session_state.get("email_input", "").strip()
    saved_email = st.session_state.get("submitted_email", "").strip()

    email = email_input or saved_email

    if not email:
        st.session_state["run_error"] = "Please enter a valid email"
        return
    
    st.session_state["submitted_email"] = email
    st.session_state["email_input"] = email
    st.session_state["smiles_input"] = ""
    st.session_state["run_error"] = ""
    clear_analysis_results()



def run_analysis():
    smiles_codes = st.session_state.get("submitted_smiles", [])
    email = st.session_state.get("submitted_email", "")

    # Tengo que meter errores TODO
    if not email:
        st.session_state["run_error"] = "Please enter a valid email before running the analysis."
        st.session_state["analysis_ready"] = False
        return
    if not smiles_codes:
        st.session_state["run_error"] = "Please enter at least one SMILES code before running the analysis."
        st.session_state["analysis_ready"] = False
        return
    try:
        with st.spinner("Running analysis... This may take a few minutes..."):
            results = run_full_pipeline(smiles_codes, email, ui = {
                "status_box": status_box,
                "progress_bar": progress_bar,
                "compound_box": compound_box,
                "interactions_box": interactions_box,
                "pathway_box": pathway_box,
                "summary_box": summary_box

            },
        )

        st.session_state["results"] = results
        # So the DataFrame stays synchronized after analysis
        if "compound_results" in results:
            st.session_state["compound_results"] = results["compound_results"].copy()
            st.session_state["compound_results_source"] = tuple(smiles_codes)
        
        st.session_state["analysis_ready"] = True
        st.session_state["run_error"] = ""
    except Exception as e:
        st.session_state["run_error"] = f"An error occurred during the analysis: {str(e)}"
        st.session_state["analysis_ready"] = False


email = st.text_input( # this is the email value that is going to be saved
        "Email",
        key = "email_input",
        placeholder="name@example.com",
        help = "Required for analysis",
    )

    # Creates a multi-line text box
smiles_input = st.text_area(
        "SMILES codes",
        key = "smiles_input", # gives this widget an internal name in st.session_state
        height = 250,
        placeholder=(
            "Paste one SMILES code per line\n"
            "Example:\n"
            "CC(C[N+](C)(C)C)OC(=O)N)"
        ),
    )


def compound_hyperlink(cid, compound_name):
    if pd.notna(cid) and str(cid).strip() not in {"", "None", "nan"}:
        url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{int(cid)}"
        return f'<a href="{url}" target="_blank">{compound_name}</a>'
    return compound_name
          
def uniprot_hyperlink(accession, ):
    accession = str(accession).strip()
    if accession and accession.lower() not in {"none", "nan"}:
        url = f"https://www.uniprot.org/uniprotkb/{accession}/entry"
        return f'<a href="{url}" target="_blank">{accession}</a>'
    return accession

def pathway_hyperlink(pwacc):
    pwacc = str(pwacc).strip()
    if pwacc and pwacc.lower() not in {"none", "nan"}:
        url = f"https://pubchem.ncbi.nlm.nih.gov/pathway/{pwacc}"
        return f'<a href="{url}" target="_blank">{pwacc}</a>'
    return pwacc 

def goterm_hyperlink(go_id, go_name=None):
    go_id = str(go_id).strip()
    go_name = str(go_name).strip() if go_name else go_id
    if go_id and go_id.lower() not in {"none", "nan"}:
        url = f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"
        return f'<a href="{url}" target="_blank">{go_name}</a>'
    return go_name  

def geneid_hyperlink(geneid):
    geneid = str(geneid).strip()
    if geneid and geneid.lower() not in {"none", "nan", ""}:
        url = f"https://www.ncbi.nlm.nih.gov/gene/{geneid}"
        return f'<a href="{url}" target="_blank">{geneid}</a>'
    return geneid
    
    

def compounds_text_with_links(compounds_text, compound_results):
    if pd.isna(compounds_text):
        return ""
    
    compound_to_cid = dict(
        zip(
            compound_results["compound_name"].astype(str).str.strip(),
            compound_results["cid"]
        )
    )
    links = []
    for compound in str(compounds_text).split(";"):
        compound = compound.strip()
        if not compound:
            continue
        cid = compound_to_cid.get(compound)
        if pd.notna(cid) and str(cid).strip() not in {"", "None", "nan"}:
            url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{int(cid)}"
            links.append(f'<a href="{url}" target="_blank">{compound}</a>')
        else:
            links.append(compound)

    return "; ".join(links)

def uniprot_accessions_text_with_links(accessions_text):
    if pd.isna(accessions_text):
        return ""
    
    links = []
    for acc in str(accessions_text).split(";"):
        acc = acc.strip()
        if not acc:
            continue
        links.append(uniprot_hyperlink(acc))

    return "; ".join(links)

def pathway_text_with_links(pathway_text):
    if pd.isna(pathway_text):
        return ""
    
    links = []
    for pwacc in str(pathway_text).split(";"):
        pwacc = pwacc.strip()
        if not pwacc:
            continue
        links.append(pathway_hyperlink(pwacc))

    return "; ".join(links)


def download_excel_analysis():
    results = st.session_state.get("results", None)
    submitted_smiles = st.session_state.get("submitted_smiles", [])
    
    output = io.BytesIO()

    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # 1. Input SMILES
        df_smiles = pd.DataFrame({"input_smiles": submitted_smiles})
        df_smiles.to_excel(writer, sheet_name = "1.Input_SMILES", header = True, index=False)

        # 2. Compounds
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

        # If analysis has not been processed TODO -> won't be able to download

        if not results:
            pd.DataFrame().to_excel(writer, sheet_name = "3.Interactions", index=False)
            pd.DataFrame().to_excel(writer, sheet_name = "4.Pathways", index=False)
            pd.DataFrame().to_excel(writer, sheet_name = "5.GO_Enrichment", index=False)
            pd.DataFrame().to_excel(writer, sheet_name = "6.Protein_Summary", index=False)
            output.seek(0)
            return output.getvalue()
        
        # 3. Interactions

        df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
        df_interactions.to_excel(writer, sheet_name = "3.Interactions", header = True, index=False)

        # 4. Pathways

        df_pathways = results.get("df_pathways", pd.DataFrame()).copy()
        df_groupedpathways = results.get("df_groupedpathways", pd.DataFrame()).copy()
        
        num_row = 0

        if not df_groupedpathways.empty:
            df_groupedpathways.to_excel(writer, sheet_name = "4.Pathways", index=False, startrow=num_row)
            num_row += len(df_groupedpathways) +3
        else:
            pd.DataFrame(columns=["pathway", "n_proteins", "n_compounds", "compounds"]).to_excel(
                writer,
                sheet_name="4.Pathways",
                index=False,
                startrow=num_row,
            )
            num_row += 3

        if not df_pathways.empty and not df_groupedpathways.empty:
            for _, row in df_groupedpathways.iterrows(): #TODO check this line
                pathway = row["pathway"]
                pd.DataFrame({"Pathway": [pathway]}).to_excel(
                    writer,
                    sheet_name = "4.Pathways",
                    index=False,
                    header=False,
                    startrow=num_row,
                )
                num_row+=1

                detailed_table = app.pathways.group_compounds(df_pathways, pathway)
                if not detailed_table.empty:
                    detail_cols = ["uniprot_accession", "protein_name"]
                    if "count" in detailed_table.columns:
                        detail_cols.append("count")
                    detail_cols.append("compounds")

                    detailed_table = detailed_table[detail_cols]
                    detailed_table.to_excel(
                        writer,
                        sheet_name="4.Pathways",
                        index=False,
                        startrow=num_row,
                    )
                    num_row += len(detailed_table) + 3
                else:
                    pd.DataFrame(columns=[
                        "uniprot_accession",
                        "protein_name",
                        "count",
                        "compounds"
                    ]).to_excel(writer, sheet_name = "4.Pathways", index=False, startrow=num_row)
                    num_row += 3

        else:
            pd.DataFrame(columns=[
                        "uniprot_accession",
                        "protein_name",
                        "count",
                        "compounds"
                    ]).to_excel(writer, sheet_name = "4.Pathways", index=False, startrow=num_row)
                    

        # 5. GO Enrichment            

        df_go_bp_grouped = results.get("df_go_bp_grouped", pd.DataFrame()).copy()
        df_go_mf_grouped = results.get("df_go_mf_grouped", pd.DataFrame()).copy()
        df_go_cc_grouped = results.get("df_go_cc_grouped", pd.DataFrame()).copy()
        
        df_go_bp = results.get("df_go_bp", pd.DataFrame()).copy()
        df_go_mf = results.get("df_go_mf", pd.DataFrame()).copy()
        df_go_cc = results.get("df_go_cc", pd.DataFrame()).copy()

        num_row = 0
        df_confirmation = pd.concat(
            [
                df_interactions[["uniprot_accession", "compound"]]
                if not df_interactions.empty and {"uniprot_accession", "compound"}.issubset(df_interactions.columns)
                else pd.DataFrame(columns=["uniprot_accession", "compound"]),
                df_pathways[["uniprot_accession", "compound"]]
                if not df_pathways.empty and {"uniprot_accession", "compound"}.issubset(df_pathways.columns)
                else pd.DataFrame(columns=["uniprot_accession", "compound"])
            ],
            ignore_index=True,
        )

        def write_go_aspect(sheet, title, grouped_df, detail_df, num_row):
            pd.DataFrame({"GO Aspect": [title]}).to_excel(writer, sheet_name = sheet, index=False, header=True, startrow = num_row)
            num_row += 2

            if not grouped_df.empty:
                grouped_df.to_excel(writer, sheet_name = sheet, index=False, startrow=num_row)
                num_row += len(grouped_df) + 2
            else:
                pd.DataFrame(columns=["go_name", "go_id", "n_compounds", "n_proteins"]).to_excel(
                    writer, sheet_name = sheet, index=False, header=True, startrow=num_row)
                num_row += 3
            
            if not grouped_df.empty and not detail_df.empty:
                for _, row in grouped_df.iterrows():
                    go_name = row["go_name"]
                    go_id = row["go_id"]

                    header_df = pd.DataFrame({"GO Term": [f"{go_name} ({go_id})"]})
                    header_df.to_excel(writer, sheet_name = sheet, index=False, header=False, startrow=num_row)
                    num_row += 1

                    subset = app.proteins.count_protgoaspect(detail_df, go_name, go_id, df_confirmation)
                    if not subset.empty:
                        detail_cols = ["symbol", "uniprot_accession"]
                        if "count" in subset.columns:
                            detail_cols.append("count")

                        detail_cols.append("compounds")

                        subset = subset[detail_cols]
                        subset.to_excel(writer, sheet_name = sheet, index=False, startrow = num_row)
                        num_row += len(subset) +3
                    else:
                        pd.DataFrame(columns=["symbol", "uniprot_accession", "compounds"]).to_excel(
                            writer, sheet_name=sheet, index=False, startrow=num_row,
                        )
                        num_row += 3

            return num_row
        
        num_row = write_go_aspect("5.GO_Enrichment", "Biological Process", df_go_bp_grouped, df_go_bp, num_row)
        num_row = write_go_aspect("5.GO_Enrichment", "Molecular Function", df_go_mf_grouped, df_go_mf, num_row)
        num_row = write_go_aspect("5.GO_Enrichment", "Cellular Component", df_go_cc_grouped, df_go_cc, num_row)

        df_proteinsummary = results.get("final_summaryGO", pd.DataFrame()).copy()
        df_proteinsummary.to_excel(writer, sheet_name = "6.Protein_Summary", header = True, index=False)

    output.seek(0)
    return output.getvalue()




# Hacer que este botón NO modifique los SMILES adicionales
st.button("Save input", on_click=save_input, type="primary", width='stretch')

submitted_smiles = st.session_state.get("submitted_smiles", [])

if "submitted_smiles":
    st.subheader("Saved SMILES & email")
    st.write(f"Email: {st.session_state['submitted_email']}")
    st.write(f"Number of SMILES code stored: {len(submitted_smiles)}")
    

    # Recalculate names if SMILES are changed TODO what is tuple
    if st.session_state.get("compound_results_source") != tuple(submitted_smiles):
        st.session_state["compound_results"] = app.chem.identify_compounds(submitted_smiles)
        # This is an indicator to know whether the SMILES changed and the results need to be updated
        st.session_state["compound_results_source"] = tuple(submitted_smiles)

    # Enseña una tab que clicas y te enseña los SMILES
    with st.expander("Show SMILES list"):
        compound_df = st.session_state.get("compound_results", pd.DataFrame()).copy()

        if not compound_df.empty:
            # Make a dictionary with the input_smiles as keys and the compound_name as values.
            # It pairs each SMILES with its compound name
            # It is a dictionary so that it can be paired
            #name = dict(zip(compound_df["input_smiles"], compound_df["compound_name"]))
            
            display_df = compound_df.copy()
            display_df["compound_name"] = display_df.apply(lambda row: compound_hyperlink(row["cid"], row["compound_name"]), axis=1)

            display_df = display_df[["smiles", "compound_name", "molecular_formula", "molecular_weight", "status"]]
            st.markdown(display_df.to_html(escape=False, index=False), unsafe_allow_html=True)
        else:
            st.info("No valid SMILES code entered yet.")
else:
    st.session_state["compound_results"] = pd.DataFrame(columns=[
        "smiles",
        "compound_name",
        "cid",
        "molecular_formula",
        "molecular_weight",
        "status",
    ])
    st.session_state["compound_results_source"] = tuple()

# BEFORE RUNNING PIPELINE -> create placeholders

status_box = st.empty()
progress_bar = st.empty()
compound_box = st.empty()
interactions_box = st.empty()
pathway_box = st.empty()
summary_box = st.empty()

st.button("Run Analysis", on_click=run_analysis, type="primary", width='stretch')
# Printing the boolean result of the analysis
if st.session_state.get("analysis_ready", False):
    st.success("Analysis completed successfully!")
if st.session_state.get("run_error", ""):
    st.error(st.session_state["run_error"])

results = st.session_state.get("results", None)
if results:
    excel = download_excel_analysis()

    st.download_button(
        label="Download analysis results as Excel",
        data=excel,
        file_name="analysis_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

    st.subheader("Analysis Overview")
    col1, col2, col3 = st.columns(3)
    col1.metric("Compounds", len(results.get("compound_names", [])))
    col2.metric("Interactions", len(results.get("df_interactions", pd.DataFrame())))
    col3.metric("Pathways", len(results.get("df_groupedpathways", pd.DataFrame())))
    
    df_go_bp_grouped = results.get("df_go_bp_grouped", pd.DataFrame())
    df_go_mf_grouped = results.get("df_go_mf_grouped", pd.DataFrame())
    df_go_cc_grouped = results.get("df_go_cc_grouped", pd.DataFrame())

    bp_count = df_go_bp_grouped["go_id"].nunique() if not df_go_bp_grouped.empty else 0
    mf_count = df_go_mf_grouped["go_id"].nunique() if not df_go_mf_grouped.empty else 0
    cc_count = df_go_cc_grouped["go_id"].nunique() if not df_go_cc_grouped.empty else 0
    
    st.markdown('#### GO enrichment')
    go_total, go1, go2, go3 = st.columns(4)
    go_total.metric("GO terms", bp_count + mf_count + cc_count)
    go1.metric("Biological process", bp_count)
    go2.metric("Molecular function", mf_count)
    go3.metric("Cellular component", cc_count)

    
    st.subheader("Analysis results")

    st.markdown("## Compounds")
    displaycomp_df = st.session_state.get("compound_results", pd.DataFrame()).copy()
    if not displaycomp_df.empty:
        displaycomp_df["compound_name"] = displaycomp_df.apply(
            lambda row: compound_hyperlink(row["cid"], row["compound_name"]),
            axis = 1,
        )
        displaycomp_df = displaycomp_df[["smiles", "compound_name", "cid", "molecular_formula", "molecular_weight", "status"]]
        st.markdown(displaycomp_df.to_html(escape=False, index=False), unsafe_allow_html=True)
    else:
        st.info("No compounds identified")

    st.divider()

    st.markdown("## Compound-Protein Interactions")

    df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
    compound_results = st.session_state.get("compound_results", pd.DataFrame())

    if not df_interactions.empty:
        display_df = df_interactions.copy()
        display_df["compound"] = display_df.apply(
            lambda row: compound_hyperlink(row["cid"], row["compound"]), 
            axis=1
        )
        display_df["uniprot_accession"] = display_df["uniprot_accession"].apply(
            lambda acc: uniprot_hyperlink(acc))
        
        display_df["geneid"] = display_df["geneid"].apply(
            lambda geneid: geneid_hyperlink(geneid)
        )
       
        display_df = display_df[
            ["compound", "geneid", "symbol", "uniprot_accession", "description"]
        ]

        st.markdown(display_df.to_html(escape=False, index=False), unsafe_allow_html=True)
    else:
        st.info("No valid SMILES code entered yet.")

    st.divider()

    st.markdown("## Pathways")

    df_groupedpathways = results.get("df_groupedpathways", pd.DataFrame()).copy()
    df_pathways = results.get("df_pathways", pd.DataFrame()).copy()

    if not df_groupedpathways.empty and not df_pathways.empty:
        # 1) Pathways table with the pathway as a selector column in plain text
        event = st.dataframe(
            df_groupedpathways[["pathway","n_compounds", "n_proteins","compounds"]],
            width="stretch",
            on_select="rerun",
            selection_mode="single-row",
            key="pathway_table",
        )

        selected_rows = event.selection["rows"]
        # 2) When clicking on a pathway, it shows pathway title + all associated proteins 
        if selected_rows:
            selected_index = selected_rows[0]
            selected_pathway = df_groupedpathways.iloc[selected_index]["pathway"]
            
            # Clickable title pathway with hyperlink to PubChem
            st.markdown(f"### Proteins in pathway: {pathway_hyperlink(selected_pathway)}",
                        unsafe_allow_html=True,
            )

            pathway_details = app.pathways.group_compounds(df_pathways, selected_pathway)

            if not pathway_details.empty:
                display_details = pathway_details.copy()

                display_details["uniprot_accession"] = display_details["uniprot_accession"].apply(
                    uniprot_hyperlink)
                
                display_details["compounds"] = display_details["compounds"].apply(
                    lambda row: compounds_text_with_links(row, compound_results)
                )
                # Show details table
                display_details = display_details[
                    ["uniprot_accession", "protein_name", "count","compounds"]
                ]

                st.markdown(display_details.to_html(escape=False, index=False), unsafe_allow_html=True)
            else:
                st.info("No proteins found in this pathway.")
        else:
            st.info("Select a pathway to see the details.")

    st.divider()

    st.markdown("## GO enrichment")
    compound_results = st.session_state.get("compound_results", pd.DataFrame())
    df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
    df_pathways = results.get("df_pathways", pd.DataFrame()).copy()

    df_confirmation = pd.concat(
        [
            df_interactions[["uniprot_accession", "compound"]],
            df_pathways[["uniprot_accession", "compound"]],
        ],
        ignore_index=True,
    )

    #BIOLOGICAL PROCESS
    with st.expander("Biological process"):
        if (
            "df_go_bp_grouped" in results and not results["df_go_bp_grouped"].empty and "df_go_bp" in results and not results["df_go_bp"].empty
        ):
            df_go_bp_grouped = results["df_go_bp_grouped"].copy()
            df_go_bp = results["df_go_bp"].copy()
            
            event_bp = st.dataframe(
                df_go_bp_grouped[["go_name", "go_id", "n_compounds", "n_proteins"]],
                width="stretch",
                on_select="rerun",
                selection_mode="single-row",
                key="go_bp_table",
            )

            selected_rows_bp = event_bp.selection["rows"]
            if selected_rows_bp:
                selected_index_bp = selected_rows_bp[0]
                selected_go_name_bp = df_go_bp_grouped.iloc[selected_index_bp]["go_name"]
                selected_go_id_bp = df_go_bp_grouped.iloc[selected_index_bp]["go_id"]
                
            # Clickable GO term with hyperlink to QuickGO
                st.markdown(f"### Proteins in GO term: {goterm_hyperlink(selected_go_id_bp, f'{selected_go_name_bp} ({selected_go_id_bp})')}",
                        unsafe_allow_html=True)
                go_details_bp = app.proteins.count_protgoaspect(df_go_bp, selected_go_name_bp, selected_go_id_bp, df_confirmation)

                if not go_details_bp.empty:
                    display_godetails_bp = go_details_bp.copy()
                    display_godetails_bp["uniprot_accession"] = display_godetails_bp["uniprot_accession"].apply(
                        uniprot_hyperlink)
                    display_godetails_bp["compounds"] = display_godetails_bp["compounds"].apply(
                        lambda row: compounds_text_with_links(row, compound_results),
                        )
                    # Show details table
                    display_godetails_bp = display_godetails_bp[
                        ["symbol", "uniprot_accession", "count","compounds"]
                    ]
                    st.markdown(display_godetails_bp.to_html(escape=False, index=False), unsafe_allow_html=True)
                else:
                    st.info("No biological process GO terms found.")
            else:
                st.info("Select a GO term to see the details.")

    #MOLECULAR FUNCTION
    with st.expander("Molecular function"):
        if (
            "df_go_mf_grouped" in results and not results["df_go_mf_grouped"].empty and "df_go_mf" in results and not results["df_go_mf"].empty
        ):
            df_go_mf_grouped = results["df_go_mf_grouped"].copy()
            df_go_mf = results["df_go_mf"].copy()
            event_mf = st.dataframe(
                df_go_mf_grouped[["go_name", "go_id", "n_compounds", "n_proteins"]],
                width="stretch",
                on_select="rerun",
                selection_mode="single-row",
                key="go_mf_table",
            )

            selected_rows_mf = event_mf.selection["rows"]
            if selected_rows_mf:
                selected_index_mf = selected_rows_mf[0]
                selected_go_name_mf = df_go_mf_grouped.iloc[selected_index_mf]["go_name"]
                selected_go_id_mf = df_go_mf_grouped.iloc[selected_index_mf]["go_id"]
                df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
                df_pathways = results.get("df_pathways", pd.DataFrame()).copy()

                df_confirmation = pd.concat(
                    [
                        df_interactions[["uniprot_accession", "compound"]],
                        df_pathways[["uniprot_accession", "compound"]],
                    ],
                    ignore_index=True,
                )

                st.markdown(f"### Proteins in GO term: {goterm_hyperlink(selected_go_id_mf, f'{selected_go_name_mf} ({selected_go_id_mf})')}",
                        unsafe_allow_html=True)
                go_details_mf = app.proteins.count_protgoaspect(df_go_mf, selected_go_name_mf, selected_go_id_mf, df_confirmation)
                if not go_details_mf.empty:
                    display_godetails_mf = go_details_mf.copy()
                    display_godetails_mf["uniprot_accession"] = display_godetails_mf["uniprot_accession"].apply(
                        uniprot_hyperlink)
                    display_godetails_mf["compounds"] = display_godetails_mf["compounds"].apply(
                        lambda row: compounds_text_with_links(row, compound_results),
                        )
                    # Show details table
                    display_godetails_mf = display_godetails_mf[
                        ["symbol", "uniprot_accession", "count","compounds"]
                    ]
                    st.markdown(display_godetails_mf.to_html(escape=False, index=False), unsafe_allow_html=True)
                else:
                    st.info("No molecular function GO terms found.")
            else:
                st.info("Select a GO term to see the details.")

    
    #CELLULAR COMPONENT
    with st.expander("Cellular component"):
        if (
            "df_go_cc_grouped" in results and not results["df_go_cc_grouped"].empty and "df_go_cc" in results and not results["df_go_cc"].empty
        ):
            df_go_cc_grouped = results["df_go_cc_grouped"].copy()
            # Detail table after clicking one GO term in the grouped table
            df_go_cc = results["df_go_cc"].copy()
            event_cc = st.dataframe(
                df_go_cc_grouped[["go_name", "go_id", "n_compounds", "n_proteins"]],
                width="stretch",
                on_select="rerun",
                selection_mode="single-row",
                key="go_cc_table",
            )

            selected_rows_cc = event_cc.selection["rows"]
            if selected_rows_cc:
                selected_idex_cc = selected_rows_cc[0]
                selected_go_name_cc = df_go_cc_grouped.iloc[selected_idex_cc]["go_name"]
                selected_go_id_cc = df_go_cc_grouped.iloc[selected_idex_cc]["go_id"]
                df_interactions = results.get("df_interactions", pd.DataFrame()).copy()
                df_pathways = results.get("df_pathways", pd.DataFrame()).copy()

                df_confirmation = pd.concat(
                    [
                        df_interactions[["uniprot_accession", "compound"]],
                        df_pathways[["uniprot_accession", "compound"]],
                    ],
                    ignore_index=True,
                )

                st.markdown(f"### Proteins in GO term: {goterm_hyperlink(selected_go_id_cc, f'{selected_go_name_cc} ({selected_go_id_cc})')}", unsafe_allow_html=True)
                go_details_cc = app.proteins.count_protgoaspect(df_go_cc, selected_go_name_cc, selected_go_id_cc, df_confirmation)
                if not go_details_cc.empty:
                    display_godetails_cc = go_details_cc.copy()
                    display_godetails_cc["uniprot_accession"] = display_godetails_cc["uniprot_accession"].apply(
                        uniprot_hyperlink)
                    display_godetails_cc["compounds"] = display_godetails_cc["compounds"].apply(
                        lambda row: compounds_text_with_links(row, compound_results),
                        )
                    # Show details table
                    display_godetails_cc = display_godetails_cc[
                        ["symbol", "uniprot_accession", "count","compounds"]
                    ]
                    st.markdown(display_godetails_cc.to_html(escape=False, index=False), unsafe_allow_html=True)
                else:
                    st.info("No cellular component GO terms found.")
            else:
                st.info("Select a GO term to see the details.")
    
    st.divider()

    st.markdown("## Protein Summary")
    
    df_proteinsummary = results.get("final_summary",pd.DataFrame()).copy()
    compound_results = st.session_state.get("compound_results", pd.DataFrame())

    if not df_proteinsummary.empty:
        display_dfsum = df_proteinsummary.copy()
        display_dfsum["uniprot_accession"] = display_dfsum["uniprot_accession"].apply(
            uniprot_hyperlink)  
        
        display_dfsum["compounds"] = display_dfsum["compounds"].apply(
            lambda row: compounds_text_with_links(row, compound_results),
        )

        display_dfsum["pathways"]  = display_dfsum["pathways"].apply(
            pathway_text_with_links
        )

        display_dfsum["pathway_compounds"] = display_dfsum["pathway_compounds"].apply(
            lambda row: compounds_text_with_links(row, compound_results),
        )

        display_dfsum = display_dfsum[
            ["uniprot_accession", "total_count", "n_compounds", "compounds", "symbol", "n_pathways", "pathways", "pathway_compounds", "source"]
        ]

        st.markdown(display_dfsum.to_html(escape=False, index=False), unsafe_allow_html=True)

    else:   
        st.info("No proteins found in the summary.")
