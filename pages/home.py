import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
from app.pipeline import run_full_pipeline
import app.chem
import app.downloads
import html
import re
from urllib.parse import quote

# For further load into the app
example_smiles = [
    "CC(C[N+](C)(C)C)OC(=O)N",
]

TAXONOMY_OPTIONS = {
    "Homo sapiens (human)": "9606",
    "Mus musculus (mouse)": "10090",
    "Rattus norvegicus (rat)": "10116",
    "Danio rerio (zebrafish)": "7955",
    "Cavia porcellus (guinea pig)": "10141",
}

# Triple quotes allow multiple lines
# st.markdown() : for text you want formatted with Markdown
# st.write() is more general-purpose and can display text, dataFrames, lists and other Python objects
st.title("Drug Discovery Analysis")
st.markdown(
    """
    Enter one or more **SMILES codes** to prepare the analysis. 
    Use **one SMILES code per line**
    """)


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
        "compound_results_source",
        "selected_tax_ids",
        "selected_taxonomy_labels",
        "input_warning",
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

if "input_warning" not in st.session_state:
    st.session_state["input_warning"] = ""

# This method is for clearing the results of the analysis only when the user modifies the SMILES input (should take into account the deletion of a SMILES code)
def clear_analysis_results():
    for key in ["results", "analysis_ready", "run_error", "prepared_excel", "prepared_excel_selection"]:
        if key in st.session_state:
            del st.session_state[key]

# These methods modify the st.session_state of the widgets that were already initialized with a key
# It DELETES the analysis results previously made
def load_example():
    if "submitted_smiles" not in st.session_state:
        st.session_state["submitted_smiles"] = []

    st.session_state["smiles_input"] = "\n".join(example_smiles)
    clear_analysis_results()

def check_email(email):
    email = str(email).strip()
    return bool(re.match(r"^[\w\.-]+@[\w\.-]+\.\w+$", email))

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

    canonical, error = app.chem.validate_smiles(value)

    if error:
        st.session_state["input_warning"] = f"Could not add invalid SMILES: {value}"
        st.session_state["new_smiles_input"] = ""
        return 

    existing_canonical = set()

    for smiles in st.session_state["submitted_smiles"]:
        existing_canonical_smiles, existing_error = app.chem.validate_smiles(smiles)
        if existing_canonical and not existing_error:
            existing_canonical.add(existing_canonical_smiles)
    
    if canonical in existing_canonical:
        st.session_state["input_warning"] = f"Duplicate SMILES ignores: {value}"
        st.session_state["new_smiles_input"] = ""
        return

    st.session_state["submitted_smiles"].append(canonical)
    st.session_state["input_warning"] = ""
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
    existing_smiles = st.session_state.get("submitted_smiles", [])

    checked_smiles = app.chem.check_smiles_box(st.session_state.get("smiles_input", ""), existing_smiles=existing_smiles)

    cleaned_smiles = checked_smiles["valid"]
    invalid_smiles = checked_smiles["invalid"]
    duplicate_smiles = checked_smiles["duplicates"]
    
    if not cleaned_smiles and not existing_smiles:
        if invalid_smiles:
            st.session_state["run_error"] = (
                "No valid SMILES were added. Invalid entried: "+", ".join(invalid_smiles)
            )
        else:
            st.session_state["run_error"] = "Please enter at least one SMILES code"
        return

    for smiles in cleaned_smiles:
        st.session_state["submitted_smiles"].append(smiles)
    
    warning_messages = []

    if invalid_smiles:
        warning_messages.append("Invalid SMILES ignored: "+", ".join(invalid_smiles))

    if duplicate_smiles:
        warning_messages.append("Duplicate SMILES ignored: "+", ".join(duplicate_smiles))

    st.session_state["input_warning"] = '\n\n'.join(warning_messages)

    email_input = st.session_state.get("email_input", "").strip()
    saved_email = st.session_state.get("submitted_email", "").strip()

    email = email_input or saved_email

    if not email:
        st.session_state["run_error"] = "Please enter a valid email"
        return
    
    if not check_email(email):
        st.session_state["run_error"] = "Please enter a valid email address"
        return
    
    taxonomy = st.session_state.get("selected_tax_ids", [])
    if not taxonomy:
        st.session_state["run_error"] = "Please select at least one protein taxonomy"
        return
    
    st.session_state["submitted_email"] = email
    st.session_state["email_input"] = email
    st.session_state["smiles_input"] = ""
    st.session_state["run_error"] = ""
    st.session_state["taxonomy"] = taxonomy
    clear_analysis_results()

# TAXONOMY SELECTOR
selected_taxonomy_labels = st.multiselect(
    "Select which protein taxonomies to include",
    options = list(TAXONOMY_OPTIONS.keys()),
    default = ["Homo sapiens (human)"],
    key = "selected_taxonomy_labels",
    help = "Only proteins from the selected organisms will be included in the analysis.",
    on_change=clear_analysis_results,
)

st.session_state["selected_tax_ids"] = [
    TAXONOMY_OPTIONS[label] for label in selected_taxonomy_labels
]

def run_analysis():
    smiles_codes = st.session_state.get("submitted_smiles", [])
    email = st.session_state.get("submitted_email", "")
    selected_tax_ids = st.session_state.get("selected_tax_ids", ["9606"])
    
    if not selected_tax_ids:
        st.session_state["run_error"] = "Please select at least one protein taxonomy before running the analysis."
        st.session_state["analysis_ready"] = False
        return

    if not email:
        st.session_state["run_error"] = "Please enter a valid email before running the analysis."
        st.session_state["analysis_ready"] = False
        return
    if not check_email(email):
        st.session_state["run_error"] = "Please enter a valid email address before running the analysis."
        st.session_state["analysis_ready"] = False
        return
    if not smiles_codes:
        st.session_state["run_error"] = "Please enter at least one SMILES code before running the analysis."
        st.session_state["analysis_ready"] = False
        return
    try:
        with st.spinner("Running analysis... This may take a few minutes..."):
            results = run_full_pipeline(smiles_codes, email, selected_tax_ids=selected_tax_ids, ui = {
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

# CODE INJECTION RISK MANAGEMENT -------
CID_re = re.compile(r"^\d+$")
GENEID_re=re.compile(r"^\d+$")
UNIPROT_re = re.compile(r"^[A-Z0-9]+(?:-[0-9]+)?$")
GO_re = re.compile(r"^GO:\d{7}$")
PATHWAY_re = re.compile(r"^[A-Za-z0-9_.:-]+$")

def safe_text(value):
    return html.escape("" if value is None else str(value), quote=True)

# It is built only from validated IDs
def safe_anchor(url, label):
    safe_url = html.escape(url, quote=True)
    safe_label = safe_text(label)
    return f'<a href="{safe_url}" target="_blank" rel="noopener noreferrer">{safe_label}</a>'

def compound_hyperlink(cid, compound_name):
    cid = str(cid).strip()
    if CID_re.fullmatch(cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{int(cid)}"
        return safe_anchor(url, compound_name)
    return safe_text(compound_name)
          
def uniprot_hyperlink(accession, ):
    accession = str(accession).strip()
    if UNIPROT_re.fullmatch(accession):
        url = f"https://www.uniprot.org/uniprotkb/{accession}/entry"
        return safe_anchor(url, accession)
    return safe_text(accession)
   
def pathway_hyperlink(pwacc, name = None):
    pwacc = str(pwacc).strip()
    name = str(name).strip() if name is not None else pwacc

    if PATHWAY_re.fullmatch(pwacc):
        url = f"https://pubchem.ncbi.nlm.nih.gov/pathway/{quote(pwacc, safe='')}"
        return safe_anchor(url, name)
    return safe_text(name) 

def goterm_hyperlink(go_id, go_name=None):
    go_id = str(go_id).strip()
    go_name = str(go_name).strip() if go_name else go_id
    if GO_re.fullmatch(go_id):
        url = f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"
        return safe_anchor(url, go_name)
    return safe_text(go_name)  

def geneid_hyperlink(geneid):
    geneid = str(geneid).strip()
    if GENEID_re.fullmatch(geneid):
        url = f"https://www.ncbi.nlm.nih.gov/gene/{geneid}"
        return safe_anchor(url, geneid)
    return safe_text(geneid)
    

def compounds_text_with_links(compounds_text, compound_results):
    if pd.isna(compounds_text):
        return ""
    
    if compound_results is None or compound_results.empty:
        return safe_text(compounds_text)
    
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
        cid_str = str(cid).strip() if pd.notna(cid) else ""
        if CID_re.fullmatch(cid_str):
            url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_str}"
            links.append(safe_anchor(url, compound))
        else:
            links.append(safe_text(compound))

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

# TABLE DISPLAY --------------

def show_interactive_table(df, table_id, height=None):
    if df is None or df.empty:
        st.info("No data available.")
        return

    # Optional: remove completely empty rows just in case
    df = df.dropna(how="all").copy()

    html_table = df.to_html(
        escape=False,
        index=False,
        table_id=table_id,
    )

    html_code = f"""
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

    <style>
    body {{
        margin: 0;
        padding: 0;
        background-color: white;
    }}

    .dataTables_wrapper, .dataTables_wrapper * {{
        font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
        color: #31333f !important;
        font-size: 15px !important;
    }}

    #{table_id} {{
        width: 100% !important;
        border-collapse: collapse !important;
        background-color: white !important;
    }}

    #{table_id},
    #{table_id} thead,
    #{table_id} tbody,
    #{table_id} tr,
    #{table_id} td,
    #{table_id} th {{
        background-color: white !important;
        color: #31333f !important;
        border-color: #ccc !important;
    }}

    #{table_id} thead th {{
        font-weight: 600 !important;
        border: 1px solid #999 !important;
        padding: 10px 12px !important;
        text-align: left !important;
        white-space: nowrap !important;
    }}

    #{table_id} tbody td {{
        border: 1px solid #999 !important;
        padding: 9px 12px !important;
        vertical-align: top !important;
        white-space: normal !important;
    }}

    #{table_id} a {{
        color: #0068c9 !important;
        text-decoration: underline !important;
        font-weight: 500 !important;
    }}

    /* Search row */
    #{table_id} thead tr.filters th {{
        padding: 6px 10px !important;
        font-weight: normal !important;
    }}

    #{table_id} thead tr.filters input {{
        width: 95% !important;
        box-sizing: border-box !important;
        background-color: white !important;
        color: #31333f !important;
        font-size: 14px !important;
        border: 1px solid #ccc !important;
        padding: 5px 7px !important;
    }}

    .dataTables_length select {{
        background-color: white !important;
        color: #31333f !important;
        border: 1px solid #ccc !important;
        font-size: 14px !important;
        padding: 4px 6px !important;
    }}

    .dataTables_filter input {{
        background-color: white !important;
        color: #31333f !important;
        border: 1px solid #ccc !important;
        font-size: 14px !important;
        padding: 5px 7px !important;
    }}

    .dataTables_info {{
        color: #31333f !important;
        font-size: 14px !important;
    }}

    .dataTables_paginate a {{
        color: #31333f !important;
        font-size: 14px !important;
        font-weight: normal !important;
        padding: 6px 12px !important;
        border-radius: 6px !important;
        text-decoration: none !important;
        margin: 0 2px !important;
    }}

    .dataTables_paginate a.current {{
        background-color: #f0f0f0 !important;
        font-weight: bold !important;
        border: 1px solid #aaa !important;
    }}

    .dataTables_paginate a:hover {{
        background-color: #e3e3e3 !important;
    }}

    .dataTables_scrollBody thead tr {{
        height: 0 !important;
    }}

    .dataTables_scrollBody thead th {{
        height: 0 !important;
        padding-top: 0 !important;
        padding-bottom: 0 !important;
        border-top: none !important;
        border-bottom: none !important;
        line-height: 0 !important;
        visibility: hidden !important;
    }}
    </style>

    <script>
    $(document).ready(function() {{
        var tableSelector = '#{table_id}';

        // Create second header row for column filters
        $(tableSelector + ' thead tr')
            .clone(false)
            .addClass('filters')
            .appendTo(tableSelector + ' thead');

        // Replace cloned header titles with search boxes BEFORE DataTables starts
        $(tableSelector + ' thead tr.filters th').each(function(i) {{
            var title = $(this).text();

            $(this).html(
                '<input type="text" placeholder="Search ' + title + '" />'
            );
        }});

        var table = $(tableSelector).DataTable({{
            pageLength: 10,
            lengthMenu: [[10, 20, 50, 100, -1], [10, 20, 50, 100, "All"]],
            paging: true,
            searching: true,
            ordering: true,
            orderCellsTop: true,
            autoWidth: false,
            fixedHeader: false
        }});

        // Column-specific search
        $(tableSelector + ' thead tr.filters input').on('keyup change clear', function(e) {{
            e.stopPropagation();

            var colIdx = $(this).parent().index();

            if (table.column(colIdx).search() !== this.value) {{
                table
                    .column(colIdx)
                    .search(this.value)
                    .draw();
            }}
        }});

        // Prevent sorting when clicking inside search boxes
        $(tableSelector + ' thead tr.filters input').on('click', function(e) {{
            e.stopPropagation();
        }});
    }});
    </script>

    <div style="overflow-x:auto; margin-bottom: 0px;">
        {html_table}
    </div>
    """

    if height is None:
        height = min(900, 300 + min(len(df), 20) * 45)

    components.html(
        html_code,
        height=height,
        scrolling=True,
    )

# Hacer que este botón NO modifique los SMILES adicionales
st.button("Save input", on_click=save_input, type="primary", width='stretch')

if st.session_state.get("input_warning"):
    st.warning(st.session_state["input_warning"])

submitted_smiles = st.session_state.get("submitted_smiles", [])

if submitted_smiles:
    st.subheader("Saved SMILES & email")
    st.write(f"Email: {st.session_state['submitted_email']}")
    st.write(f"Number of SMILES code stored: {len(submitted_smiles)}")
    

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

    st.markdown("### Download selected result tables")
    selected_pages = []
    col1, col2 = st.columns(2)

    for i, (section_key, section_label) in enumerate(app.downloads.page_selection.items()):
        column = col1 if i % 2 == 0 else col2 # design alternating between column 1 and column 2
        # Checkbox
        with column:
            checked = st.checkbox(
                section_label,
                value = True,
                key = f"download_{section_key}",
            )
        if checked:
            selected_pages.append(section_key)

    if selected_pages:
        if st.button("Prepare Excel file", type = "primary"):
            with st.spinner("Preparing excel file for download..."):
                st.session_state["prepared_excel"] = app.downloads.download_excel_analysis(selected_pages)
                st.session_state["prepared_excel_selection"] = tuple(selected_pages)
            # Button
        if "prepared_excel" in st.session_state:
            st.download_button(
                label = "Download selected tables",
                data = st.session_state["prepared_excel"],
                file_name = "Analysis_results.xlsx",
                mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                on_click = "ignore",
            )
    else:
        st.info("Select at least one table to download.")

    st.subheader("Analysis results")

# COMPOUNDS -----------------
    st.markdown("## Compounds")
    displaycomp_df = st.session_state.get("compound_results", pd.DataFrame()).copy()
    if not displaycomp_df.empty:
        displaycomp_df["compound_name"] = displaycomp_df.apply(
            lambda row: compound_hyperlink(row["cid"], row["compound_name"]),
            axis = 1,
        )
        displaycomp_df = displaycomp_df[["smiles", "compound_name", "cid", "molecular_formula", "molecular_weight", "status"]]
        show_interactive_table(
            displaycomp_df,
            table_id="compounds_table",
        )
    else:
        st.info("No compounds identified")

    st.divider()

# COMPOUND-PROTEIN INTERACTIONS -------------
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
            ["compound", "geneid", "symbol", "protein_name", "uniprot_accession", "description", "taxid", "taxname",]
        ]
        show_interactive_table(
            display_df,
            table_id="interactionsTable",
        )
    else:
        st.info("No valid SMILES code entered yet.")

    st.divider()

# PATHWAYS -------------
    st.markdown("## Pathways")

    df_groupedpathways = results.get("df_groupedpathways", pd.DataFrame()).copy()
    df_pathways = results.get("df_pathways", pd.DataFrame()).copy()

    if not df_groupedpathways.empty and not df_pathways.empty:
        # 1) Pathways table with the pathway as a selector column in plain text
        for col in ["pathway_name", "proteins", "taxid", "taxname"]:
            if col not in df_groupedpathways.columns:
                df_groupedpathways[col] = ""
        
        event = st.dataframe(
            df_groupedpathways[["pathway", "pathway_name", "n_compounds", "n_proteins","compounds"]],
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
            selected_pathway_name = df_groupedpathways.iloc[selected_index].get("pathway_name", "")
            
            if pd.notna(selected_pathway_name) and str(selected_pathway_name).strip():
                pathway_label = f"{selected_pathway_name} ({selected_pathway})"
            else:
                pathway_label = selected_pathway
            # Clickable title pathway with hyperlink to PubChem
            st.markdown(f"### Proteins in pathway: {pathway_hyperlink(selected_pathway, pathway_label)}",
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
                    ["uniprot_accession", "protein_name", "symbol", "count","compounds", "taxid", "taxname"]
                ]
                show_interactive_table(
                    display_details,
                    table_id="pathwayDetailsTable",
                )
            else:
                st.info("No proteins found in this pathway.")
        else:
            st.info("Select a pathway to see the details.")

    st.divider()

# GO ENRICHMENT -------------
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
                    show_interactive_table(
                        display_godetails_bp,
                        table_id="goDetailsBPTable",
                    )
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
                    show_interactive_table(
                        display_godetails_mf,
                        table_id="goDetailsMFTable",
                    )
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
                    show_interactive_table(
                        display_godetails_cc,
                        table_id="goDetailsCCTable",
                    )
                else:
                    st.info("No cellular component GO terms found.")
            else:
                st.info("Select a GO term to see the details.")
    
    st.divider()

# PROTEIN SUMMARY -------------
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

        display_dfsum = display_dfsum[[
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
        ]]

        show_interactive_table(
            display_dfsum,
            table_id="proteinSummaryTable",
        )

    else:   
        st.info("No proteins found in the summary.")