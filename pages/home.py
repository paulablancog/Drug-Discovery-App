import streamlit as st
import pandas as pd
from app.pipeline import run_full_pipeline

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
    
    if not cleaned_smiles and "submitted_smiles" not in st.session_state:
        st.session_state["run_error"] = "Please enter at least one SMILES code"
        return

    for smiles in cleaned_smiles:
        if smiles not in st.session_state.get("submitted_smiles", []):
            st.session_state["submitted_smiles"].append(smiles)
    

    email = st.session_state.get("email_input", "").strip()

    if not email:
        st.error("Please enter an email.")
        return
    
    st.session_state["submitted_email"] = email
    st.session_state["smiles_input"] = ""
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
            results = run_full_pipeline(smiles_codes, email)

        st.session_state["results"] = results
        st.session_state["analysis_ready"] = True
        st.session_state["run_error"] = ""
        st.success("Analysis completed successfully!")
    except Exception as e:
        st.session_state["run_error"] = f"An error occurred during the analysis: {str(e)}"
        st.session_state["analysis_ready"] = False

left, right = st.columns([2,1])

with left:
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

with right:
    st.subheader("Options")
    st.button("Load example",on_click=load_example, width='stretch')
    st.button("Clear", on_click=clear_inputs, width='stretch')
    
    st.button("Add SMILES code", on_click=show_additional_input, width="stretch")
    if st.session_state["show_add_smiles"]:
        st.text_input(
            "Add one more SMILES code",
            key = "new_smiles_input",
            placeholder="Please enter a new SMILES code"
        )
        st.button("Store new SMILES", on_click=store_additional_smiles, width='stretch')

    st.button("Delete SMILES code", on_click=show_delete_smiles, width='stretch')
    if st.session_state["show_delete_smiles"]:
        submitted = st.session_state.get("submitted_smiles", [])

        if submitted:
            st.multiselect(
                "Select SMILES code to delete",
                options = submitted,
                key = "smiles_to_delete",
                placeholder="Choose one or more SMILES codes",
            )
            st.button("Confirm", on_click=update_smiles, width="stretch")
        
        else:
            st.warning("There are no saved SMILES codes to delete")

st.info(
    "Enter your SMILES, run the analysis and view the results below on this page!"
)

st.divider()

# Hacer que este botón NO modifique los SMILES adicionales
st.button("Save input", on_click=save_input, type="primary", width='stretch')

if "submitted_smiles" in st.session_state:
    st.subheader("Saved SMILES & email")
    st.write(f"Email: {st.session_state['submitted_email']}")
    st.write(f"Number of SMILES code stored: {len(st.session_state['submitted_smiles'])}")
    # Enseña una tab que clicas y te enseña los SMILES
    with st.expander("Show SMILES list"):
        for i, smiles in enumerate(st.session_state["submitted_smiles"],start=1):
            st.write(f"{i}. {smiles}")

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
    col2.metric("Interactions", len(results.get("df_interactions", [])))
    col3.metric("Pathways", len(results.get("df_pathways", [])))
    
    bp_count = len(results.get("df_go_bp", pd.DataFrame().get("go_id", pd.Series(dtype=str)).nunique()))
    mf_count = len(results.get("df_go_mf", pd.DataFrame().get("go_id", pd.Series(dtype=str)).nunique()))
    cc_count = len(results.get("df_go_cc", pd.DataFrame().get("go_id", pd.Series(dtype=str)).nunique()))

    st.markdown('#### GO enrichment')
    go_total, go1, go2, go3 = st.columns(4)
    go_total.metric("GO terms", bp_count + mf_count + cc_count)
    go1.metric("Biological process", bp_count)
    go2.metric("Molecular function", mf_count)
    go3.metric("Cellular component", cc_count)

    
    st.subheader("Analysis results")
    st.markdown("## Compounds")
    st.write(results["compound_names"])

    st.divider()

    st.markdown("## Interactions")
    st.dataframe(results["df_interactions"], width="stretch")

    st.divider()

    st.markdown("## Pathways")
    st.dataframe(results["df_pathways"], width="stretch")   

    st.divider()

    st.markdown("## GO enrichment")
    with st.expander("Biological process"):
        st.dataframe(results["df_go_bp"], width="stretch")
    with st.expander("Molecular function"):
        st.dataframe(results["df_go_mf"], width="stretch")
    with st.expander("Cellular component"):
        st.dataframe(results["df_go_cc"], width="stretch")

    st.divider()

    st.markdown("## Final summary")
    st.dataframe(results["final_summary"], width="stretch")

    st.divider()
    st.markdown("## Final summary + GO terms")
    st.dataframe(results["final_summaryGO"], width="stretch")

    """tabs = st.tabs([
        "Interactions",
        "Pathways",
        "GO biological process",
        "GO molecular function",
        "GO cellular component",
        "Final summary",
        "Final summary + GO"
    ])

    with tabs[0]:
        st.dataframe(results["df_interactions"], width="stretch")
    with tabs[1]:
        st.dataframe(results["df_pathways"], width="stretch")
    with tabs[2]:
        st.dataframe(results["df_go_bp"], width="stretch")
    with tabs[3]:
        st.dataframe(results["df_go_mf"], width="stretch")
    with tabs[4]:
        st.dataframe(results["df_go_cc"], width="stretch")
    with tabs[5]:
        st.dataframe(results["final_summary"], width="stretch")
    with tabs[6]:
        st.dataframe(results["final_summaryGO"], width="stretch")"""

#ADD DOWNLOADS TODO