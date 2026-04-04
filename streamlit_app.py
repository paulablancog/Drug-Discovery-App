import streamlit as st

st.set_page_config(
    page_title="Drug Discovery Web",
    layout="wide"
)
 
home_page = st.Page("pages/home.py", title = "Home")
analysis_page = st.Page("pages/analysis.py", title="Analysis")
compounds_page = st.Page("pages/compound.py", title = "Compounds")
interactions_page = st.Page("pages/interactions.py", title = "Interactions")
pathways_page = st.Page("pages/pathways.py", title = "Pathways")
go_page = st.Page("pages/go_enrichment.py", title = "GO enrichment")
summary_page = st.Page("pages/final_summary.py", title = "Final summary")

pg = st.navigation(
    [
        home_page,
        analysis_page,
        compounds_page,
        interactions_page,
        pathways_page,
        go_page,
        summary_page
    ]
)

with st.sidebar:
    st.title("Drug Discovery")
    st.caption("SMILES -> proteins -> pathways -> GO")

    if "submitted_smiles" in st.session_state:
        st.write(f"SMILES loaded: {len(st.session_state['submitted_smiles'])}")
    if st.session_state.get("analysis_ready", False):
        st.success("Analysis completed! Check the results in the pages.")
    else:
        st.info("Run the analysis from Home")

pg.run()