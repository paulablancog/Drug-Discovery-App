import streamlit as st
import pandas as pd

st.title("Pathways")

if("results") not in st.session_state:
    st.warning("No analysis has been run yet. Run the analysis from Home page first.")
    st.stop()

results = st.session_state.get("results", None)
if results is None:
    st.warning("No analysis results found in the session. Run the analysis from Home page first.")
    st.stop()

df_pathways = results.get("df_pathways", pd.DataFrame())

st.write(f"Pathways identified: {len(df_pathways)}")
st.dataframe(df_pathways, width="stretch")