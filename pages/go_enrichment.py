import streamlit as st
import pandas as pd

st.title("GO enrichment")

if("results") not in st.session_state:
    st.warning("No analysis has been run yet. Run the analysis from Home page first.")
    st.stop()

results = st.session_state.get("results", None)
if results is None:
    st.warning("No analysis results found in the session. Run the analysis from Home page first.")
    st.stop()

df_bp = results.get("df_go_bp", pd.DataFrame())
df_mf = results.get("df_go_mf", pd.DataFrame())
df_cc = results.get("df_go_cc", pd.DataFrame())

# USE EXPANDERS?? TODO
st.subheader("Biological process GO terms")
st.dataframe(df_bp, width="stretch")

st.subheader("Molecular function GO terms")
st.dataframe(df_mf, width="stretch")

st.subheader("Cellular component GO terms")
st.dataframe(df_cc, width="stretch")