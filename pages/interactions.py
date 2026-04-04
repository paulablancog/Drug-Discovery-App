import streamlit as st
import pandas as pd

st.title("Interactions")

if("results") not in st.session_state:
    st.warning("No analysis has been run yet. Run the analysis from Home page first.")
    st.stop()

results = st.session_state.get("results", None)

if results is None:
    st.warning("No analysis results found in the session. Run the analysis from Home page first.")
    st.stop()

df_interactions = results.get("df_interactions", pd.DataFrame())

st.write(f"Interaction proteins identified: {len(df_interactions)}")
st.dataframe(df_interactions, width="stretch")