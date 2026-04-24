import streamlit as st
import pandas as pd

st.title("Final summary")

if("results") not in st.session_state:
    st.warning("No analysis has been run yet. Run the analysis from Home page first.")
    st.stop()

results = st.session_state.get("results", None)
if results is None:
    st.warning("No analysis results found in the session. Run the analysis from Home page first.")
    st.stop()

df_summary = results.get("final_summary", pd.DataFrame())
df_summary_go = results.get("final_summaryGO", pd.DataFrame())

st.subheader("Final summary of compounds, proteins and pathways")
st.dataframe(df_summary, width="stretch")

# Include in download section that GO terms are included in the final summary