import streamlit as st
import pandas as pd
import app.chem

st.title("Compounds")

# For this page, as soon as the user enters the SMILES code, it should be available for the user
# to see information about the compound before the full downstream analysis run

# 1. Compound identification (DISPLAY) from submitted_smiles in home.py??

# Pensamiento: que sirva como una page para que el usuario pueda elegir libremente qué compuestos elegir
# para el analisis -> Opción para ingresar un SMILES code (meter PubChem link por ejemplo)

# Compounds debe leer de los submitted_smiles y compound_results (de analisis page) y no necesita
# resultados finales del pipeline

if "submitted_smiles" not in st.session_state:
    st.warning("No SMILES code input in the session. Go to Home Page first")
    st.stop()

submitted_smiles = st.session_state["submitted_smiles"]

def identify_compound(smiles_list):
    rows = []
    for smiles in smiles_list:
        try:
            print(smiles)
            compound = app.chem.compound_retrieval(smiles)
            if compound is None:
                rows.append(
                    {
                        "input_smiles": smiles,
                        "cid": None,
                        "compound_name": None,
                        "molecular_formula": None,
                        "molecular_weight": None,
                        "status": "No compound found",
                    }
                )
                continue
        
            compound_info = app.chem.compound_information(compound)
            compound_name = app.chem.compound_display_name(compound)

            rows.append(
                    {
                        "input_smiles": smiles,
                        "cid": compound_info.get("cid"),
                        "compound_name": compound_name,
                        "iupac_name": compound_info.get("name"),
                        "molecular_formula": compound_info.get("molecular_formula"),
                        "molecular_weight": compound_info.get("molecular_weight"),
                        "status": "Identified",
                    }
                )
        except Exception as ex:
            rows.append(
                {
                    "input_smiles": smiles,
                    "cid": None,
                    "compound_name": None,
                    "molecular_formula": None,
                    "molecular_weight": None,
                    "status": f"error: {ex}"
                }
            )
    return pd.DataFrame(rows)

if "compound_results" not in st.session_state or st.session_state.get("compound_results_source") != tuple(submitted_smiles):
    # Shows a temporary working... message while the code runs
    with st.spinner("Identifying compounds from submitted SMILEs"):
        compound = identify_compound(submitted_smiles)
        st.session_state["compound_results"] = compound
        # So that it refreshes when the submitted_smiles change, but does not change when the user navigates through pages
        st.session_state["compound_results_source"] = tuple(submitted_smiles)
df_compounds = st.session_state["compound_results"]

st.write(f"Compounds identified from {len(submitted_smiles)} submitted SMILES:")
st.dataframe(df_compounds, width='stretch')
