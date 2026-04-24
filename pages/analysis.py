# I dont know what to do with this page TODO
"""st.subheader("Options")
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

st.divider()"""