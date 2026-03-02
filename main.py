# Streamlit code
import pandas as pd
import app.chem
import app.interactions
import app.utils
import app.proteins
import app.pathways

# Asking for compound SMILES code
smiles_code = input("Enter the SMILES code: ")
compound = app.chem.retrieve_compound(smiles_code)

if compound:
    print("Compound found!")
    compound_info = app.chem.compound_information(compound)
    print("\n")
    print(compound_info)
else:
    print("No compound found, check the SMILES code entered again")
    raise SystemExit(1)

# 1. Creates the URL
index_url = f"{app.utils.URL_base}/rest/pug_view/index/compound/{compound.cid}/JSON"
# 2. Retrieves the index JSON
index_json = app.utils.get_json(index_url)
# 3. Saves the index JSON to a file
app.utils.save_json(index_json, f"compound_{compound.cid}_{compound.synonyms[0]}__index.json", "indexes")
# 4. Gets all sections in the index JSON
# 5. Finds the Interactions and Pathways section and retrieves the data of interactions
out, data = app.interactions.get_all_sections(compound)

if data is None:
    print("No interactions data for this compound")
    raise SystemExit(0)

# 6. Get the tables externally
tables = app.interactions.retrieve_externaltable(data)

# 7. Download and save each table
where_pathways = {"ands": [{"cid": str(compound.cid)}, {"core": "1"}]}
rows = None

for subsection, table_list in tables:
    for table_name in table_list:
        subsection_1 = (subsection or "").lower()
        # -- SPECIAL FOR PATHWAYS TABLES --
        if subsection_1 == "pathways":
            print("Special pathways table: ", table_name)
            rows = app.interactions.get_interactions_table(compound, "pathway", where = where_pathways, order="pathwayid,asc")   
            print(subsection, "pathways", "rows: ", len(rows))
        # -- NORMAL SDQ EXTERNAL TABLES -- 
        elif subsection_1 == "chemical-target interactions": 
            rows = app.interactions.get_interactions_table(compound, table_name,order="geneid,asc")
            print(subsection, table_name, "rows: ", len(rows))

        # -- SAVING INTERACTION TABLES IN SEPARATED SUBSECTIONS FOLDERS --
        json_path = app.utils.save_rows_json(rows, compound, f"1.{subsection}/compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.json")
        print("\nSuccessfully saved interactions table as JSON for compound: "+compound.synonyms[0])
        csv_path = app.utils.save_rows_csv(rows, compound, f"1.{subsection}/compound_{compound.cid}_{compound.synonyms[0]}_interactionstable.csv")
        print("\nSuccessfully saved interactions table as CSV for compound: "+compound.synonyms[0])

# -- RETRIEVE PROTEINS IN INTERACTIONS (PROTEIN COUNT) --
app.proteins.retrieve_targets_1(compound)

# -- RETRIEVE THE PATHWAY (PATHWAYS COUNT) --
app.pathways.retrieve_pathways(compound)

# -- RETRIEVE PATHWAY PROTEINS (PROTEIN COUNT INSIDE A PATHWAY) --
#app.pathways.retrieve_proteins_from_pathway(compound)

# Bethanechol: CC(C[N+](C)(C)C)OC(=O)N  
# Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C  
# Carbachol: C[N+](C)(C)CCOC(=O)N.[Cl-]  
# Ethanolamine: C(CO)N  
# Forskolin: CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O  
# Ginsenoside Rb1: CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C  
# Maprotiline: CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24  
# Pilocarpine: CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C  
