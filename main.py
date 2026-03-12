# Streamlit code
import pandas as pd
import app.chem
import app.interactions
import app.utils
import app.proteins
import app.pathways
from pathlib import Path

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
# 1) Write protein_data.txt with geneids
app.proteins.retrieve_targets_1(compound)

# 2) NCBI: geneid -> symbol + description
email = "paulablglez@gmail.com"
out = app.proteins.translate_geneid_to_protein(email,"protein_data.txt", compound)

# 3) UniProt: geneid -> uniprot accession
df_map = app.proteins.map_genes_to_uniprot("protein_data.txt", compound)

# 4) Merge both data Frames and save
out_nodups = out.drop_duplicates(subset=["geneid"])
out2 = out_nodups.merge(df_map, on="geneid", how="left")
out2.to_csv(f"{compound.synonyms[0]}_geneid_symbols_uniprot.csv", index=False)

# -- RETRIEVE protein_data OF EACH COMPOUND IN A DATAFRAME (Chemical-Target Interactions)
compounds = ["bethanechol", "caffeine", "Ethanolamine", "carbachol", "forskolin", "Ginsenoside rb1", "maprotiline", "pilocarpine"]
csv_files = [f"{compound}_geneid_symbols_uniprot.csv" for compound in compounds]

dfs_int = []
for csv in csv_files:
    if not Path(csv).exists():
        print("Not a csv for Protein Interactions: "+csv)
        continue
    df = pd.read_csv(csv)
    dfs_int.append(df)

df_interactions = pd.concat(dfs_int, ignore_index=True) if dfs_int else pd.DataFrame(columns = ["uniprot_accession", "compound", "symbol", "geneid"])
df_interactions.to_csv("ProteinInteractions.csv")

#summary = app.proteins.results_summary_count(csv_files)
#summary.to_csv("Protein_Mapping_Interactions.csv", index = False)
#print(summary.head(20))

# --RETRIEVE ALL PROTEINS FROM PATHWAY COMPOUNDS (also creates pathways.txt)
compounds = ["CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C(CO)N"]
dfs_proteins = []

# TODO: VER SI SE DUPLICAN PATHWAYS
for c in compounds:
    compound = app.chem.retrieve_compound(c)
    filename = f"{compound.synonyms[0]}_ProteinsAllPathways.csv"
    compound_file = Path(filename)
    if compound_file.exists():
        print("Protein Pathways already exists")
        df = pd.read_csv(compound_file)
        dfs_proteins.append(df)
        continue

    df = app.pathways.retrieve_pathways(compound)
    if df is None or df.empty:
        continue

    df.to_csv(filename, index=False)
    dfs_proteins.append(df)

df_all = pd.concat(dfs_proteins, ignore_index=True) if dfs_proteins else pd.DataFrame(columns=["uniprot_accession", "protein_name", "pathway", "compound"])
df_all.to_csv("AllProteinsPathways.csv", index=False)

# --READ BOTH PROTEIN INTERACTIONS (chemical & pathways)

filename = "Protein_mapping.csv"
proteinmapping_file = Path(filename)
if proteinmapping_file.exists():
    print("Protein_mapping.csv already exists, loading it...")
    final_summary = pd.read_csv(proteinmapping_file)
else:
    df_interactions = pd.read_csv("ProteinInteractions.csv")
    df_pathways = pd.read_csv("AllProteinsPathways.csv")

    df_interactions.columns = df_interactions.columns.str.strip()
    df_pathways.columns = df_pathways.columns.str.strip()

    df_interactions["uniprot_accession"] = (df_interactions["uniprot_accession"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NaN": ""})
    )
    df_interactions = df_interactions[df_interactions["uniprot_accession"] != ""].copy()

    df_pathways["uniprot_accession"] = (df_pathways["uniprot_accession"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NaN": ""})
    )
    df_pathways = df_pathways[df_pathways["uniprot_accession"] != ""].copy()


    interactions_summary = (
        df_interactions.groupby("uniprot_accession", as_index=False).agg(
            total_count = ("uniprot_accession", "size"),
            n_compounds = ("compound", "nunique"),
            compounds=("compound", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            symbol=("symbol", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
            geneid=("geneid", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        )
    )

    pathway_summary = (
        df_pathways.groupby("uniprot_accession", as_index=False).agg(
            n_pathways = ("pathway", "nunique"),
            pathways = ("pathway", lambda x: ";".join(sorted(set(x)))),
        )
    )

    final_summary = (
        interactions_summary.merge(pathway_summary, on="uniprot_accession", how="outer").fillna(
            {"n_pathways":0, "pathways": ""}).sort_values(["total_count", "n_compounds"], ascending=[False, False])
        )

    final_summary.to_csv("Protein_mapping.csv", index=False)

print(final_summary.head(20))

df_go = app.proteins.fetch_goterms("Protein_mapping.csv", aspects=["biological_process", "molecular_function", "cellular_component"])
print(df_go.head(20))
print("GO Terms")
df_go.to_csv("Go_terms.csv", index=False)

go_name = app.proteins.fetch_gonames(df_go["go_id"].dropna().unique())
df_go = df_go.merge(go_name, on="go_id", how="left")
df_go.to_csv("Go_terms_Names.csv", index = False)

if df_go.empty:
    final_summaryGO = final_summary.copy()
    final_summaryGO["n_go_terms"] = 0
    final_summaryGO["go_ids"] = ""
    final_summaryGO["go_names"] = ""
    final_summaryGO["evidence_code"] = ""

    final_summaryGO["go_bp_ids"] = ""
    final_summaryGO["go_bp_names"] = ""
    
    final_summaryGO["go_mf_ids"] = ""
    final_summaryGO["go_mf_names"] = ""
    
    final_summaryGO["go_cc_ids"] = ""
    final_summaryGO["go_cc_names"] = ""

else:
    df_go = df_go.drop_duplicates(subset=["uniprot_accession", "go_id", "aspect", "evidence_code"]).copy()

    go_summary = (
        df_go.groupby("uniprot_accession", as_index=False).agg(
        n_go_terms = ("go_id", lambda x: x.dropna().nunique()),
        go_ids = ("go_id", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        go_names = ("go_name", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),
        evidence_code = ("evidence_code", lambda x: ";".join(sorted(set(str(v).strip() for v in x if pd.notna(v) and str(v).strip())))),  
        )
    )

    bp_summary = app.proteins.summarize_goaspect(df_go, "biological_process", "bp")
    mf_summary = app.proteins.summarize_goaspect(df_go, "molecular_function", "mf")
    cc_summary = app.proteins.summarize_goaspect(df_go, "cellular_component", "cc")


    final_summaryGO = final_summary.merge(go_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(bp_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(mf_summary, on="uniprot_accession", how="left")
    final_summaryGO = final_summaryGO.merge(cc_summary, on="uniprot_accession", how="left")
    
    final_summaryGO["n_go_terms"] = final_summaryGO["n_go_terms"].fillna(0).astype(int)

    for col in ["go_ids", "go_names", "evidence_code", 
                "go_bp_ids", "go_bp_names", 
                "go_mf_ids", "go_mf_names", 
                "go_cc_ids", "go_cc_names"]:
        final_summaryGO[col] = final_summaryGO[col].fillna("")

final_summaryGO.to_csv("Protein_mappingGO.csv", index=False)
excelsummary = pd.read_csv("Protein_mappingGO.csv")
excelsummary.to_excel("Protein_mapping.xlsx", index=False)


#print(final_summaryGO.head(20))


# Bethanechol: CC(C[N+](C)(C)C)OC(=O)N  
# Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C  
# Carbachol: C[N+](C)(C)CCOC(=O)N.[Cl-]  
# Ethanolamine: C(CO)N  
# Forskolin: CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O  
# Ginsenoside Rb1: CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C  
# Maprotiline: CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24  
# Pilocarpine: CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C  
