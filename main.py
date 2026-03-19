# Streamlit code
import app.pipeline
import pandas as pd

interaction_compounds = [
    "bethanechol",
    "caffeine",
    "ethanolamine",
    "carbachol",
    "forskolin",
    "ginsenoside rb1",
    "maprotiline",
    "pilocarpine"
]

smiles_codes = [
    "CC(C[N+](C)(C)C)OC(=O)N",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "C[N+](C)(C)CCOC(=O)N.[Cl-]"
    "CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O",
    "CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C",
    "CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24",
    "CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C",
]

pathway_compounds = [
    "CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "C(CO)N",
]

email = "paulablglez@gmail.com"

for smiles in smiles_codes:
    smiles_code = smiles
    compound, compound_info, df_proteins = app.pipeline.fetch_pubchem_compound(smiles_code, email)
    print(compound_info)

    df_interactions = app.pipeline.fetch_interactions_summary(interaction_compounds)
    df_pathways = app.pipeline.fetch_pathway_summary(pathway_compounds)
    final_summary = app.pipeline.build_final_summary(df_interactions, df_pathways)
    df_go, final_summaryGO = app.pipeline.build_go_enrichment(final_summary)

#smiles_code = input("Enter the SMILES code: ")

print(final_summaryGO.head(20))


final_summaryGO.to_csv("Protein_final_summaryGO.csv", index=False)
excelsummary = pd.read_csv("Protein_final_summaryGO.csv")
excelsummary.to_excel("Protein_final_summaryGO.xlsx", index=False)


# Bethanechol: CC(C[N+](C)(C)C)OC(=O)N  
# Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C  
# Carbachol: C[N+](C)(C)CCOC(=O)N.[Cl-]  
# Ethanolamine: C(CO)N  
# Forskolin: CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O  
# Ginsenoside Rb1: CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C  
# Maprotiline: CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24  
# Pilocarpine: CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C  
