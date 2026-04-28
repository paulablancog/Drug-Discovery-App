# Streamlit code
import app.pipeline


smiles_codes = [
    "CC(C[N+](C)(C)C)OC(=O)N",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "C[N+](C)(C)CCOC(=O)N.[Cl-]",
    "C(CO)N",
    "CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O",
    "CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C",
    "CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24",
    "CC[C@H]1C@HCC2=CN=CN2C",
]

email = "paulablglez@gmail.com"

result = app.pipeline.run_full_pipeline(smiles_codes, email)
print(result["final_summaryGO"].head(20))

# Bethanechol: CC(C[N+](C)(C)C)OC(=O)N  
# Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C  
# Carbachol: C[N+](C)(C)CCOC(=O)N.[Cl-]  
# Ethanolamine: C(CO)N  
# Forskolin: CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O  
# Ginsenoside Rb1: CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C  
# Maprotiline: CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24  
# Pilocarpine: CC[C@H]1C@HCC2=CN=CN2C

#CC(C[N+](C)(C)C)OC(=O)N  
#CN1C=NC2=C1C(=O)N(C(=O)N2C)C  
#C[N+](C)(C)CCOC(=O)N.[Cl-]  
#C(CO)N  
#CC(=O)O[C@H]1[C@H]([C@@H]2[C@]([C@H](CCC2(C)C)O)([C@@]3([C@@]1(O[C@@](CC3=O)(C)C=C)C)O)C)O  
#CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1[C@@H](C[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O)C)C)O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O)O)O)C  
#CNCCCC12CCC(C3=CC=CC=C31)C4=CC=CC=C24  
#CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C