import requests

URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data"

interaction_json = "Users\pblan\TFG_code\TFG-\1.Chemical-Target Interactions\compound_5831_carbachol_interactionstable.json"

# Appear in main so that the user can fetch the proteins in the interactions of the compound?
def retrieve_protein(compound,interaction_json, protacxn, foldername):
    out = []
    section = interaction_json.get([])
    for sections in section:
        if sections == "protacxn" and protacxn in sections:
            print(sections)
            out.extend(sections)

    print("Successfully retrieved Proteins IDs")
    filename = f"{compound.synonyms[0]}/protein_{protacxn}"

    return filename

filename = retrieve_protein("carbachol", interaction_json, "P20309")
print(filename)