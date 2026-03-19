import pubchempy as pcp

# Step 1: Recognize the drug by SMILES code with PubChem database

def compound_retrieval(smiles_code):
    """Retrieve PubChem compound from its SMILES code"""
    smiles_code = str(smiles_code).strip()
    if not smiles_code:
       print("No SMILES code entered")
       return None 
    
    compounds = pcp.get_compounds(smiles_code, namespace = "smiles") 
    if not compounds: 
        return None
    
    return compounds[0]

def compound_information(compound):
    """Return compound metadata dictionary from a PubChem compound"""
    synonyms = getattr(compound, "synonyms", None) or []
    first_synonym = synonyms[0] if synonyms else ""

    return {
        'cid': getattr(compound, "cid", None),
        'synonyms': str(first_synonym).strip(),
        'name': getattr(compound, "iupac_name", None),
        'molecular_formula': getattr(compound, "molecular_formula", None),
        'molecular_weight': getattr(compound, "molecular_weight", None),
    }


def compound_display_name(compound):
    """Return a clear compound name"""
    synonyms = getattr(compound, "synonyms", None) or []
    if synonyms:
        return str(synonyms[0]).strip()
    
    iupac_name = getattr(compound, "iupac_name", None)
    if iupac_name:
        return str(iupac_name).strip()
    
    cid = getattr(compound, "cid", None)
    if cid is not None:
        return f"compound_{cid}"

    return "compound"
