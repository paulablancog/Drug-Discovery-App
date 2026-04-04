import pubchempy as pcp
import requests

# Step 1: Recognize the drug by SMILES code with PubChem database

def compound_retrieval(smiles_code):
    """Retrieve PubChem compound from its SMILES code"""
    smiles_code = str(smiles_code).strip()
    if not smiles_code:
       print("No SMILES code entered")
       return None 
    
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
    
    def fetch_cids(identity_type=None):
        data = {"smiles": smiles_code}
        if identity_type:
            data["identity_type"] = identity_type
        
        response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()
        text = response.text.strip()
        if not text:
            return None
        return [int(line.strip()) for line in text.splitlines() if line.strip()]
        
    def load_compounds(cids):
        candidates = []

        for cid in cids:
            try:
                compound = pcp.Compound.from_cid(cid)
                candidates.append(compound)
            except Exception as ex:
                print(f"Error retrieving compound for CID {cid}: {ex}")
        return candidates

    def score(compound):
        synonyms = getattr(compound, "synonyms", []) or []
        cid = getattr(compound, "cid", None)
        iupac_name = getattr(compound, "iupac_name", None)

        return(
            int(cid is not None),
            int(bool(iupac_name)),
            len(synonyms),
            -(cid if cid is not None else 10**18),
        )
    try:
        #1. try exact stereochemical identity
        cids = fetch_cids(identity_type="same_stereo")
        candidates = load_compounds(cids)
        if candidates:
            return max(candidates, key=score)
        
        #2. try exact connectivity
        cids = fetch_cids(identity_type="same_connectivity")
        candidates = load_compounds(cids)
        if candidates:
            return max(candidates, key=score)
        
        return None
    
    except Exception as ex:
        print(f"No PubChem retrieval: {ex}")
        return None
    

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
    
    for syn in synonyms[:10]:
        syn = str(syn).strip()
        if syn and len(syn)<60 and not syn.startswith("("):
            return syn
        
    iupac_name = getattr(compound, "iupac_name", None)
    if iupac_name:
        return str(iupac_name).strip()
    
    cid = getattr(compound, "cid", None)
    if cid is not None:
        return f"compound_{cid}"

    return "compound"


def compound_retrieval_anotherone(smiles_code):
    """Retrieve PubChem compound from its SMILES code"""
    smiles_code = str(smiles_code).strip()
    if not smiles_code:
       print("No SMILES code entered")
       return None 
    
    input_norm = "".join(smiles_code.split())
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
    try:
        response = requests.post(url, data={"smiles" : smiles_code}, timeout=30)
        response.raise_for_status()

        text = response.text.strip()
        if not text:
            return None
        cids = [int(line.strip()) for line in text.splitlines() if line.strip()]
        if not cids:
            return None
        
        candidates = []

        for cid in cids:
            try:
                compound = pcp.Compound.from_cid(cid)
                candidates.append(compound)
            except Exception as ex:
                print(f"Error retrieving compound for CID {cid}: {ex}")
        
        if not candidates:
            return None
        
        def score(compound):
            compound_smiles = "".join(str(getattr(compound, "smiles", "") or "").split())
            connectivity_smiles = "".join(str(getattr(compound, "connectivity_smiles", "") or "").split())

            synonym_count = len(getattr(compound, "synonyms", []) or [])
            cid = int(getattr(compound, "cid", 11**12) or 11**12) # The large number will serve for compounds missing a CID, just use a very large number

            return(
                int(input_norm == compound_smiles and compound_smiles != ""),
                int(input_norm == connectivity_smiles and connectivity_smiles != ""),
                synonym_count,
                -cid,
            )
        
        return max(candidates, key=score)
    
    except Exception as ex:
        print(f"No PubChem retrieval: {ex}")
        return None