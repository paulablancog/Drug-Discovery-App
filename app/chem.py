import pubchempy as pcp
import pandas as pd
import requests
from rdkit import Chem

# Step 1: Recognize the drug by SMILES code with PubChem database
# Step 2: Retrieve the compound information (CID, molecular formula, molecular weight, synonyms, etc.) by candidate ranking list.

def validate_smiles(smiles_code):
    """Validates a SMILES code and returns its canonical form if valid, or an error if invalid."""
    
    smiles_code = str(smiles_code).strip()

    if not smiles_code:
        return None, "Empty SMILES"

    mol = Chem.MolFromSmiles(smiles_code)

    if mol is None:
        return None, f"Invalid SMILES: {smiles_code}"
    
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return canonical_smiles, None

def check_smiles_box(text, existing_smiles=None):
    """Rejects multiline SMILES input:
        -   Rejects invalid SMILES
        -   Converts valid SMILES into canonical SMILES
        -   Removes duplicates, including equivalent SMILES written differently"""
    
    existing_smiles = existing_smiles or []
    
    existing_canonical = set() 
    
    for smiles in existing_smiles:
        canonical, error = validate_smiles(smiles)
        if canonical and not error:
            existing_canonical.add(canonical)

    valid_smiles = []
    invalid_smiles = []
    duplicate_smiles = []

    seen = set()

    for line in str(text).splitlines():
        raw_smiles = line.strip()
        if not raw_smiles:
            continue

        canonical, error = validate_smiles(raw_smiles)
        
        if error:
            invalid_smiles.append(raw_smiles)
            continue
        
        if canonical in existing_canonical or canonical in seen:
            duplicate_smiles.append(raw_smiles)
            continue
        seen.add(canonical)
        valid_smiles.append(canonical)

    return{
        "valid": valid_smiles,
        "invalid": invalid_smiles,
        "duplicates": duplicate_smiles,
    }


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
        
        response = requests.post(url, data=data, timeout=60)
        if response.status_code == 404:
            return []
        
        response.raise_for_status()

        text = response.text.strip()
        if not text:
            return []
        
        cids = []
        for line in text.splitlines():
            line = line.strip()
            if line.isdigit():
                cids.append(int(line))

        return cids
        
    def load_compounds(cids):
        """Given a candidate list of CIDs, retrieve the corresponding PubChem compounds by stereochemical identity or by connectivity, and return the best candidate"""
        
        candidates = []

        for cid in cids:
            try:
                compound = pcp.Compound.from_cid(cid)
                candidates.append(compound)
            except Exception as ex:
                print(f"Error retrieving compound for CID {cid}: {ex}")
        
        return candidates

    def score(compound):
        """Score a compound based on the presence of metadata, to select the best candidate among the retrieved compounds"""
        
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
        print(f"No PubChem retrieval for SMILES {smiles_code}: {ex}")
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

def identify_compounds(smiles_list):
    """For a list of SMILES code, show the identified compound name"""
    rows = []

    for smiles in smiles_list:
        smiles = str(smiles).strip()
        if not smiles:
            continue
        
        canonical_smiles, smiles_error = validate_smiles(smiles)

        if smiles_error:
            rows.append({
                "smiles": smiles,
                "compound_name": "Invalid SMILES",
                "cid": None,
                "molecular_formula": None,
                "molecular_weight": None,
                "status": smiles_error,
            })
            continue

        try:
            compound = compound_retrieval(smiles)
            if compound is None:
                rows.append({
                    "smiles": smiles, 
                    "compound_name": "Not identified", 
                    "cid": None,
                    "molecular_formula": None,
                    "molecular_weight": None,
                    "status": "Valid SMILES, not identified in PubChem.",
                })
                continue

            compound_name = compound_display_name(compound)
            compound_info = compound_information(compound)

            rows.append({
                "smiles": smiles, 
                "compound_name": compound_name, 
                "cid": compound_info.get("cid"),
                "molecular_formula": compound_info.get("molecular_formula"),
                "molecular_weight": compound_info.get("molecular_weight"),
                "status": "Identified",
                })
            
        except Exception as ex:
            rows.append({
                    "smiles": smiles, 
                    "compound_name": "Not identified", 
                    "cid": None,
                    "molecular_formula": None,
                    "molecular_weight": None,
                    "status": f"Error: {ex}",
                })
            
    return pd.DataFrame(rows, columns=[
        "smiles", "compound_name", "cid", "molecular_formula", "molecular_weight", "status",
    ])
