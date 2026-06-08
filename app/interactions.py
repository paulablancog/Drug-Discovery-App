import json

import app.utils


# PIPELINE to retrieve interactions and pathways of a compound given its cid

# 1) Download PUG View index of the compound (UTILS)
# 2) Checks if the compound has an "Interactions and Pathways" section
# 3) If it does, download that section as a JSON file
# 4) Scans JSON that points to proteins, genes and pathways (for later use PUG-REST information geneID, proteinID and pathway info + interactions)
# 5) Returns those IDs
 
def load_interactions_and_pathways_data(compound):
    """Loads the Interactions and Pathway section from compound and returns a JSON"""
    # 1. Creates the URL
    url = f"{app.utils.URL_BASE}/rest/pug_view/data/compound/{compound.cid}/JSON?heading=Interactions%20and%20Pathways"
    # 2. Retrieves the index JSON data of interactions and pathways
    return app.utils.get_json(url)


def flatten_record_sections(record_json):
    """Returns all Sections from a Record in JSON"""
    record = record_json.get("Record", {})
    sections = record.get("Section", []) or []

    out = []

    for section in sections:
        out.extend(flatten_section(section))
    return out  

def flatten_section(section):
    """Returns each subsection from a bigger Section"""
    out = [section]
    for subsection in section.get("Section", []) or []:
        out.extend(flatten_section(subsection))
    return out # returns a list of subsections inside a single Section

def find_section_by_heading(record_json, heading):
    """Looks for the first section whose TOCHeading matches heading"""
    heading = heading.strip().lower()

    for section in flatten_record_sections(record_json):
        if section.get("TOCHeading", "").strip().lower() == heading:
            return section
    
    return None

def get_all_sections(index_json):
    """Flatten and return all sections from the initial index JSON"""
    return flatten_record_sections(index_json)

def has_interactions_and_pathways(index_json):
    """Checks if the compound has Interactions and Pathways section"""
    return find_section_by_heading(index_json, "Interactions and Pathways") is not None

def retrieve_externaltable(interactions_json):
    """Return the subsection name (Chemical-Target or Pathways) and its external name table"""
    section = find_section_by_heading(interactions_json, "Interactions and Pathways")
    if section is None:
        return []
    
    results = []

    # Get subsection (Chemical-Target Interactions or Pathways) and the ExternalTableName of each Information in each subsection (if it exists)
    for subheader in section.get("Section", []):
        subsection = subheader.get("TOCHeading", "")
        externaltables = [] # external table in each subsection (Chemical-Target or Pathways)
        
        for information in subheader.get("Information", []):
            externaltablename = information.get("Value", {}).get("ExternalTableName")
            if externaltablename:
                externaltables.append(externaltablename) 

        tables = list(dict.fromkeys(externaltables)) # remove duplicates
        results.append((subsection, tables)) # add the subsection and its tables to the results list
    
    return results 


# Build a SDQ query for PubChem for Chemical-Target Interactions External Tables
def sdq_query_externaltable(collection, where, select = "*", start = 1, limit = 1000, order = "cid,asc"):
    """Run a SDQ query on PubChem and return the resulting JSON"""
   
    if order is None:
        order_list = []
    elif isinstance(order, str):
        order_list = [order]
    else: 
        order_list = list(order)
    
    query = {
        "select": select, # which columns you want back
        "collection": collection, # which table to query
        "order": order_list, # sort order ("cid, asc" means sort by cid ascending)
        "start": start, # pagination start row (1=start row)
        "limit": limit, # how many rows to return 
        "where": where, # the filter condition (example: cid must equal 5831)
        "width": 10000000  # allow wide fields (for long text citations for example)
    }

    # endpoint for SDQ queries
    url = f"{app.utils.URL_BASE}/sdq/sphinxql.cgi"
    params = {
        "infmt": "json",
        "outfmt": "json",
        "query": json.dumps(query)
    }

    return app.utils.get_json(url, params=params)

def normalize_taxonomy_ids(selected_tax_ids):
    """Normalizes the selected taxonomy IDs by the user."""
    
    if selected_tax_ids is None:
        return set()
    if isinstance(selected_tax_ids, str):
        selected_tax_ids = [selected_tax_ids]
    
    return {
        str(taxid).strip() for taxid in selected_tax_ids if str(taxid).strip()
    }

def match_taxonomy(row, selected_tax_ids=None):
    """Checks if the taxonomy ID of the row matched any of the selected taxonomies by the user."""
    
    selected_tax_ids = normalize_taxonomy_ids(selected_tax_ids)

    if not selected_tax_ids:
        return True

    taxid = str(row.get("taxid", "")).strip()
    return taxid in selected_tax_ids

def filter_taxonomy(rows, selected_tax_ids=None):
    """Filter rows by taxonomy ID, if selected taxonomy is provided, if not returns all rows."""
    
    return [row for row in rows
            if isinstance(row,dict) and match_taxonomy(row, selected_tax_ids)]

def get_interactions_table(compound, collection, page_size = 1000, where = None, order="cid,asc", selected_tax_ids=None): 
    """Retrieve all rows from PubChem SDQ External table with pagination"""
    
    if where is None:
        where = {"ands": [{"cid": str(compound.cid)}]} # default where clause if not provided
    
    request = sdq_query_externaltable(collection, where, start = 1, limit = min(page_size,1000), order=order)
    
    if request is None:
        return[]
    
    out_set = request.get("SDQOutputSet", [])
    if not out_set:
        print("No data found in the SDQ query.")
        return []
    
    block = out_set[0]

    status = block.get("status", {}) or {}
    status_code = status.get("code", 0)

    if str(status_code) != "0":
        print(f"SDQ query returned non-zero status code: {status_code}, failed for collection {collection}: {status}")
        return []
    
    rows = block.get("rows", []) or []
    if not isinstance(rows, list):
        rows = []

    try:    
        total = int(block.get("totalCount", len(rows)))
    except (TypeError, ValueError):
        total=len(rows)
    
    all_rows = list(rows)
    start = 1+len(rows)

    while len(all_rows) < total:
        data = sdq_query_externaltable(collection, where, start = start, limit = min(page_size,1000), order=order)
        
        if data is None:
            break
        
        rows_page = data.get("SDQOutputSet", [{}])[0].get("rows", []) or []
        if not rows_page:
            break

        all_rows.extend(rows_page)
        start += len(rows_page)
    
    if selected_tax_ids is not None:
        all_rows = filter_taxonomy(all_rows, selected_tax_ids)

    return all_rows
