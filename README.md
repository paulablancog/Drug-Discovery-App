## Drug Discovery App

This project consists of a web-based application for **compound-centered drug discovery analysis.** The application is designed to start from user-provided **SMILES codes**, identify the corresponding compounds and retrieve associated **compound-protein interactions** and **pathway-related proteins**, and performs **Gene Ontology (GO) enrichment** by accessing external databases such as **PubChem**, **NCBI** and **QuickGO.**
This application enables users to carry out interactive exploratory analysis in a web environment.

## Key Features

### Compound identification
- PubChem compound identification via SMILES code
- Retrieval of basic compound information 

### Compound-Protein interaction
- Direct compound-protein interactions retrieval from PubChem SDQ external tables

### Pathway-related proteins
- Retrieval of compound-associated pathways
- Retrieval of pathway-associated proteins from PubChem pathway endpoints

### GO term annotation
- Annotation of protein-associated GO terms using QuickGO
- Grouping of GO terms by aspect:
    - Biological process
    - Molecular function
    - Cellular component
    
### GeneID translation & mapping
- Translation of protein-associated GeneIDs using NCBI Entrez
- Mapping of protein-associated GeneIDs to UniProt accession codes using UniProt ID service

### Additional features
- Interactive Streamlit interface
- Excel export of analysis results
- No intermediate in disk file writing


## Workflow Overview
The main workflow of the application is:

1. **SMILES** -> PubChem compound retrieval
2. Compound metadata extraction from PubChem
3. Interactions and pathways retrieval
4. GeneID extraction, translation through NCBI Entrez, and UniProt mapping
5. Pathway protein retrieval
6. GO annotation retrieval and aspect-like grouping
7. Final protein summary construction
8. Excel export of the analysis results

## Project Structure (Simplified Overview)
```markdown
app/
├── chem.py
├── interactions.py
├── pathways.py
├── proteins.py
├── pipeline.py
├── utils.py
pages/
└── home.py
streamlit_app.py
```

## Getting started
1. Clone this repository
```bash
git clone https://github.com/paulablancog/TFG-.git
cd TFG-
```

2. Install dependencies
```bash
pip install -r requirements.txt
```

3. Run Streamlit application
```bash
streamlit run streamlit_app.py
```
Make sure you have Internet access before running the application, since queries to external databases and APIs are present.

## Author
This project has been developed as the **Final Thesis project at CEU San Pablo University** by:
@paulablancog