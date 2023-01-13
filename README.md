# Lymphoma KB

### Setting up knowledge base
1. Install OrientDB and run from Docker Desktop
3. Run `graph.py` to initialize and populate graph knowledge base with variants and statements
4. To access the database, go to http://localhost:2480/ and login to `lymphomakb` database with username `root` and password `root`

Currently incorporated sources:
- NCIt: https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip
- ClinVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Currently supported variant types:
- Single nucleotide variant
- Copy gain/loss
- Indel
- Insertion
- Deletion

No. of diseases: 889

No. of variants: 201

No. of statements: 399