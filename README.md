# bioinformatics_group_project_1
Introduction to Computational Biology - CS 4233 Bioinformatics - CS 5263 Fall 2019 Group Project 1

scrape_pdb.py first fetches a list of PDBs to scrape for from cullpdb_pc30_res3.0_R1.0_d191017_chains18877.gz
and tries to scrape the atom, helix, and sheet data in the following dataframes in pdb_parser class respectively: df_atom, df_helix, df_sheet


parser = pdb_parser() # To initialize parser. Empty argument will try to download all protein data in the gz file onto pdb_data directory and parse them into dataframes...

parser = pdb_parser(3) # To only parse first 3 proteins in the list for testing...
