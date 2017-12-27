2 scripts:
1) pdb_from_frags takes parameters from frags.100.5mers files (which are created by Rosetta fragment picker), that are neccessary for excising fragments from PDB structures, and write them to a new file for further proccessing. 
The new file (frags_params) will contain a pdb name, a chain, start and stop residues, and a fragment sequence (for creating a resfile)
2) excisepdb_run - excise fragments from pdb structures. This scripts uses several another scripts. It call prody to fetch a structure, then uses excise_pdb_working.pl script for fragment extraction, and then delete a structure and store only the fragment. If the pdb file is problematic (wrong numbering) - calls also pdb_renumber.py script.
