2 scripts:
1) pdb_from_frags takes parameters from frags.100.5mers files (which are created by Rosetta fragment picker), 
that are neccessary for excising fragments from PDB structures, and write them to a new file for further proccessing. 
The new file (frags_params) will contain a pdb name, a chain, start and stop residues, and a fragment sequence
(for creating a resfile)
2) excise_frag_run_fixbb.py - excise fragments from pdb structures, create resfile for each fragment separately
and run fixbb design with this resfile.
