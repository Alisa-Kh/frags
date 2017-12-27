import subprocess
import os
import sys

# Get filename as an argument to the script
parameters_file_path = sys.argv[1]

"""Update PYTHONPATH to be /vol/ek/share/rosetta/rosetta_src_2017.45.59812_bundle/tools/protein_tools/"""

# Open the frags_parameters, extract and append parameters to different lists
f = open(parameters_file_path, 'r')
list_of_lines = f.readlines()
pdbs = []
chains = []
start_res = []
end_res = []
for line in list_of_lines:
    pdbs.append(line.split()[0])
    chains.append(line.split()[1])
    start_res.append(line.split()[2])
    end_res.append(line.split()[3])

# Create a list with output names for excisePdb
# output_names = []
# for i in range(len(pdbs)):
#     output_names.append(pdbs[i]+'.'+chains[i]+'.'+start_res[i]+'.'+end_res[i]+'.pdb')

# Fetch PDBs with prody, call excisePdb and then delete the pdb file
for pdb, chain, start, end in zip(pdbs, chains, start_res, end_res):
    outfile = pdb +'.'+chain +'.'+start +'.'+end+'.pdb'
    pdb_full = pdb + '.pdb'

    subprocess.run(['/vol/ek/share/csbw/tools/env2017/bin/prody', 'fetch', pdb])

    subprocess.run(['/vol/ek/Home/alisa/scripts/check_my_scripts/excisePdb_working.pl',
                        pdb_full, chain, start, end, outfile])

    if not os.path.exists(outfile):
        print("renumbering:")
        new_pdb = pdb_full
        renumbering = \
            '/vol/ek/share/rosetta/rosetta_src_2017.45.59812_bundle/tools/protein_tools/scripts/pdb_renumber.py %s %s' \
                      % (pdb_full, new_pdb)
        os.system(renumbering) # Renumber PDBs which don't start from 1

        subprocess.run(['/vol/ek/Home/alisa/scripts/check_my_scripts/excisePdb_working.pl',
                        new_pdb, chain, start, end, outfile])

    os.remove('/vol/ek/Home/alisa/scripts/check_my_scripts/' + pdb_full)

