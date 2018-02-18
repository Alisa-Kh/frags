import subprocess
import os
import sys

parameters_file_path = sys.argv[1] # Get filename as an argument to the script
original_seq = sys.argv[2] # peptide_sequence file

rosetta_database = '/vol/ek/share/rosetta/rosetta_src_2017.45.59812_bundle/main/database'

"""Update PYTHONPATH to be /vol/ek/share/rosetta/rosetta_src_2017.45.59812_bundle/tools/protein_tools/"""

# Read the original peptide sequence
with open(original_seq, 'r') as s:
    ori_seq = s.readline()

# Open the frags_parameters, extract and append parameters to different lists
with open(parameters_file_path, 'r') as f:
    fragments = f.readlines()
pdbs = []
chains = []
start_res = []
end_res = []
sequences = []
for frag in fragments:
    pdbs.append(frag.split()[0])
    chains.append(frag.split()[1])
    start_res.append(frag.split()[2])
    end_res.append(frag.split()[3])
    sequences.append(frag.split()[4])

# Fetch PDBs with prody, call excisePdb and then delete the pdb file
for pdb, chain, start, end, sequence in zip(pdbs, chains, start_res, end_res, sequences):
    outfile = pdb + '.' + chain + '.' + start + '.' + end + '.pdb'
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

    os.remove(pdb_full)

# Create resfile and run fixbb
    resfile = open('resfile', 'w')
    resfile.write('NATRO\nstart')
    for i, res in enumerate(sequence):
        if res == ori_seq[i]:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' NATRO')
        else:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' PIKAA ' + ori_seq[i] +
                          ' EX 1 EX 2')
    resfile.close()

    #fixbb run
    fixbb = '/vol/ek/share/rosetta/rosetta_src_2017.52.59948_bundle//main/source/bin/fixbb.linuxgccrelease' \
            ' -database %s -in:file:s %s -resfile resfile -ex1 -ex2 -use_input_sc -scorefile design_score.sc' \
            ' >design.log' % (rosetta_database, outfile)
    os.system(fixbb)
