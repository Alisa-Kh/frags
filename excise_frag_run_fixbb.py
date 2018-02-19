import os
import sys

ROSETTA_DIR = '/vol/ek/share/rosetta/'
ROSETTA_VERSION = 'rosetta_src_2017.45.59812_bundle'
ROSETTA_DATABASE = ROSETTA_DIR + ROSETTA_VERSION + '/main/database'

# Commands

PRODY = '/vol/ek/share/csbw/tools/env2017/bin/prody fetch %s'  # used to fetch PDB, but any other method can be used
RENUMBERING = ROSETTA_DIR + ROSETTA_VERSION + '/tools/protein_tools/scripts/pdb_renumber.py {} {}'
EXCISE_PDB = ROSETTA_DIR + '/pdbUtil/excisePdb_v2.pl {} {} {} {} {}'
FIXBB = ROSETTA_DIR + ROSETTA_VERSION + '/main/source/bin/fixbb.linuxgccrelease' \
        ' -database %s -in:file:s %s -resfile resfile_tmp -ex1 -ex2 -use_input_sc' \
        ' -scorefile design_score.sc -ignore_zero_occupancy false >>design.log'

# Get filename as an argument to the script
parameters_file_path = sys.argv[1]
original_seq = sys.argv[2]  # peptide_sequence

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


def delete_frag(fragment):  # delete bad fragments
        os.remove(fragment)
        print("Wrong length fragment has been deleted")

# Fetch PDBs (here done with prody), call excisePdb and then delete the pdb file
for pdb, chain, start, end, sequence in zip(pdbs, chains, start_res, end_res, sequences):
    outfile = pdb + '.' + chain + '.' + start + '.' + end + '.pdb'
    pdb_full = pdb + '.pdb'

    os.system(PRODY % pdb)
    os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))

    if not os.path.exists(outfile):
        new_pdb = pdb_full
        print("renumbering:")
        os.system(RENUMBERING.format(pdb_full, new_pdb))  # Renumber PDB if it doesn't start from 1

        os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))
    os.remove(pdb_full)

# check correctness of fragments
    with open(outfile) as frag:
        bad_pdb = False
        residues = set()
        cur_line = frag.readline().split()
        while cur_line[0] != 'ATOM':  # find ATOM lines
            cur_line = frag.readline().split()

            # if there is no atoms...
            if not cur_line:
                bad_pdb = True
                break
        if bad_pdb:
            delete_frag(outfile)
            continue

        cur_line = frag.readline()
        while cur_line[22:27].strip() != start:
            cur_line = frag.readline()

            # if there is no start residue
            if not cur_line:
                bad_pdb = True
                break
        if bad_pdb:
            delete_frag(outfile)
            continue

        cur_line = frag.readline()
        while cur_line[22:27].strip() != end:
            residues.add(cur_line[22:27])
            cur_line = frag.readline()

            # if there is no end residue
            if not cur_line:
                bad_pdb = True
                break
        if bad_pdb:
            delete_frag(outfile)
            continue

        residues.add(cur_line[22:27])
    if len(residues) != (int(end) - int(start) + 1):
        delete_frag(outfile)
        continue

# Create resfile and run fixbb
    resfile = open('resfile_tmp', 'w')
    resfile.write('NATRO\nstart')
    if chain == '_':
        chain = 'A'
    for i, res in enumerate(sequence):
        if res == ori_seq[i]:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' NATRO')
        else:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' PIKAA ' + ori_seq[i] +
                          ' EX 1 EX 2')
    resfile.close()

    # fixbb run
    print("running fixbb design")
    os.system(FIXBB % (ROSETTA_DATABASE, outfile))

    os.remove('resfile_tmp')
    os.remove(outfile)
