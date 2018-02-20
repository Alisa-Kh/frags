import sys
import os

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


def create_params_file(frags):

    # Read only needed values from frags_file and store them in frags_file_values
    parameters_sets = []
    frags_parameters = open('frags_parameters', 'w+')
    with open(frags) as frags_file:
        all_file_lines = frags_file.readlines()

    # Count a number of lines in one block (one fragment)
    first_blank_line = None
    second_blank_line = None
    for i in range(len(all_file_lines)):
        if all_file_lines[i] in ['\n', '\r\n']:
            if first_blank_line is None:
                first_blank_line = i
            elif second_blank_line is None:
                second_blank_line = i
                break

    lines_in_block_num = second_blank_line - first_blank_line
    # Get parameters and save them in parameters_sets list
    j = 0
    for i in range(2, len(all_file_lines)):
        if j % lines_in_block_num == 0:
            first_line_in_set = all_file_lines[i]
            last_line_in_set = all_file_lines[i+(lines_in_block_num-2)]
            seq = ""
            for k in range(i, i+lines_in_block_num-1):
                line = all_file_lines[k].split()
                seq += line[3]
            first_line_words = first_line_in_set.split()
            last_line_words = last_line_in_set.split()
            pdb = first_line_words[0]
            chain = first_line_words[1]
            start_res = first_line_words[2]
            end_res = last_line_words[2]
            parameters_sets.append([pdb, chain, start_res, end_res, seq])
        j += 1

    for item in parameters_sets:
        for par in item:
            frags_parameters.write("%s\t" % par)
        frags_parameters.write("\n")
    frags_parameters.close()


def delete_frag(fragment):  # delete bad fragments
    os.remove(fragment)
    print("Wrong length fragment has been deleted")


def fixbb_design(ori_seq, chain, start, sequence, outfile):

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


def review_frags(outfile, start, end):
    # check correctness of fragments
    with open(outfile) as frag:
        residues = set()
        cur_line = frag.readline().split()
        while cur_line[0] != 'ATOM':  # find ATOM lines
            cur_line = frag.readline().split()

            # if there is no atoms...
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()
        while cur_line[22:27].strip() != start:
            cur_line = frag.readline()

            # if there is no start residue
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()
        while cur_line[22:27].strip() != end:
            residues.add(cur_line[22:27])
            cur_line = frag.readline()

            # if there is no end residue
            if not cur_line:
                delete_frag(outfile)
                return False

        residues.add(cur_line[22:27])
    if len(residues) != (int(end) - int(start) + 1):
        delete_frag(outfile)
        return False

    else:
        return True


def extract_fragments(orig_seq):

    # Read the original peptide sequence
    with open(orig_seq, 'r') as s:
        ori_seq = s.readline()

    # Open the frags_parameters, extract and append parameters to different lists
    with open('frags_parameters', 'r') as f:
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

        is_frag_ok = review_frags(outfile, start, end)

        if is_frag_ok:
            fixbb_design(ori_seq, chain, start, sequence, outfile)


if __name__ == "__main__":

    # Get filename as an argument to the script
    frags_file_path = sys.argv[1]
    create_params_file(frags_file_path)

    original_seq = sys.argv[2]  # peptide_sequence

    extract_fragments(original_seq)
