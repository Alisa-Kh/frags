import sys
import os

ROSETTA_DIR = '/vol/ek/share/rosetta/'
ROSETTA_VERSION = 'rosetta_src_2017.45.59812_bundle'
ROSETTA_DATABASE = ROSETTA_DIR + ROSETTA_VERSION + '/main/database'
PARAMS_FILENAME = 'frags.100.{}mers'

# Commands

PRODY = '/vol/ek/share/csbw/tools/env2017/bin/prody fetch %s'  # used to fetch PDB, but any other method can be used
RENUMBERING = ROSETTA_DIR + ROSETTA_VERSION + '/tools/protein_tools/scripts/pdb_renumber.py {} {}'
EXCISE_PDB = ROSETTA_DIR + '/pdbUtil/excisePdb_v2.pl {} {} {} {} {}'
FIXBB = ROSETTA_DIR + ROSETTA_VERSION + '/main/source/bin/fixbb.linuxgccrelease' \
                                        ' -database %s -in:file:s %s -resfile resfile_tmp -ex1 -ex2 -use_input_sc' \
                                        ' -scorefile design_score.sc -ignore_zero_occupancy false >>design.log'

BUILD_PEPTIDE = ROSETTA_DIR + ROSETTA_VERSION + '/main/source/bin/BuildPeptide.linuxgccrelease -in:file:fasta {}' \
                                                ' -database ' + ROSETTA_DATABASE + ' -out:file:o peptide.pdb ' \
                                                                                   '> build_peptide.log'
MAKE_FRAGMENTS = 'perl make_fragments.pl -verbose -id xxxxx {} 2>log'
FRAG_PICKER = ROSETTA_DIR + ROSETTA_VERSION + '/main/source/bin/fragment_picker.linuxgccrelease' \
              ' -database ' + ROSETTA_DATABASE + ' @flags >makeFrags.log'


def make_pick_fragments(pep_seq):
    # if not os.path.exists('resfiles'):
    #     os.makedirs('frag_picker')
    with open('xxxxx.fasta', 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')

    os.system(MAKE_FRAGMENTS.format('xxxxx.fasta'))
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n'
                     '#SequenceIdentity\t150\t0.0\t-\n'
                     '#RamaScore\t100\t6.0\t-\tpsipred\n'
                     '#FragmentCrmsd\t30\t0.0\t-\n'
                     '#FragmentAllAtomCrmsd\t20\t0.0\t-')
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + ROSETTA_DIR + 'rosetta_fragments_latest/nnmake_database'
                                                           '/vall.dat.2006-05-05\n')
        flags_file.write('-in:file:checkpoint\txxxxx.checkpoint\n')
        flags_file.write('-frags:describe_fragments\tfrags.fsc\n')
        flags_file.write('-frags:frag_sizes\t' + str(len(pep_seq)) + '\n')
        flags_file.write('-frags:n_candidates\t2000\n')
        flags_file.write('-frags:n_frags\t100\n')
        flags_file.write('-out:file:frag_prefix\tfrags\n')
        flags_file.write('-frags:ss_pred\txxxxx.psipred_ss2 psipred\n')
        flags_file.write('-frags:scoring:config\tpsi_L1.cfg\n')
        flags_file.write('-frags:bounded_protocol\ttrue\n')
        flags_file.write('-mute\tcore.util.prof\n')
        flags_file.write('-mute\tcore.conformation\n')
        flags_file.write('-mute\tcore.chemical\n')
        flags_file.write('-mute\tprotocols.jumping')
    os.system(FRAG_PICKER)


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


def extract_fragments(pep_sequence):

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

        if os.path.exists(pdb_full):
            os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))
        else:
            continue  # if failed to fetch PDB

        if not os.path.exists(outfile):
            new_pdb = pdb_full
            print("renumbering:")
            os.system(RENUMBERING.format(pdb_full, new_pdb))  # Renumber PDB if it doesn't start from 1

            os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))
        os.remove(pdb_full)

        is_frag_ok = review_frags(outfile, start, end)

        if is_frag_ok:
            fixbb_design(pep_sequence, chain, start, sequence, outfile)


def create_resfile(ori_seq, chain, start, sequence):

    if not os.path.exists('resfiles'):  # create directory for storing resfiles
        os.makedirs('resfiles')

    path_to_resfile = 'resfiles/resfile_%s'

    # Create resfile for each fragment

    resfile = open(path_to_resfile % sequence, 'w')
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
    return resfile


def fixbb_design(ori_seq, chain, start, sequence, outfile):

    if not os.path.exists('resfiles'):  # create directory for storing resfiles
        os.makedirs('resfiles')
    path_to_resfile = 'resfiles/resfile_%s'

    # Create resfile for each fragment
    resfile = open(path_to_resfile % sequence, 'w')
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

    # os.remove('resfile_tmp')
    os.remove(outfile)


def build_peptide(pep_seq):

    # Build extended peptide
    os.system(BUILD_PEPTIDE.format(pep_seq))

    # Change chain ID to 'B'
    renamed_peptide = []
    with open('peptide.pdb', 'r') as peptide:
        pdb_lines = peptide.readlines()
    for line in pdb_lines:
        if line[21].isalpha():
            new_line = list(line)
            new_line[21] = 'B'
            renamed_peptide.append("".join(new_line))
    os.remove('peptide.pdb')
    with open('peptide.pdb', 'w') as new_peptide:
        for bchain_line in renamed_peptide:
            new_peptide.write(bchain_line)


if __name__ == "__main__":

    with open(sys.argv[1], 'r') as peptide:
        peptide_seq = peptide.readline().strip()

    make_pick_fragments(peptide_seq)

    # extract parameters from fragment picker output for extracting fragments
    create_params_file(PARAMS_FILENAME.format(str(len(peptide_seq))))

    extract_fragments(peptide_seq)  # extract fragments, create resfiles and run fixbb design

    build_peptide(peptide_seq)  # build extended peptide and change it's chain id to 'B'
