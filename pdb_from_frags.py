import sys

# Get filename as an argument to the script
frags_file_path = sys.argv[1]

# Read only needed values from frags_file and store them in frags_file_values
parameters_sets = []
frags_parameters = open('frags_parameters', 'w+')
with open(frags_file_path) as frags_file:
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

# Get parameters from first line and last line and save them in values_sets
    j = 0
    for i in range(2, len(all_file_lines)):
        if j % lines_in_block_num == 0:
            first_line_in_set = all_file_lines[i]
            last_line_in_set = all_file_lines[i+4]

            first_line_words = first_line_in_set.split()
            last_line_words = last_line_in_set.split()
            pdb = first_line_words[0]
            chain = first_line_words[1]
            start_res = first_line_words[2]
            end_res = last_line_words[2]
            parameters_sets.append([pdb, chain, start_res, end_res])
        j += 1


for item in parameters_sets:
    for par in item:
        frags_parameters.write("%s " % par)
    frags_parameters.write("\n")

frags_parameters.close()
