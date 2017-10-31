import sys

# Get filename as an argument to the script
frags_file_path = sys.argv[1]

# Read only needed values from frags_file and store them in frags_file_values
parameters_sets = []
frags_parameters = open('frags_parameters', 'w+')
with open(frags_file_path) as frags_file:
    all_file_lines = frags_file.readlines()
    j = 0
    for i in range(2, len(all_file_lines)):
        if j % 6 == 0:
            first_line_in_set = all_file_lines[i]
            last_line_in_set = all_file_lines[i+4]

            # Get parameters from first line and last line and save them in values_sets
            first_line_words = first_line_in_set.split()
            last_line_words = last_line_in_set.split()

            pdb = first_line_words[0]
            chain = first_line_words[1]
            start_res = first_line_words[2]
            end_res = last_line_words[2]
            # write the
            parameters_sets.append([pdb, chain, start_res, end_res])
        j += 1
        
print(parameters_sets)

for item in parameters_sets:
    frags_parameters.write(item, '\n')
frags_parameters.readlines()
print(frags_parameters)
frags_parameters.close()


