import re
import glob

try:
    for filename in glob.glob('frags.500.*mers'):
        frags_handler = open(filename, 'r')
except:
    print "Fragments file doesn't exist\n"

for pdb in frags_handler:
    re.findall('[a-z0-9]{4}') #return a list of strings
    # and write it to the file to the first column, while in the first line write pdb
for chain in frags_handler:
    re.findall('[A-Z]{1}')
    # this one es write to the second column
for start in frags_handler:
    re.findall() #how can I find start and end??

#    1is3 A    13 L L  -48.489  126.715  174.923
#    1is3 A    14 G L   99.634  -10.635  174.479
#    1is3 A    15 M L  -83.802  158.472  176.188
#    1is3 A    16 Y L -119.787  136.568  175.233
#    1is3 A    17 L E -113.031  120.499  172.454
#    1is3 A    18 T E -116.697  124.163  177.320
#
#    1jtd B    65 S L -170.956  164.451  178.838
#    1jtd B    66 G L   84.201    7.782  179.939
#    1jtd B    67 V L  -75.572  131.989 -178.827
#    1jtd B    68 D L -108.504  -26.298  179.569
#    1jtd B    69 A E -151.614  159.915 -178.618
#    1jtd B    70 I E -142.300  159.400  178.150



