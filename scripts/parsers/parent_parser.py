import pandas as pd
from collections import defaultdict
import ast
import re
import argparse

parser = argparse.ArgumentParser(
    description='Processes a taxid + parent csv file and modifies the species names on a newick tree based on common ancestry ')
parser.add_argument('-i', '--infile',
                    help='File where the 2nd column has taxids and the third column has a list of parents for that taxid')
parser.add_argument('-n, --newick', help='Newick tree')
parser.add_argument('-o', '--outfile', help='Output newick tree')
parser.add_argument('-c', '--cutoff',
                    help='Negative number, used to define how far back you want to move in the tree of life. -6 = cellular organisms etc, -5 is more resolution etc',
                    default=-4, type=int)

args = parser.parse_args()

df = pd.read_csv(args.infile)
unique_df = df.drop_duplicates(subset="taxids")

parent_dict = defaultdict(list)
#Make a dictionary where each key is a component of the evolutionary tree
#And the value is a list of taxids that have that key as a parent
for index, row in unique_df.iterrows():
    parent = row["parents"]
    taxid = row['taxids']
    #Used to convert the string list in the csv to an object that is actually a list
    parents = ast.literal_eval(parent)
    for item in parents:
        parent_dict[item].append(taxid)


taxid_check = []
taxid_dict = defaultdict(list)
#Using the dictionary generated above, generate a dictionary where:
#Each key is a taxid
#The value consists of a list of lists where each list is a list of parents of each taxid from the dictionary above

for clade, taxid_list in reversed(parent_dict.items()):
    #We now want to get every key to which a taxid belongs
    for taxid in taxid_list:
        if taxid not in taxid_check:
            key_list = [key for key, val in reversed(parent_dict.items()) if taxid in val]
            taxid_dict[taxid].append(key_list)
            taxid_check.append(taxid)
final_dict = {}
for key, val in taxid_dict.items():
    initialisation_set = set(val[0])
    common_list = []
    for sublist in val[1:]:
        initialisation_set &= set(sublist)
        common_list = [element for element in val[0] if element in initialisation_set]
    if not common_list:
        common_list = val[0]
    final_dict.update({key:common_list[args.cutoff]})

with open("newick_test", "r") as file:
    beep = file.readlines()
beep = beep[0]
pattern = r'(tr.*?)(?=OX)'
#beep = re.sub(pattern, "", beep)
for key, val in final_dict.items():
    find_str = "OX_" + str(key)
    beep = beep.replace(find_str, val)
    beep = beep.replace(" ", "")
    beep = beep.replace("(in:glomeromycetes)", "")

with open(args.outfile, "w") as outfile:
    outfile.write(beep)
