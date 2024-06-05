from taxid_class import ncbi_object
from Bio import Entrez
import json
import re
import argparse
from collections import defaultdict

Entrez.email = ""
Entrez.api_key = ""
parser = argparse.ArgumentParser(description='Can fetch taxonomic parents from NCBI taxids and then map the parent of users choice onto a newick tree')
parser.add_argument('-f', '--fasta',help='Any file where the ID lines are in fasta format and contain a taxid eg. OX=id', required=True)
parser.add_argument('-o', '--output', help = 'Output path', required=True)
parser.add_argument('-c', '--count', help = 'Count species in each taxonomic group', default=True)

grouped_args = parser.add_mutually_exclusive_group(required=True)
grouped_args.add_argument('-g', '--get', help = 'True = determining whether software connects to ncbi. USE AN API KEY', action='store_true')
grouped_args.add_argument('-j', '--json', help = 'json file of id_lines + parents (if you dont have this use the g flag first)', type=str)
args = parser.parse_args()
dict_list = []
if args.get:
    infile = args.fasta
    with open(infile, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                new_obj = ncbi_object(line)
                new_obj.fetch()
                if new_obj.scientific_rank:
                    dict_list.append(new_obj.scientific_rank)
    json_name = args.output + ".json"
    with open(json_name, "w") as json_out:
        json.dump(dict_list, json_out)

json_file = json_name if not args.json else args.json

with open(json_file, "r") as json_in:
    file = json.load(json_in)
count_dict = defaultdict(int)
for outer_dict in file:
    #Wll be the ID line
    for key in outer_dict.keys():
        inner_dict = outer_dict[key]
        for val in inner_dict.items():
            count_dict[val] += 1

for key, val in count_dict.items():
    if val != 1:
        print(key, val)
# with open("/home/sjb179/PhD/scripts/parsers/beast_processed.txt", "r") as map:
#     beep = map.readlines()
#     beep = beep[0]
#     for dict in pls:
#         for key in dict.keys():
#             new_dict = dict[key]
#             try:
#                 key = key.split("=")[-1]
#                 reg = f"\w*OX[_| |=]{key}\w*"
#                 key_list = [i for i in new_dict.keys() if not i.startswith("clade")]
#                 beep = re.sub(reg, new_dict[key_list[1]], beep)
#             except Exception as e:
#                 continue
# print(beep)