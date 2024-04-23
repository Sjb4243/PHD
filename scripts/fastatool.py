import argparse
import matplotlib.pyplot as plt
from itertools import zip_longest
from helper_functions import *
from collections import defaultdict
import statistics
import re
class fasta_obj:
    def __init__(self, id, species, protein, fasta_string, full_id):
        self.id = id
        self.species = species
        self.protein = protein
        self.fasta_string = "".join(fasta_string)
        self.full_id = full_id


def main():
    parser = argparse.ArgumentParser(description='FASTA parser')
    parser.add_argument('-i', '--ID', help='directory to a file containing uniprot IDs separated by newline')
    parser.add_argument('-f', '--fasta', help='directory to a fasta file', required=True)
    parser.add_argument('-p', '--plots', help='True/false for if you want plots', default=False)
    parser.add_argument('-s', '--search', help='Search the fasta line for text (supports regex)')
    parser.add_argument('-o', '--output', help='Output file containing curated fastas')
    parser.add_argument('-c', '--cutoff', help = 'Length cutoff for fragments, has a minimum and maximum length, sep by comma', default=[0,200000000])

    args = parser.parse_args()
    fasta_list = get_fastas(args.fasta)
    input_id_list = []
    if args.ID:
        with open(args.ID, "r") as idfile:
            input_id_list = [line.strip() for line in idfile]
    obj_list = []
    for fasta in fasta_list:
        id_line = fasta[0]
        if args.search:
            if re.search(r"^>" + ".*" + args.search, id_line, flags=re.IGNORECASE):
                obj_list.append(create_obj(fasta, input_id_list, args.cutoff))
        else:
            if id_line.startswith(">"):
                obj_list.append(create_obj(fasta, input_id_list, args.cutoff))
    length_list, species_dict = print_summary(obj_list)
    if args.output:
        save_output(obj_list, args.output)
    if args.plots:
        print_plots(length_list, species_dict)

def save_output(fasta_objects, output_filename):
    id_lines = []
    fasta_strings = []
    for item in fasta_objects:
        if item:
            id_lines.append(item.full_id)
            fasta_strings.append(item.fasta_string)
    with open(output_filename, "w") as output_fasta:
        for id, fasta in zip_longest(id_lines, fasta_strings):
            output_fasta.write(id + "\n")
            output_fasta.write(fasta + "\n")
    with open(output_filename + "_ids", "w") as id_output:
        for id in id_lines:
            id_output.write(id.split("|")[1] + "\n")


def print_summary(fasta_objects):
    species_dict = defaultdict(int)
    length_list = []
    for item in fasta_objects:
        if item:
            species_dict[item.species] += 1
            length_list.append(len(item.fasta_string))
    print("#######################SPECIES#######################")
    for key,val in species_dict.items():
        print(f"{key}: {val}")
    print("#######################SIZE#######################")
    print(f"Mean length of fragment size is: {statistics.mean(length_list)}")
    return length_list, species_dict

def print_plots(length, species):
    plt.figure()
    plt.hist(length, bins=10, color='blue', edgecolor='black')
    plt.title("Fragment length distribution")

    species_keys = list(species.keys())
    values = list(species.values())

    plt.figure()
    plt.bar(species_keys, values, color="skyblue")
    plt.title("Species distribution")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def create_obj(fasta_entry, id_list, length_cutoff):
    id_line = fasta_entry[0]
    data = extract_data(id_line)
    id, species, protein = data[0], data[1], data[2]
    fasta_lines = fasta_entry[1:]
    mini,maxi = length_cutoff.split(",")
    mini = int(mini)
    maxi = int(maxi)
    print(maxi)
    if len("".join(fasta_lines)) > mini and len("".join(fasta_lines)) < maxi:
        if len(id_list) > 0:
            if id in id_list:
                finished_obj = fasta_obj(id, species, protein, fasta_lines, fasta_entry[0])
            else:
                return None
        else:
            finished_obj = fasta_obj(id, species, protein, fasta_lines, fasta_entry[0])
        return finished_obj

main()