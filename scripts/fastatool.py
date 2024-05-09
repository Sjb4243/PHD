#!/usr/bin/env python3
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
    parser.add_argument('-c', '--cutoff', help = 'Length cutoff for fragments, has a minimum and maximum length, sep by comma', default="200,200000000")
    parser.add_argument('-m', '--metadata', help= 'If writing to output, create a metadata file', default=False)
    parser.add_argument('-d', '--diff', help = 'Additional fasta file  to compare. Will return every fasta not in the second file but in the first.')
    args = parser.parse_args()
    #Uses the get_fastas function from helper_functions.py
    #This returns a 2d list, where each element consists of a list of the ID line followed by the fasta strings
    fasta_list= get_fastas(args.fasta)
    fasta_list = remove_duplicates(fasta_list)
    if args.diff:
        second_fasta = get_fastas(args.diff)
        second_fasta = remove_duplicates(second_fasta)
    fasta_list = get_difference(fasta_list, second_fasta)
    #Removes duplicate ID lines
    input_id_list = []
    #If the user wants to filter by IDs, generate a list of IDs based on an input file
    if args.ID:
        with open(args.ID, "r") as idfile:
            input_id_list = [line.strip() for line in idfile]
    #Init a list for storing fasta objects
    obj_list = []
    for fasta in fasta_list:
        id_line = fasta[0]
        #Loop through the fastas - If the user has a string to search with, append that to a regex
        #This regex could probably be made better but for every usecase so far it's proved to be enough
        if args.search:
            if re.search(r"^>" + ".*" + args.search, id_line, flags=re.IGNORECASE):
                obj_list.append(create_obj(fasta, input_id_list, args.cutoff))
        else:
            if id_line.startswith(">"):
                obj_list.append(create_obj(fasta, input_id_list, args.cutoff))
    #Sends off the species + amino acid length to be printed to terminal
    length_list, species_dict = print_summary(obj_list)
    #Saving files
    if args.output:
        save_output(obj_list, args.output, args)
    #Plotting output
    if args.plots:
        print_plots(length_list, species_dict)

def get_difference(first_fasta, second_fasta):
    second_id_lines = [line[0] for line in second_fasta]
    difference = [line[0:] for line in first_fasta if line[0] not in second_id_lines]
    print(difference)
    return difference

def remove_duplicates(fasta_list):
    #init a set - now that I think about it this may not need to be a set
    seen = set()
    non_dupes = []
    for line in fasta_list:
        #If ID isnt in the set append it to another list
        if line[0] not in seen:
            non_dupes.append(line)
            seen.add(line[0])
    return non_dupes

def save_output(fasta_objects, output_filename, args):
    #Init empty lists
    id_lines = []
    fasta_strings = []
    #Loop through all fasta objects and append the fasta strings to a list
    for item in fasta_objects:
        if item:
            data = extract_data(item.full_id)
            id, species, protein, taxid = data[0], data[1], data[2], data[3]
            new_id = ">tr|" + id + "|" + protein + "|" + "OS=" + species + "|OX=" + taxid
            id_lines.append(new_id)
            fasta_strings.append(item.fasta_string)
    #Zip longest interleaves the lists which makes it easier to write
    with open(output_filename, "w") as output_fasta:
        for id, fasta in zip_longest(id_lines, fasta_strings):
            output_fasta.write(id + "\n")
            output_fasta.write(fasta + "\n")
    #Write all the IDs to a file for further checking
    #This has become less necessary as time has passed so I may remove this eventually
    with open(output_filename + "_ids", "w") as id_output:
        for id in id_lines:
            id_output.write(id.split("|")[1] + "\n")
    if args.metadata:
        with open(output_filename + "_metadata", "w") as meta_data_output:
            #using default dict because it handles the initial addition to the dictionary
            species_dict = defaultdict(int)
            full_fasta_strings = []
            #append each species and fasta string to dictionary and list respectively
            for item in fasta_objects:
                if item:
                    species_dict[item.species] += 1
                    full_fasta_strings.append(item.fasta_string)
            #Write all the info out to the metadata file
            meta_data_output.write(">Species" + "\n")
            for key,val in species_dict.items():
                meta_data_output.write(f"{key}:{val}" + "\n")
            meta_data_output.write(">Length" + "\n")
            mini,maxi = args.cutoff.split(",")
            meta_data_output.write(f"Cutoff:{mini},{maxi}" + "\n")
            meta_data_output.write(f"Smallest/largest frag sizes:{len(min(full_fasta_strings, key=len))},{len(max(full_fasta_strings, key=len))}")


def print_summary(fasta_objects):
    #Fairly simple printing to terminal
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
    try:
        print(f"Mean length of fragment size is: {statistics.mean(length_list)}")
    except:
        print("Not enough data for mean!")
    return length_list, species_dict

def print_plots(length, species):
    #Simple matplotlib plots
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
    #Uses the extract data function from helper_functions
    #It just uses some regexes on the ID line eo extract the needed data
    data = extract_data(id_line)
    id, species, protein, taxid = data[0], data[1], data[2], data[3]
    fasta_lines = fasta_entry[1:]
    #split the cutoff
    mini,maxi = length_cutoff.split(",")
    if len("".join(fasta_lines)) > int(mini) and len("".join(fasta_lines)) < int(maxi):
        if len(id_list) > 0:
            if id in id_list:
                finished_obj = fasta_obj(id, species, protein, fasta_lines, fasta_entry[0])
            else:
                return None
        else:
            finished_obj = fasta_obj(id, species, protein, fasta_lines, fasta_entry[0])
        return finished_obj

main()