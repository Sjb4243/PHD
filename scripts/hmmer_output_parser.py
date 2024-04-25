from helper_functions import *
fasta_path = "/home/sjb179/PhD/results/protist_long_frags.txt"
fasta_list = get_fastas(fasta_path)

with open("/home/sjb179/PhD/results/long.txt", "r") as infile:
    infile_lines = infile.readlines()

counter = 0
id_cut_dict = {}
for line in infile_lines:
    if line.startswith(">> tr"):
        id = line.split("|")[1]
        minmax = infile_lines[counter+3].split("..")[1]
        minmax = minmax.strip()
        if len(minmax) > 11:
            minmax = minmax[0:10]
        minmax = minmax.replace("     ", ",")
        id_cut_dict[id] = minmax
    counter +=1


id_list = []
trimmed_list = []
for line in fasta_list:
    if line[0].startswith(">"):
        id = line[0].split("|")[1]
        if id in id_cut_dict.keys():
            fasta_string = "".join(line[1:])
            mini,maxi = id_cut_dict[id].split(",")
            mini = int(mini)
            maxi = int(maxi)
            maxi += 10
            trimmed = fasta_string[mini:maxi]
            id_list.append(line[0])
            trimmed_list.append(trimmed)

from itertools import zip_longest
with open("/home/sjb179/PhD/results/trimmed_longest_fastas.txt", "w") as outfile:
    for id, fasta in zip_longest(id_list, trimmed_list):
        outfile.write(id + "\n")
        outfile.write(fasta + "\n")



