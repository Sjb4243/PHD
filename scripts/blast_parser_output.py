file = "/home/sjb176/sprot_output.txt"
import re

with open(file, "r") as f:
    for line in f:
        full_file = f.readlines()

counter = 0
query_list = []
hit_list = []
for line in full_file:
    if line.split("\t")[0] not in query_list:
        query_list.append(line.split("\t")[0])
        hit = full_file[counter:counter + 5]
        full_protein_string = ""
        full_hit_string = ""
        for i in hit:
            temp = i.split("\t")[1]
            score  = i.split("\t")[2]
            protein = re.findall(r".*(?=OS=)", temp)[0]
            species = re.findall(r"(?<=OS=).*(?= OX)", temp)[0]
            full_hit_string += species + "|"
            full_protein_string += protein[10:] + "__" + score + "__"
        test = line.split("\t")[0]
        with open("output.txt", "a") as output:
            output.write(test + "\t" + "\t" + full_protein_string + "\n")
        counter += 1

with open("output.txt", "r") as checkfile1:
    file1 = checkfile1.readlines()
with open("/home/sjb176/PhD/results/procced_for_blast.txt", "r") as checkfile2:
    file2 = checkfile2.readlines()

counter = 0
species_from_input = [line.split("\t")[0] for line in file2]
new_species = [line.split("|")[1] for line in species_from_input if line.startswith(">")]

out_species = []
for line in file1:
    species_from_output = line.split("\t")[0]
    species_from_output = species_from_output.split("|")[1]
    out_species.append(species_from_output)

for line in new_species:
    if line not in out_species:
        print("AAAAAAA")