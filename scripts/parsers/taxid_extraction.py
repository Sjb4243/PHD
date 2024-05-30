import argparse
import os
import re
from Bio import Entrez
import time


def main():
    parser = argparse.ArgumentParser(description='HMMER output parser')
    parser.add_argument('-d', '--directory', help='directory containing files', required=True)
    args = parser.parse_args()
    dir = args.directory
    for file in os.listdir(dir):
        if is_nr(file):
            fullpath = dir + file
            ncbi_extraction(fullpath)
        else:
            swissprot_extraction(dir + file)


def is_nr(infile):
    if "nr.gz" in infile:
        return True


def get_existing_ids():
    with open("output.txt", "r") as output_file:
        ids = [line.split("\t")[0] for line in output_file]
    return ids


def swissprot_extraction(infile):
    print(infile)
    ids = get_existing_ids()
    with open(infile, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                segline = re.split(r".{2}=", line)
                id = segline[0].split("|")[1]
                print("Writing to file:", id, segline[1], segline[2])
                write_file(id, segline[1], segline[2])


def ncbi_extraction(infile):
    Entrez.email = ""
    Entrez.api_key = ""
    ids = get_existing_ids()
    with open(infile, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                spec_id = line.strip(">").split(" ")[0]
                if spec_id in ids:
                    print(f"Skipping {spec_id}")
                    continue
                name = re.search(r"\[[a-z A-Z]*\]", line)
                if name:
                    name = name.group()
                else:
                    name = "UNKNOWN"
                try:
                    request = Entrez.efetch(db="protein", id=spec_id, rettype="gb", retmode="text")
                    xml = request.read()
                    taxon = re.findall(r"taxon:\d*", xml)[0]
                    taxon = taxon.replace("taxon:", "")
                except:
                    taxon = "UNKNOWN"
                write_file(spec_id, name, taxon)
                print(f"Species ID: {spec_id} \n Species name: {name} \n Species taxon: {taxon}")


def write_file(species_id_set, species_name_set, taxid_set):
    with open('output.txt', 'a') as output:
        output.write(species_id_set + "\t" + species_name_set + "\t" + taxid_set + "\n")


main()
