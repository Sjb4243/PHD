import os
from pymol import cmd
import sys
from helper_functions import *
import re
import math

#Because Pymol is annoying and iterates through every atom, when you just want residues you have to unduplicate the lists it produces
def unduplicate(aa_List):
    unduped = []
    for item in aa_List:
        if item not in unduped:
            unduped.append(item)
    return unduped
#Dictionary to convert from 3 aa codes to 1
one_letter = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q','ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y','ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A','GLY':'G', 'PRO':'P', 'CYS':'C'}
#Passed from the calling script
cmd.load(os.environ["file_name"])
#Pymol commands to get the first ATP in the .pdb and get either residues around the adenine or phosphate group of the ATP
#Need a way for a user to specify this, probably more os.environs but I'm lazy because this has been a nightmare to get working
cmd.select("firstATP", "chain D")
cmd.select("aden", "firstATP and (name N6 or name N1)")
#cmd.select("aden", "firstATP and (name O2G or name O3G)")
cmd.select("contacts", "aden around 3")
cmd.select("new_contacts", "contacts and not resn ATP")
#Here we get every residue around atp in a certain radius, and deduplicate it
aroundATP = []
cmd.iterate("new_contacts", "aroundATP.append([resi, resn, chain])")
aroundATP = unduplicate(aroundATP)

#Spacing is used to get n amount of residues either side of our target one in aroundATP
spacing = 6
AA_strings = []
for resi in aroundATP:
    residues = []
    num = int(resi[0])
    low, high  = num - spacing, num + spacing
    #Get a range of residues as defined by spacing
    cmd.select("all_aa_in_range", f"chain A and resi {low}-{high}")
    cmd.iterate("all_aa_in_range", "residues.append([resi,resn])")
    #Again we have to deduplicate this due to pymol
    new_list = unduplicate(residues)
    #After getting the corresponding 1 letter code, append to a string and store as a list
    AA_strings.append(''.join(one_letter[item[1]] for item in new_list))
#Lots of debugging print statements, its nice to see whats going on but I should tidy these up at some point
print(aroundATP)
sys.stdout.flush()
print(AA_strings)
sys.stdout.flush()

#Get the original alignment file and convert into a 2d list of seqIDs and fasta strings
fasta_list = get_fastas("/home/sjb179/PhD/msa_fastas/second_alignment/separated_fastas/new_alignments/gradual_msa/alphafold/separated_alignment_copy.fa")
#Used for forming a regex to match seqIDs
file_name = os.environ["file_name"]
file_name = file_name.split("/")[-2]
file_name_for_regex = file_name[5:-4]

out_list = []
for item in fasta_list:
    if re.search(file_name_for_regex, item[0], flags=re.IGNORECASE):
        for aa_string in AA_strings:
            target_amino = aa_string[math.floor(len(aa_string)/2)]
            aa_regex = []
            #Construct a regex to find the sequence separated by any number of gaps
            aa_regex = ''.join(amino + "[-]*" for amino in aa_string)
            #More debug print statements, its nice to see matches + the regex string itself, will remove eventually
            print(aa_regex)
            sys.stdout.flush()
            print(aa_string)
            sys.stdout.flush()
            print(re.search(aa_regex, item[1], flags=re.IGNORECASE))
            sys.stdout.flush()
            #First match the regex in the total fasta string, then replace any of the target letters in the matched group with the letter but lowercase
            #It's almost impossible to find which letter it is in the target string we actually want to be lowercase so we do them all
            #It isnt too disruptive so it has to stay like this for now (I will probably never fix this)
            item[1] = re.sub(aa_regex, lambda m: m.group(0).replace(target_amino, target_amino.lower()),item[1], flags=re.IGNORECASE)
        out_list.append([item[0], item[1]])
#I could make this into a write (which i might do) by returning outlist and then passing it back and forth and writing in the other script but
#It's kinda painful doing this with subprocess so it can stay like this for now
with open("/home/sjb179/output_file.txt", "a") as output_file:
    for item in out_list:
        output_file.write(item[0] + "\n")
        output_file.write(item[1] + "\n")
cmd.quit()
