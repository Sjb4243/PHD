import re
#Series of helper functions for general parsing of swissprot files

def get_fastas(fasta_file):
    """
    Obtain a list of fastas from a fasta file
    Parameters:
    param (type): fasta file DIRECTORY

    Returns:
    type: List of fasta files + a header
    """
    with open(fasta_file, "r") as infile:
        fasta_file = infile.readlines()
    skipped = 0
    fasta_list = []
    for i in range(len(fasta_file)):
        index = i + skipped
        if index >= len(fasta_file):
            break
        if fasta_file[index].startswith(">"):
            done = False
            newPos = index + 1
            temp_fastas = []
            temp_fastas.append(fasta_file[index].strip())
            runs = 0
            while not done:
                if newPos + 1 > len(fasta_file):
                    fasta_list.append(temp_fastas)
                    break
                if fasta_file[newPos].startswith(">"):
                    fasta_list.append(temp_fastas)
                    skipped = skipped + runs
                    done = True
                if not done:
                    temp_fastas.append(fasta_file[newPos].strip())
                    newPos += 1
                    runs += 1
    return fasta_list

def extract_data(fasta_header):
    """
    :param fasta_header: fasta header line (uniprot)
    :return: important information
    """
    split_header = fasta_header.split("|")
    reg_split = re.split(r"\||_", fasta_header)
    id = split_header[1]
    species_name = re.findall(r"(?<=OS=).*(?=[ |_]OX)", split_header[2])[0]
    protein = re.findall(r"(?<=[0-9][A-Z]{4}_)\w*(?=[\||_]OS)|(?<=\|)\w*(?=\|OS)", fasta_header)
    taxid = re.findall(r"(?<=OX=)[0-9]*", fasta_header)[0]
    if not protein:
        protein = "UNKNOWN"
    else:
        protein = protein[0]
    return([id,species_name, protein, taxid])

def make_id(*comps):
    empt = ""
    for arg in comps:
        empt += arg + "_"
    empt = empt[:-1]
    empt += "\n"
    return empt