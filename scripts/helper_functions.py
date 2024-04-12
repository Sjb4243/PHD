import re
#Series of helper functions for general parsing of swissprot files

def get_fastas(fasta_file):
    """
    Obtain a list of fastas from a fasta file
    Parameters:
    param (type): fasta file

    Returns:
    type: List of fasta files + a header
    """
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

def get_species_id_name(fasta_file):
    """
    :param fasta_file:
    :return: Returns a dictionary of IDs plus the species names
    """
    trembl_ids = {line.split("|")[1]: re.findall("(?<=OS=)\w* \w*", line.split("|")[2]) for line in fasta_file if line.startswith(">")}
    return trembl_ids


def extract_data(fasta_header):
    """
    :param fasta_header: fasta header line (uniprot)
    :return: important information
    """
    split_header = fasta_header.split("|")
    reg_split = re.split(r"\||_", fasta_header)
    id = split_header[1]
    species_name = re.findall(r"(?<=OS=)\w* \w*", split_header[2])[0]
    protein = re.findall(r"(?<=.{5} ).*(?= OS=)", reg_split[3])
    if not protein:
        protein = "UNKNOWN"
    else:
        protein = protein[0]
    return([id,species_name, protein])