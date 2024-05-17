import Bio.PDB
import numpy as np
import os
import argparse
import subprocess
import ast
resi_dict = {}

def main():
    parser = argparse.ArgumentParser(description='FASTA parser')
    parser.add_argument('-d', '--directory', help='Directory that somewhere contains files ending with model_0.cif')
    parser.add_argument('-r', '--radius', help='Angstrom radius to search for ligands around ATP', type = float)
    global args
    args = parser.parse_args()
    dir = args.directory
    files = [os.path.join(root, file) for root, _, files in os.walk(dir) for file in files if file.endswith("model_0.cif")]
    for file in files:
        filename = os.path.basename(file).split(".")[0]
        cysteines = get_cysteines(file)
        distances = get_distances(cysteines)
        print_outputs(f"DISULPHIDE BOND LOCATIONS FOR {filename} ", distances)
        start_pymol(file)

#I hate working with pymol so we read the pdb and get them directly
def get_cysteines(file):
    resi_dict = {}
    parser = Bio.PDB.MMCIFParser()
    test = parser.get_structure("", file)
    #Get every sulphur atom of every cysteine in the FIRST chain (there arent any cross chain disulphide bridges in p2x)
    for model in test:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if residue.get_resname() == "CYS" and atom.element == "S":
                        resi_dict[residue.get_id()[1]] = atom.coord
            break
        break
    return resi_dict

def get_distances(cysteine_dictionary):
    pairings = []
    for resi1 in cysteine_dictionary:
        for resi2 in cysteine_dictionary:
            # Get the euclidian distance between the two sulphur atoms, if the distance is disulphide size, add it
            if resi1 != resi2 and np.linalg.norm(cysteine_dictionary[resi1]-cysteine_dictionary[resi2]) <= 2.5:
                #Sort the list to be able to filter out duplicate bonds
                pair = sorted([resi1,resi2])
                if pair not in pairings:
                    pairings.append(pair)
    return(pairings)

def print_outputs(header,  list):
    print(header)
    for item in list:
        print(item)

#Need to fix the output eventually but the pymol stuff was a nightmare to get working so it can stay like this for now
def start_pymol(file):
    os.environ["file_name"] = file
    pymol_command = "pymol -cpq PymolWrapper.py"
    process = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    AAoutput = stdout.decode().split("\n")[1]
    # test = ast.literal_eval(AAoutput)
    print(stdout.decode())
    process.wait()
    process.terminate()
    # print(f"ATOMS AROUND ATP for {file}")
    # printed = []
    # for item in test:
    #     if item not in printed:
    #         print(f"{item[0]}-{item[1]}")
    #         printed.append(item)

main()