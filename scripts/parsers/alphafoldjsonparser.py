import json
import os
import re
import argparse
class alphafold_obj:
    def __init__(self, summary_jsonpath, full_json_path, cifpath, target_atoms):
        self.cifpath = cifpath
        self.summary_jsonfile = json.load(open(summary_jsonpath, "r"))
        self.full_jsonfile = json.load(open(full_json_path, "r"))
        self.target_atoms = target_atoms
        self.ciffile = open(cifpath, "r").readlines()
        self.name = self.get_name()
        self.summary_confidence = self.get_summary_data()
        self.atom_confidence = self.get_atom_confidence()

    def get_name(self):
        file_name = re.findall("(?<=fold_).*(?=_model)", self.cifpath.split("/")[-1])[0].replace("3atp", "").strip("_")
        print(f"#################################{file_name}##############################")
        return file_name

    def get_summary_data(self):
        chain_confidence = self.summary_jsonfile["chain_iptm"]
        print(f"ATP confidence is {chain_confidence[-1]} vs human ATP's confidence of 0.66")
        return chain_confidence

    def get_atom_confidence(self):
        pdb_found = False
        pdb_file = []
        target_atoms_indexes = {}
        for line in self.ciffile:
            #Capture pdb part
            if pdb_found and line.strip() != "#":
                pdb_file.append(line)
                #For every line, check if any of the target atoms are in the current line
                for name in self.target_atoms:
                    if name in line:
                        #Get the first instance of each name
                        if name not in target_atoms_indexes.keys():
                            #Regex to find residue position(each ATP atom is considered a residue)
                            search = re.search(r"HETATM ([0-9]*)", line)
                            match = search.group(1)
                            #Construct dictionary of ATP atoms: resi index
                            target_atoms_indexes[name] = int(match)
            #Start of PDB line
            if line.startswith("_atom_site.pdbx_PDB_model_num"):
                pdb_found = True
        #Using the residue index, get the confidence score from the atom_plddts index
        for key, val in target_atoms_indexes.items():
            conf_score = self.full_jsonfile["atom_plddts"][val]
            print(f"({val})    Atom:{key} --- Confidence score: {conf_score}")
        return None




parser = argparse.ArgumentParser(description='Parses alphafold files and reports the general ligand confidence (assuming trimer and 3 ligands) Can also report individual atom confidences')
parser.add_argument('-d', '--directory', help='Start directory, will walk recursively from here and find all model_0 alphafold files')
parser.add_argument('-t', '--targets', help='Target atom names separated by comma', type = str, default=["O2G", "O3G", "N6", "N1"])
args = parser.parse_args()
dir = args.directory
targets = args.targets.split(",")
full_json_files = [os.path.join(root, file) for root, _, files in os.walk(dir) for file in files if file.endswith("full_data_0.json")]
summary_json_files = [os.path.join(root, file) for root, _, files in os.walk(dir) for file in files if file.endswith("summary_confidences_0.json")]
cif_files = [os.path.join(root, file) for root, _, files in os.walk(dir) for file in files if file.endswith("model_0.cif")]

for summary, full, ciffil in zip(summary_json_files, full_json_files, cif_files):
    new_obj = alphafold_obj(summary_jsonpath=summary, full_json_path=full, cifpath=ciffil, target_atoms=targets)

