import re
import time
import sys
import json
class fasta_obj:
    def __init__(self, fasta_dir):
        self.base_fasta = self.parse(fasta_dir)
        self.lines = [line[0] for line in self.base_fasta]
        self.base_ids = self.get_ids(self.base_fasta)
        self.taxids = self.get_taxids(self.base_fasta)
        self.species = self.get_species(self.base_fasta)
        self.proteins = self.get_proteins(self.base_fasta)
        self.families = []
        self.__reset_count = 0
    #Wrapper for printing out the return of functions, probably obsolete now
    def print_return(func):
        def wrapper(*args, **kwargs):
            print_result = kwargs.pop('print_res', False)
            result = func(*args, **kwargs)
            if print_result:
                print(f"return_val: {result}")
            return result
        return wrapper
    #Parses fasta files, if the object is a list (aka its a fasta_list provided from somewhere else)
    #It just returns the list
    def parse(self, obj):
        if isinstance(obj, list):
            return obj
        with open(obj, "r") as infile:
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
    #If the "fasta" is a single string, aka just the header just split that
    #Otherwise do the same thing for the whole fasta file
    #This is the same for the next few functions as well
    def get_ids(self, fasta):
        if isinstance(fasta, str):
            if "|" in fasta:
                ids = fasta.split("|")[1]
            else:
                ids = fasta.split(" ")[0].strip(">")
        else:
            for line in fasta:
                print(line)
            ids = [line[0].split("|")[1] if "|" in line[0] else line[0].split(" ")[0].strip(">") for line in fasta]
        return ids

    def get_taxids(self, fasta):
        if isinstance(fasta, str):
            taxids = re.findall(r"(?<=OX[=_])[0-9]*", fasta)[0]
        else:
            taxids = [re.findall(r"(?<=OX[=_])[0-9]*", id_line[0])[0] for id_line in fasta]
        return taxids

    def get_species(self, fasta):
        species = [
            re.findall(r"(?<=OS=).*(?=[ |_]OX)|\[[a-z A-Z0-9\.]*\]", id_line[0])[0] if re.findall(r"(?<=OS=).*(?=[ |_]OX)|\[[a-z A-Z0-9\.]*\]", id_line[0])
            else "UNKNOWN"
            for id_line in fasta
            ]
        return species

    def get_proteins(self, fasta):
        proteins_list = []
        for id_line in fasta:
            protein = re.findall(r"(?<=_[A-Z0-9]{5}[_| ]).*(?=[\||_| ]OS)|(?<=\|)\w*(?=\|OS)", id_line[0])
            protein = protein[0] if protein else "UNKNOWN"
            proteins_list.append(protein)
        return proteins_list
    #Return sequences where the length is between min and max
    def get_length(self, min, max):
        return_fasta = [
            item for item in self.base_fasta
            if int(min) <= len(item[1]) <= int(max)
        ]
        return fasta_obj(return_fasta)

    def fetch_family(self, secrets, save = "", read = ""):
        #Set up credentials found in file
        if not read:
            self.__set_creds(secrets)
            #for each full line, start the fetching process
            for header in self.lines:
                self._start_fetch(header)
            if save:
                with open(save, "w") as out:
                    json.dump(self.families, out)
        else:
            with open(read, "r") as infile:
                self.families = json.load(infile)
        return self

    def __set_creds(self, secrets):
        from Bio import Entrez
        with open(secrets) as creds:
            file = creds.readlines()
            Entrez.email = file[0].split("=")[-1].strip()
            Entrez.api_key = file[1].split("=")[-1].strip()
    ##begin fetching
    def _start_fetch(self, header):
        from Bio import Entrez
        taxid = self.get_taxids(header)
        try:
            #standard entrez stuff
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            records = records[0]
            parents = records["LineageEx"]
            print(f"Fetched {taxid}")
            self.__reset_count = 0
            clade_counter = 1
            name_rank = {}
            end_dict = {}
            for dict in parents:
                rank = dict["Rank"]
                if rank != "no rank":
                    name = dict["ScientificName"]
                    if rank == "clade":
                        rank = rank + str(clade_counter)
                        clade_counter += 1
                    name_rank[rank] = name
            end_dict[header] = name_rank
            self.families.append(end_dict)
            #Basically always failed connections
        except Exception as e:
            print(e)
            if self.__reset_count == 5:
                print(f"Retried too many times. NCBI could be down or you're not using an API key.")
                sys.exit(1)
            time.sleep(2)
            print("Connection failure: retrying...")
            self.__reset_count += 1
            self._start_fetch(header)


    #Take the user's ID list
    #If anti = false, it will return a fasta where the IDs between the two match
    #If anti = true, it will return a fasta where the IDs are in self but NOT in the ID list (or other fasta object)
    def _compare_ids(self, id_list, anti):
        compare_dict = {self.get_ids(line[0]): line for line in self.base_fasta}
        return_fasta = [value for key, value in compare_dict.items()
                        if (not anti and key in id_list) or (anti and key not in id_list)]
        return return_fasta

    # Take the user's taxid list If anti = false, it will return a fasta where the taxids between the two match If
    # anti = true, it will return a fasta where the taxids are in self but NOT in the taxid list (or other fasta
    # object)
    def _compare_taxids(self, taxid_list, anti):
        compare_dict = {tuple(line):self.get_taxids(line[0]) for line in self.base_fasta}
        return_fasta = [key for key, value in compare_dict.items()
                        if (not anti and value in taxid_list) or (anti and value not in taxid_list)]
        return return_fasta

    #House keeping - Deciding whether the ID list is a list, or if its an object
    #If its a list, it can go straight to the functions
    #If its not a list, the corresponding desired list
    def compare(self,list_of, by = "id", anti=False):
        if not isinstance(list_of, list):
            if by == "id":
                return_fasta = self._compare_ids(list_of.base_ids, anti)
            if by == "taxid":
                return_fasta = self._compare_taxids(list_of.taxids, anti)
        else:
            if by == "id":
                return_fasta = self._compare_ids(list_of,anti)
            if by == "taxid":
                return_fasta = self._compare_taxids(list_of,anti)

        return fasta_obj(return_fasta)

    def get_by_family(self, search):
        return_fasta = []
        for fast, dict in zip(self.base_fasta, self.families):
            if fast[0] == list(dict.keys())[0]:
                #Go to inner dict
                inner = dict[list(dict.keys())[0]]
                values = [val for val in inner.values()]
                if search in values:
                    return_fasta.append(fast)
                print(values)
            else:
                print(f"Something has gone horribly wrong")
        return fasta_obj(return_fasta)

    def write_out(self, dir):
        with open(dir, "w") as out:
            for line in self.base_fasta:
                print(line)
                out.write(line[0] + "\n")
                out.write(line[1] + "\n")

