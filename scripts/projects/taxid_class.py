from Bio import Entrez
import re
import time
import sys
class ncbi_object():
    existing_ids = set()
    def __init__(self, id_line):
        self.id_line = id_line
        self.taxid = re.findall("OX[=|_][0-9]*", id_line)[0]
        self.parents = ""
        self.scientific_rank = {}

    def fetch(self):
        reset_count = 0
        if self.taxid not in ncbi_object.existing_ids:
            try:
                handle = Entrez.efetch(db="taxonomy", id=self.taxid, retmode="xml")
                records = Entrez.read(handle)
                records = records[0]
                parents = records["LineageEx"]
                self.parents = parents
                ncbi_object.existing_ids.add(self.taxid)
                print(f"Fetched {self.taxid}")
                self.extract()
            except Exception as e:
                if reset_count == 5:
                    print(f"Retried too many times. NCBI could be down or you're not using an API key.")
                    sys.exit(1)
                time.sleep(2)
                print("Connection failure: retrying...")
                reset_count +=1
                self.fetch()

    def extract(self):
        clade_counter = 1
        name_rank = {}
        for dict in self.parents:
            rank = dict["Rank"]
            if rank != "no rank":
                name = dict["ScientificName"]
                if rank == "clade":
                    rank = rank + str(clade_counter)
                    clade_counter += 1
                name_rank[rank] = name
        self.scientific_rank[self.id_line] = name_rank

