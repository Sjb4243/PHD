import re
trembl = "/home/sjb176/miniproject/results/fastas/smoltrembl.txt"
ids = "/home/sjb176/miniproject/results/fastas/ids.txt"

with open(trembl, "r") as trembl_file:
    #Dictionary comprehension to generate species IDs:species name using regex
    trembl_ids = {line.split("|")[1]:re.findall("(?<=OS=)\w* \w*", line.split("|")[2]) for line in trembl_file if line.startswith(">")}

with open(ids, "r") as id_file:
    checked_ids = [line.strip() for line in id_file]

for id in checked_ids:
    if id in trembl_ids.keys():
        print(trembl_ids[id])