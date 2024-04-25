#To be used to generate a metadata file for MSA input files
from helper_functions import *
p2xdictohmm = "/home/sjb179/miniproject/results/fastas/P2Xdicto.hmm_uniprot_trembl.fasta.gz.txt_hits.fasta"
p2xstandardhmm = "/home/sjb179/miniproject/results/fastas/standard_P2X_HMM.txt_uniprot_trembl.fasta.gz.txt_hits.fasta"

p2xdictofastas = get_fastas(p2xdictohmm)
p2xstandardhmmfastas = get_fastas(p2xstandardhmm)


ids_from_dicto = [line[0] for line in p2xdictofastas]
ids_from_standard = [line[0] for line in p2xstandardhmmfastas]

common = [value for value in ids_from_dicto if value in ids_from_standard]
not_common = [value for value in ids_from_dicto if value not in ids_from_standard]
not_common2 = [value for value in ids_from_standard if value not in ids_from_dicto]

print(len(not_common))
print(len(not_common2))

print(len(ids_from_dicto))
print(len(ids_from_standard))