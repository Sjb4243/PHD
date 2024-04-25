from helper_functions import *

small_trembl = "/home/sjb176/miniproject/results/fastas/smoltrembl.fasta"
with open(small_trembl, "r") as fastas:
    file = fastas.readlines()
fasta_list = get_fastas(file)

with open("/home/sjb176/miniproject/merged_dfs/merged/new_trembl") as scorefile:
    file = scorefile.readlines()
    file.pop(0)
    good_ids = []
    for line in file:
        split_data = line.split("\t")
        id = split_data[0]
        score1 = split_data[1].strip()
        score2 = split_data[2].strip()
        if float(score1) < 100 and float(score2) < 100:
            good_ids.append(id)


with open("/home/sjb176/PhD/results/procced_for_blast.txt", "a") as output:
    for line in fasta_list:
        id = line[0]
        if id.split("|")[1] in good_ids:
            #test_fasta = "".join(line[1:len(line)])
            #test_fasta = test_fasta.replace("-", "")
            #if len(test_fasta) < 200 or len(test_fasta) > 800:
                #continue
            fasta_string = "".join(line[1:len(line)])
            #fasta_string = fasta_string[38:480]
            fasta_string = fasta_string.replace("-", "")
            output.write(id + "\n")
            output.write(fasta_string + "\n")






# amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
# with open("trembl_combined_no_bad_scores.txt", "a") as output:
#     for line in fasta_list:
#         id = line[0]
#         id_seg = id.split("|")[1]
#         if id_seg not in bad_ids:
#             amino = "".join(line[1:len(line)])
#             if "X" in amino:
#                 continue
#             dat = extract_data(id)
#             output.write(make_id(">" + dat[0], dat[1], dat[2]))
#             output.write(amino + "\n")
#
