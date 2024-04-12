fastafile = "/home/sjb176/miniproject/results/fastas/smoltrembl.txt"
counter = 0
with open(fastafile, "r") as fasta:
    file = fasta.readlines()

skipped = 0
full_fastas = []
for i in range(len(file)):
    index = i + skipped
    if index >= len(file):
        break
    if file[index].startswith(">"):
        done = False
        newPos = index + 1
        fasta_to_append = []
        fasta_to_append.append(file[index])
        runs = 0
        while not done:
            if newPos >= len(file):
                break
            if file[newPos].startswith(">"):
                full_fastas.append(fasta_to_append)
                done = True
                skipped = skipped + runs
            if not done:
                fasta_to_append.append(file[newPos])
            newPos += 1
            runs += 1
fastafile = "/home/sjb176/miniproject/results/fastas/smoltrembl.txt"

with open(fastafile, "r") as fasta:
    file = fasta.readlines()
skipped = 0
full_fastas = []
for i in range(len(file)):
    index = i + skipped
    if index >= len(file):
        break
    if file[index].startswith(">"):
        done = False
        newPos = index + 1
        fasta_to_append = []
        fasta_to_append.append(file[index].strip())
        runs = 0
        while not done:
            if newPos + 1 > len(file):
                full_fastas.append(fasta_to_append)
                break
            if file[newPos].startswith(">"):
                full_fastas.append(fasta_to_append)
                skipped = skipped + runs
                done = True
            if not done:
                fasta_to_append.append(file[newPos].strip())
                newPos += 1
                runs += 1

long_fastas = []
print(len(full_fastas))
for i in full_fastas:
    id = i[0].split("|")[1]
    fasta_string = "".join(i[1:len(i)])
    if len(fasta_string) > 200 and len(fasta_string) < 800:
        long_fastas.append(id + "\n")
    else:
        continue

with open("long_fastas.txt", "w") as outfile:
    for fasta in long_fastas:
        outfile.write(fasta)
