#!/usr/bin/python3
#Changelog - SJB - 08/04/2024:
#Changed the output file name to allow for submission of these as jobs as they are quite long
#Previously did not work on the NR files - potentiall ymade a mix for this
##### parse HMMER tbl output and generate FASTA files from Uniprot and tab for comparisons 

import argparse
import gzip


def main():
    '''Captures arguments, parses hmmer file, generates fasta and deals with problems'''
# Deal with command line arguments
    parser = argparse.ArgumentParser(description='HMMER output parser')
    parser.add_argument('-f','--file', help='hmmer output file', required=True)
    parser.add_argument('-d','--db', help='e.g. uniprot_sprot.fasta.gz', required=True)
    parser.add_argument('-c','--cut',  help='score cut-off value (default: none)', default=0, type=float)    
    args = parser.parse_args()
    file_name = args.file.split("/")[-1]
# parse hmmer file and store ids in list
    acc_list = hmmer_parse(args.file,args.cut, file_name)
# generate fasta and catch any problematic ids
    missing = retrieve_local(args.db,acc_list, file_name)

# errors
    if len(missing)>0:
        print("The following entries were not found:\n") 
        for prob in missing:
            print(prob)

    
 
def hmmer_parse(infile, cutoff, output_name):
    '''Parse hmmer output, return list of ids, and write id /t score out'''
    hmmer_file = open(infile, 'r')
    score_out = open("new_cutoff_parse/" + output_name + "_scores.tsv", 'w')
    score_out.write("id\t" + "score\n")
    ids = []
    for line in hmmer_file:
        if line[0] != '#':
            NonNRflag = 0
            fields = line.split()
            acc = fields[0].split('|')
	    #If the length of acc after splitting is greater than 1, file is not a non redundant db
            if len(acc) > 1:
            	NonNRflag = 1
            score = float(fields[5])
            if (score > cutoff):
                if NonNRflag:
		    #If non redundant file is true, get the first entry (the actual id)
                    ids.append(acc[NonNRflag])
                    score_out.write(acc[NonNRflag] + "\t" + str(score) +"\n")
                else:
                    #If non redundant file, get the zeroth entry (the proper id)
                    ids.append(acc[NonNRflag])
                    score_out.write(acc[NonNRflag] + "\t" + str(score) +"\n")
                    
    hmmer_file.close()
    score_out.close()
    return(ids)
 


def retrieve_local(infile,hits, output_name):    
    '''retrieve sequences from local instance of uniprot'''
    acc = set(hits)
    with open("new_cutoff_parse/" + output_name + "_hits.fasta",'w') as fasta_file:
        up_file = gzip.open(infile, "rt")
        for line in up_file:
            if line.startswith('>'): 
                seq_flag = 0
                #There's probably a better way to do this with an additional function call to determine
                #If the file is a NR db or not but for now
                #Get a small sample of the line and split on |, if successful then header = the id (non-NR DB's have the |)
                #If no | found, header = the zeroth entry
                if len(line[0:4].split("|")) > 1:
                    header = line.split('|')[1]
                else:
                    header = line.strip(">").split(" ")[0]
                if header in acc:
                    fasta_file.write(line)
                    seq_flag = 1
                    acc.remove(header)
            else:
                if seq_flag == 1:
                    fasta_file.write(line)

        up_file.close()
    return(acc) 


# run things
if __name__== "__main__":
    main()

