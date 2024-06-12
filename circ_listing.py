#Programme pour lister les circ
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"o:")

for opt, arg in optlist:
    if opt in "-o":
        sequencing=open(str(arg)+"/aln_split.sam")
        sequences=sequencing.readlines()
        backsplice_junctions=open(str(arg)+"/aln_bsj_out.txt")
        bsj_list=backsplice_junctions.readlines()
        backsplice_count=open(str(arg)+"/aln_kmer_compare.csv")
        count_list=backsplice_count.readlines()
        sam=open(str(arg)+"/aln_circ_list.sam",'w')

from difflib import *

circ_list=[]

# Function to generate k-mers from a sequence
def kmerizator (seq,taille):
    tot=taille
    kmer_list=[]
    while tot <= len (seq):
        kmer=seq[tot-taille:tot]
        tot=tot+1
        kmer_list.append(kmer)
    return (kmer_list)


circ_list=[]

# Extract circular RNA IDs from the count file
for line in count_list:
    line=line.split("\t")
    line[3]=int(line[3])

    if line[3] > 1:
        line[6]=line[6].replace("[","")
        line[6]=line[6].replace(",","\t")  
        line[6]=line[6].replace(">","")
        line[6]=line[6].replace("'","")
        line[6]=line[6].replace("]","")
        line[6]=line[6].replace("\"","")
        
        if " " in line[6]:
            line[6]=line[6].replace(" ","")
        if "\n" in line[6]:
            line[6]=line[6].replace("\n","")
        line[6]=line[6].split("\t")
        for id in line[6]:
            circ_list.append(id)

appa_list=[]

# Process each read in the sequencing file and write the reads in a new alignment file
for read in sequences:
    if ((sequences.index(read)%5000) == 0) and (sequences.index(read)!=0) :
        print(str(round((sequences.index(read)/len(sequences))*100),2)+"% "+"of the sequences are stocked in a new alignment file")
    reada=read.split("\t")
    if reada[0] in circ_list and (not read in appa_list):
        appa_list.append(read)
        sam.write(read)
print("Done")