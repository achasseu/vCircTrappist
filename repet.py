#Test Programme comparaison de kmer
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"f:o:")
ref=None

for opt, arg in optlist:
    if opt in "-f":
        ref = open(arg)
    if opt in "-o":
        count = open (str(arg)+"/aln_new_count.csv")
        new_new_count = open(str(arg)+"/aln_kmer_compare.csv",'w')

if ref == None:
    print("You forgot your fasta".upper())


fastref=ref.readlines()
count_lines = count.readlines()

ref=""
for line in fastref:
    if not ">" in line:
        ref=ref+line.replace("\n","")

# Function to generate k-mers from a sequence
def kmerizator (seq,taille):
    tot=taille
    kmer_list=[]
    while tot <= len (seq):
        kmer=seq[tot-taille:tot]
        tot=tot+1
        kmer_list.append(kmer)
    return (kmer_list)

# Function to compare k-mers and return the number of common kmers
def compare_kmer_score(seqa,seqb,taille):
    score = 0
    lista=kmerizator(seqa, taille)
    listb=kmerizator(seqb, taille)
    for kmer in lista:
        if kmer in listb:
            score = score + 1
    return score

# Function to find the length of the longest common k-mer
def plus_long_kmer_commun(seqa,seqb):
    length=len(seqa)
    while (length > 0):
        if compare_kmer_score(seqa, seqb, length) == 0:
            length=length-1
        else :
            score = length
            length = 0
    return score

# Function to compare k-mers and return a list of common k-mers
def compare_kmer_kmer(seqa,seqb,taille):
    score = 0
    lista=kmerizator(seqa, taille)
    listb=kmerizator(seqb, taille)
    kmer_list=[]
    for kmer in lista:
        if (kmer in listb) and (not kmer in kmer_list):
            score = score + 1
            kmer_list.append(kmer)
    return kmer_list
 
#find the common k-mers between the donor and acceptor splice sites for each line
for line in count_lines:
    line=line.split("\t")
    acceptor = line[1].upper()
    donor = line[2].upper()
    score = plus_long_kmer_commun(acceptor, donor)
    kmer_list=compare_kmer_kmer(acceptor,donor,score)
    line[4]=line[4].replace("\n","")
    if (score <= 6) and (ref.index(line[1])<ref.index(line[2])):
        new_new_count.write(line[0]+"\t"+acceptor+"\t"+ donor+"\t"+line[3]+"\t"+str(score)+"\t"+str(kmer_list)+"\t"+str(line[4:])+"\n")
        

