#Filtering reads that have at least one primary and one secondary alignment
#Eliminating reads that map multiple times on the genome (repetitions ?)

#Importing the dependencies
import sys
import getopt

#Accepting arguments to make the program run properly
argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"a:o:")
fichier=None


for opt, arg in optlist:
    if opt in "-a":
        fichier = open(arg)
    if opt in "-o":
        out_file=open(str(arg)+"/aln_split.sam","w")

if fichier == None:
    print("You forgot your alignment file".upper())


#Initiating the search for splitted reads
lignes=fichier.readlines()
index = -1
line_count=-2
for line in lignes:
    line_count=line_count+1
for line in lignes:
    antiscore=0
    score = 0
    if index < line_count:
        index=index+1
    line2=line.split("\t")
    if "@" in line2[0]:
        out_file.write(line)
    else:
        if "0" in line2[1] and (not "AAAAAAAAAAAA" in line2[9]): #if a read is a primary alignment and do not contain a polyA stretch
            for line4 in lignes[index:index+30]: #look at the next 30 alignments
                line3=line4.split("\t")
                if (line3[0] == line2[0]) and (("256" in line3[1]) or (line3[2]!=line2[2])): #if a read is a supplementary (repetition) alignment, give it a bad score
                    antiscore = antiscore + 5
                if (line3[0] == line2[0]) and ("2048" in line3[1]) and (line3[2]==line2[2]): #if a read is a secondary alignment, give it a good score
                    score=score+1
        elif "2048" in line2[1] and (not "AAAAAAAAAAAA" in line2[9]): #if a read is a secondary alignment and do not contain a polyA stretch
            for line4 in lignes[index-10:index+10]: #look at the surrounding 20 alignments
                line4=line4.split("\t")
                if line4[0] == line2[0] and (("256" in line4[1]) or (line4[2]!=line2[2])): #if a read is a supplementary (repetition) alignment, give it a bad score
                    antiscore=antiscore+5
            if antiscore < 3:
                out_file.write(line) #if a read was not repeated, it is accepted


        #Same as above but with reverse complement
        if "16" in line2[1] and (not "AAAAAAAAAAAA" in line2[9]):
            for line4 in lignes[index:index+30]:
                line3=line4.split("\t")
                if (line3[0] == line2[0]) and (("272" in line3[1]) or (line3[2]!=line2[2])):
                    antiscore = antiscore + 5
                if (line3[0] == line2[0]) and ("2064" in line3[1]) and (line3[2]==line2[2]):
                    score=score+1
        if score > 0 and antiscore < 3:
            out_file.write(line)
    
        elif "2064" in line2[1] and (not ("AAAAAAAAAAAA" or "GGGGGGGGGGGG" or "CCCCCCCCCCCC" or "TTTTTTTTTTTT") in line2[9]):
            for line4 in lignes[index-10:index+10]:
                line4=line4.split("\t")
                if line4[0] == line2[0] and (("272" in line4[1]) or (line4[2]!=line2[2])):
                    antiscore=antiscore+5
            if antiscore < 3:
                out_file.write(line)
