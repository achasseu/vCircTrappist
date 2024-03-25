###Program to characterize the backsplice sites
###It gives a score to determine if the sites are canonical or not
###It also extracts the strandness of the sites
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"f:g:s:o:")
fastafile =None
gfffile=None

for opt, arg in optlist:
    if opt in "-f":
        fastafile = open(arg)
    if opt in "-g":
        gfffile = open (arg)
    if opt in "-s":
        strand="{}".format(arg)
    if opt in "-o":
        sequencing=open(str(arg)+"/aln_circ_list.sam")
        finalsorting=open(str(arg)+"/sites_sorting.csv",'w')
        alnpluscan = open(str(arg)+"/circ_sense_U2_a.sam","w")
        alnminuscan = open(str(arg)+"/circ_antisense_U2_a.sam","w")
        alnplusnocan = open(str(arg)+"/circ_sense_nonU2_a.sam","w")
        alnminusnocan = open(str(arg)+"/circ_antisense_nonU2_a.sam","w")
        backsplice_count=open(str(arg)+"/aln_kmer_compare.csv")

if fastafile == None or gfffile==None:
    print("You forgot your arguments".upper())


sequences=sequencing.readlines()

fasta=fastafile.readlines()

gfflines=gfffile.readlines()

count_list=backsplice_count.readlines()



finalsorting.write(str("Backjunction\t5'_Splice_Site_(Sense_Strand)\t3'_Splice_Site_(Sense_Strand)\tCount\tSense_Score\tU2_Score\t5'_Position_(Sense_Strand)\t3'_Position_(Sense_Strand)\tMax_CircRNA_Size\tColocalized_Gene\tLongest_Common_Kmer\tCommon_Sequences\tChromosome\tIdentifiers_Of_The_Reads\n"))

if strand=="R": #if the strandness of the library is reversed
    for line in count_list:
        lineb=line
        line=line.split("\t")
        line[3]=int(line[3])
        circ_list=[]


        if line[3] > 1: #if there is more than one count for the analyzed read
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
            porm=0
            for id in line[6]: 
                for read in sequences :
                    reada=read.split("\t")
                    if reada[0] == id:
                        circ_list.append(reada)
            for read in circ_list: #determine the sense of the mappings
                read[1]=int(read[1])
                if read[1] == 0 or read[1]==2048:
                    porm=porm-1
                    
                elif read[1] == 16 or read[1] == 2064 :
                    porm=porm+1
                
            score = porm/line[3]
            line=line[:4]+[score]+line[4:]
            
            acceptor=line[1]
            donor=line[2]
            if score > 0: #determine if the transcript was processed through canonical splicing by analyzing a sliding window surrounding the sites
                if "AG" in acceptor[13:15] :
                    scoreacceptor=0.5
                elif "AG" in acceptor[12:16]:
                    scoreacceptor=0.4
                elif "AG" in acceptor[11:17]:
                    scoreacceptor=0.3
                elif "AG" in acceptor[10:18]:
                    scoreacceptor=0.2
                elif "AG" in acceptor[9:19]:
                    scoreacceptor=0.1
                else :
                    scoreacceptor=0
                if "GT" in donor[15:17]:
                    scoredonor=0.5
                elif "GT" in donor[14:18]:
                    scoredonor=0.4
                elif "GT" in donor[13:19]:
                    scoredonor=0.3
                elif "GT" in donor[12:20]:
                    scoredonor=0.2
                elif "GT" in donor[11:21]:
                    scoredonor=0.1
                else :
                    scoredonor=0    
            else :
                if "AC" in acceptor[13:15] :
                    scoreacceptor=0.5
                elif "AC" in acceptor[12:16]:
                    scoreacceptor=0.4
                elif "AC" in acceptor[11:17]:
                    scoreacceptor=0.3
                elif "AC" in acceptor[10:18]:
                    scoreacceptor=0.2
                elif "AC" in acceptor[9:19]:
                    scoreacceptor=0.1
                else :
                    scoreacceptor=0
                if "CT" in donor[15:17]:
                    scoredonor=0.5
                elif "CT" in donor[14:18]:
                    scoredonor=0.4
                elif "CT" in donor[13:19]:
                    scoredonor=0.3
                elif "CT" in donor[12:20]:
                    scoredonor=0.2
                elif "CT" in donor[11:21]:
                    scoredonor=0.1
                else :
                    scoredonor=0 
            for read in sequences: #write the chromosome on which was mapped the read
                reada=read.split("\t")
                circ_list=line[-1]
                if reada[0] == circ_list[0] and (line[7] != reada[2]):
                    chromosome=reada[2]
                    line=line[:7]+[chromosome]+line[7:] 

            for chrom in fasta: #write the max potential size of the candidate circRNA
                if line [7] in chrom:
                    positiondonor=fasta[fasta.index(chrom)+1].upper().index(donor)+15
                    positionacceptor=fasta[fasta.index(chrom)+1].upper().index(acceptor)+15   
            size=positiondonor-positionacceptor
            
            scoreU2=scoreacceptor+scoredonor
            line=line[:5]+[scoreU2,positionacceptor,positiondonor,size]+line[5:] 


            appa_list=[]
            for read in sequences: #create lists of sequences according to their strandness and their splicing pattern
                reada=read.split("\t")
                circ_list=line[-1]
                
                if reada[0] in circ_list and (not read in appa_list):
                    appa_list.append(read)
                    if line[4] > 0 and line[5] > 0.6 :
                        alnpluscan.write(read)
                    if line[4] > 0 and line[5] < 0.6 :
                        alnplusnocan.write(read)                
                    if line[4] <= 0 and line[5] > 0.6 :
                        alnminuscan.write(read)
                    if line[4] <= 0 and line[5] < 0.6 :
                        alnminusnocan.write(read)    
            
            featurelist=[]
            if line[4] > 0:
                sense="+"
            else:
                sense="-"

            for feature in gfflines: #determine on which gene the circRNA is mapped
                if "#" not in feature:
                    feature=feature.split("\t")
                    identifying=feature[-1].split(";")
                    
                    if "gene" in feature[2]:
                        
                        if (int(feature[3]) <= line[6] <= int(feature[4]) or int(feature[3]) <= line[7] <= int(feature[4])) and feature[0] == line[11]:
                            if sense != feature[6]:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring="Anti-"+obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)
                            else:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring=obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)

                    elif "CDS" in feature[2] and featurelist==[]: 
                        if (int(feature[3]) <= line[6] <= int(feature[4]) or int(feature[3]) <= line[7] <= int(feature[4])) and feature[0] == line[11]:
                            if sense != feature[6]:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring="Anti-"+obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)
                            else:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring=obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)

            if featurelist==[]:
                featurelist=["Unknown"]
                
            line=line[:9]+[featurelist]+line[9:]



            scoop=str()
            for element in line:
                if int(line[8])>0:
                    if element == line[-1]:
                        scoop=scoop+str(element)+"\n"
                    else:
                        scoop=scoop+str(element)+"\t"
            finalsorting.write(scoop)

if strand=="F": #same but for Forward-stranded libraries
    for line in count_list:
        lineb=line
        line=line.split("\t")
        line[3]=int(line[3])
        circ_list=[]


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
            porm=0
            for id in line[6]:
                for read in sequences :
                    reada=read.split("\t")
                    if reada[0] == id:
                        circ_list.append(reada)
            for read in circ_list:
                read[1]=int(read[1])
                if read[1] == 0 or read[1]==2048:
                    porm=porm-1
                    
                elif read[1] == 16 or read[1] == 2064 :
                    porm=porm+1
                
            score = porm/line[3]
            line=line[:4]+[-score]+line[4:]
            
            acceptor=line[1]
            donor=line[2]
            if score > 0:
                if "AG" in acceptor[13:15] :
                    scoreacceptor=0.5
                elif "AG" in acceptor[12:16]:
                    scoreacceptor=0.4
                elif "AG" in acceptor[11:17]:
                    scoreacceptor=0.3
                elif "AG" in acceptor[10:18]:
                    scoreacceptor=0.2
                elif "AG" in acceptor[9:19]:
                    scoreacceptor=0.1
                else :
                    scoreacceptor=0
                if "GT" in donor[15:17]:
                    scoredonor=0.5
                elif "GT" in donor[14:18]:
                    scoredonor=0.4
                elif "GT" in donor[13:19]:
                    scoredonor=0.3
                elif "GT" in donor[12:20]:
                    scoredonor=0.2
                elif "GT" in donor[11:21]:
                    scoredonor=0.1
                else :
                    scoredonor=0    
            else :
                if "AC" in acceptor[13:15] :
                    scoreacceptor=0.5
                elif "AC" in acceptor[12:16]:
                    scoreacceptor=0.4
                elif "AC" in acceptor[11:17]:
                    scoreacceptor=0.3
                elif "AC" in acceptor[10:18]:
                    scoreacceptor=0.2
                elif "AC" in acceptor[9:19]:
                    scoreacceptor=0.1
                else :
                    scoreacceptor=0
                if "CT" in donor[15:17]:
                    scoredonor=0.5
                elif "CT" in donor[14:18]:
                    scoredonor=0.4
                elif "CT" in donor[13:19]:
                    scoredonor=0.3
                elif "CT" in donor[12:20]:
                    scoredonor=0.2
                elif "CT" in donor[11:21]:
                    scoredonor=0.1
                else :
                    scoredonor=0 
            for read in sequences:
                reada=read.split("\t")
                circ_list=line[-1]
                if reada[0] == circ_list[0] and (line[7] != reada[2]):
                    chromosome=reada[2]
                    line=line[:7]+[chromosome]+line[7:] 

            for chrom in fasta:
                if line [7] in chrom:
                    positiondonor=fasta[fasta.index(chrom)+1].upper().index(donor)+15
                    positionacceptor=fasta[fasta.index(chrom)+1].upper().index(acceptor)+15   
            size=positiondonor-positionacceptor
            
            scoreU2=scoreacceptor+scoredonor
            line=line[:5]+[scoreU2,positionacceptor,positiondonor,size]+line[5:] 


            appa_list=[]
            for read in sequences:
                reada=read.split("\t")
                circ_list=line[-1]
                
                if reada[0] in circ_list and (not read in appa_list):
                    appa_list.append(read)
                    if line[4] > 0 and line[5] > 0.6 :
                        alnpluscan.write(read)
                    if line[4] > 0 and line[5] < 0.6 :
                        alnplusnocan.write(read)                
                    if line[4] <= 0 and line[5] > 0.6 :
                        alnminuscan.write(read)
                    if line[4] <= 0 and line[5] < 0.6 :
                        alnminusnocan.write(read)    
            
            featurelist=[]
            if line[4] > 0:
                sense="+"
            else:
                sense="-"

            for feature in gfflines:
                if "#" not in feature:
                    feature=feature.split("\t")
                    identifying=feature[-1].split(";")
                    
                    if "gene" in feature[2]:
                        
                        if (int(feature[3]) <= line[6] <= int(feature[4]) or int(feature[3]) <= line[7] <= int(feature[4])) and feature[0] == line[11]:
                            if sense != feature[6]:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring="Anti-"+obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)
                            else:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring=obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)

                    elif "CDS" in feature[2] and featurelist==[]: 
                        if (int(feature[3]) <= line[6] <= int(feature[4]) or int(feature[3]) <= line[7] <= int(feature[4])) and feature[0] == line[11]:
                            if sense != feature[6]:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring="Anti-"+obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)
                            else:
                                for obj in identifying:
                                    if "Name=" in obj:
                                        featuring=obj[5:]
                                        if "\n" in featuring:
                                            featuring=featuring.replace("\n","")
                                        featurelist.append(featuring)

            if featurelist==[]:
                featurelist=["Unknown"]
                
            line=line[:9]+[featurelist]+line[9:]



            scoop=str()
            for element in line:
                if int(line[8])>0:
                    if element == line[-1]:
                        scoop=scoop+str(element)+"\n"
                    else:
                        scoop=scoop+str(element)+"\t"
            finalsorting.write(scoop)

