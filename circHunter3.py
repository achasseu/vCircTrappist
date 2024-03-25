#We look for backsplicing signatures
import sys
import getopt
from cigar import Cigar

#accepting the reference file as argument
argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"f:o:")
ref=None

for opt, arg in optlist:
    if opt in "-f":
        ref = open(arg)
    if opt in "-o":
        out_file=open(str(arg)+"/aln_circ.sam",'w')                             #Creating a new sequence file containing the circRNA reads
        out_bsj=open(str(arg)+"/aln_bsj_out.txt",'w')                           #Creating a list containing the backsplice junctions
        bsj_sites=open(str(arg)+"/aln_bsj_sites.txt",'w')                       #Creating a list containing the backsplice sites
        verifa=open(str(arg)+"/aln_verifa.txt",'w')                             #Opening the alignment file containing splitted reads
        fichier = open(str(arg)+"/aln_split.sam")                               #creating a text file to crosscheck the reads

if ref == None:
    print("You forgot your fasta".upper())

                               
                          

lignes = fichier.readlines()                                    
#l_fastq=fastq.readlines()                                       

refa=ref.readlines()


index = -1                                                      #Monitoring the indexes of the lines to target them easily
count = -2                                                      #Counting the number of lines

circ_list=[]                                                    #Opening a list containing the circ reads IDs
bsj_list=[]

for line in lignes:
    count = count+1                                             #Counting lines
for line in lignes:                                             #Entering the loop to find circRNAs
    if index < count:                                           #While the index is smaller than the count, increment
        index=index+1
    if "@" in line[0]:                                             #We keep the header line of the SAM file
        out_file.write(line)
    else:                                         #For all the lines of the SAM files
        line2 =line.split("\t")                                 #Each variable of a line is saved in a list
        for line4 in lignes[index+1:index+30]:
            line3 = line4.split("\t")                               #Same with the next lines, to compare them
            if line2[0] == line3[0]:                                #If the read ID is the same as the next line
                line3[3]=int(line3[3])
                line2[3]=int(line2[3])
                if line2[3]>line3[3] and (line2[3]-line3[3]<50000) and (line2[3]-(line3[3]+len(line3[9]))+15>0):      #If the read is less than 50 000 bases downstream from the next one
                    id_circ=line2[0]                               #We save the read ID
                    length=len(line3[9])+15
                    if line3[9] in line2[9][-length:]:              #If the primary alignment is a circRNA, we look for the next alignment
                        if not (id_circ in circ_list):              #If we don't have saved the ID yet (if it's a new identification), we save it
                            circ_list.append(id_circ)               
                            out_file.write(line)

                if line2[3]<line3[3] and (line3[3]-line2[3]<50000) and (line3[3]-(line2[3]+len(line2[9])-len(line3[9]))+15>0):       #If the read is less than 50 000 bases upstream from the next one
                    id_circ=line2[0]
                    length=len(line3[9])+15
                    if line3[9] in line2[9][:length]:               #If the primary alignment is a circRNA, we look for the secondary alignment at the end of the primary one
                        if not (id_circ in circ_list):              #If we don't have saved the ID yet (if it's a new identification), we save it
                            circ_list.append(id_circ)               
                            out_file.write(line)
def repet(seqa,seqb): #program to find repeated kmers by comparing two mated sequences
    kmera=[]
    kmerb=[]
    total = 0
    for a in range(len(seqa)-6):
        kmera.append(seqa[a:a+7])
    for b in range(len(seqb)-6):
        kmerb.append(seqb[b:b+7])
    for c in kmera:
        if c in kmerb:
            total=total+1
    return (total)

def plus_longue_fin_commune(seqa,seqb): #program to find the longest common end between two strings
    total = len(seqa)
    subseq = seqa[-total:]
    while (not subseq in seqb[:total+1]) and total > 0:
        subseq = seqa[-total:]
        total=total-1
    return subseq

index = -1                                                      
for read in lignes:
    if index < count:                                           
        index=index+1
        read=read.split("\t")
    if (read[0] in circ_list) and not (read[0] in bsj_list) : #If a read is a circ and not yet saved, start analyzing it       
        bsj_list.append(read[0])
        listereads = []
        read[3]=int(read[3])
        ma_sequence=read[9]
        for read2 in lignes[index+1:index+30]: #trim the reads to keep only the matching sequence
            read2=read2.split("\t")
            read2[3]=int(read2[3])
            if read2[0] == read[0]:
                listereads.append(read2)
                if read2[3] > read[3] and (read2[9] in read[9][:((len(read[9])//2)+15)]) and (read2[3]-read[3]<100000):
                    loc=(read[9].index(read2[9]))+len(read2[9])
                    read[9]=read[9][loc:]
                if read2[3] > read[3] and (not read2[9] in read[9][:((len(read[9])//2)+15)]) and (read2[3]-read[3]<100000) and (len(plus_longue_fin_commune(read2[9],read[9][:((len(read[9])//2)+15)])) > 15):
                    read2[9]=plus_longue_fin_commune(read2[9],read[9])
                    loc=(read[9].index(read2[9]))+len(read2[9])
                    read[9]=read[9][loc:]

                if (read2[3] < read[3]) and (read2[9] in read[9][-((len(read[9])//2)+15):]) and (read[3]-read2[3]<100000):
                    loc=read[9].index(read2[9])
                    read[9]=read[9][:loc]      
                elif (read2[3] < read[3]) and (not read2[9] in read[9][-((len(read[9])//2)+15):]) and (read[3]-read2[3]<100000) and (len(plus_longue_fin_commune(read[9],read2[9][-((len(read[9])//2)+15):])) > 15):
                        read2[9]=plus_longue_fin_commune(read[9],read2[9])
                        loc=read[9].index(read2[9])
                        read[9]=read[9][:loc] 
                                 
                read[9]=read[9].replace(read2[9],"")
        listereads.append(read)
        listereads=sorted(listereads,key=lambda x: x[3])#sort the reads according to their matching positions on the genome
        locb=listereads[0][3]
        c=Cigar(listereads[-1][5])
        
        
        e =list(c.items())
        d=len(c)
        
        for tortellini in e: #locate the positions of the splicing on the reads
            
            if tortellini[1] == 'S':
                d=d-tortellini[0]
            
        if d == 150:
            print (listereads)
        loca=listereads[-1][3]+d
        
        

        refb=str()
        for chrom in refa: #save the backsplice junctions and backsplice sites from the reads, using their locations
            if listereads[0][2] in chrom:
                cbsj=refa[refa.index(chrom)+1][loca-16:loca+14]
                dbsj=refa[refa.index(chrom)+1][locb-15:locb+15]
            if not ">" in chrom:
                refb=refb+chrom.replace("\n","")
                
        
        bbsj=listereads[-1][9][-15:]
        absj=listereads[0][9][:15]

        if (repet(absj, bbsj) == 0) and (not "TTAGGGTT" in (absj or bbsj or cbsj or dbsj)): #if there are no telomeric repeats nor repetitions in the backsplice junctions, start saving the sequences
            bsj=bbsj+absj
            if not bsj in ma_sequence:
                if listereads[0][9] in ma_sequence[-(len(listereads[0][9])+15):]:
                    position_read=ma_sequence.index(listereads[0][9])
                    bsj = ma_sequence[position_read-15:position_read+15]
                elif listereads[-1][9] in ma_sequence[:(len(listereads[-1][9])+15)]:
                    position_read=ma_sequence.index(listereads[-1][9])+len(listereads[-1][9])
                    bsj = ma_sequence[position_read-15:position_read+15]
                else:
                    bsj=bbsj+absj
            if (not bsj in refb) and (not "\n" in bsj) and (not "\n" in dbsj) and (not "\n" in cbsj) and (len(cbsj+dbsj) > 40) :    
                out_bsj.write(">"+read[0]+"\t"+bsj+"\n")
                bsj_sites.write(">"+read[0]+"\t"+bsj+"\t"+dbsj+"\t"+cbsj+"\n")

