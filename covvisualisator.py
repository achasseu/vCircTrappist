import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"b:e:o:")
for opt, arg in optlist:
    if opt in "-o":
        ssU2cov=open(str(arg)+"/circ_sense_U2_coverage.csv")
        ssnU2cov=open(str(arg)+"/circ_sense_nonU2_coverage.csv")
        asU2cov=open(str(arg)+"/circ_antisense_U2_coverage.csv")
        asnU2cov=open(str(arg)+"/circ_antisense_nonU2_coverage.csv")
        totcov=open(str(arg)+"/aln_virus_coverage.csv")

        ###sort the top 20 of the backsplice sites that are the most covered by reads
        sites_sorting=pd.read_csv(str(arg)+"/sites_sorting.csv",sep="\t",header=0)
sites_sorting.sort_values("Count", axis=0, ascending=False,inplace=True, na_position='first')
first_circ=sites_sorting[:20]
first_circ=first_circ.loc[:,["Count","5'_Position_(Sense_Strand)","3'_Position_(Sense_Strand)","Chromosome","Sense_Score","U2_Score"]]
first_circ=first_circ.values.tolist()


ssU2cov=ssU2cov.readlines()
ssnU2cov=ssnU2cov.readlines()
asU2cov=asU2cov.readlines()
asnU2cov=asnU2cov.readlines()
totcov=totcov.readlines()

#determining the mean coverage on the viral genome of the dataset --> normalization
total=0
for line in totcov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]

#creating the coverage table for the genereation of the coverage plot
ssU2cova=[]
ssnU2cova=[]
asU2cova=[]
asnU2cova=[]
for line in ssU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*1000000000
    ssU2cova.append(line)
for line in ssnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*1000000000
    ssnU2cova.append(line)
for line in asU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*-1000000000
    asU2cova.append(line)
for line in asnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*-1000000000
    asnU2cova.append(line)

#Generating the lists of the genomic segments that are processed
id_list=[]
for line in (ssU2cova+ssnU2cova+asU2cova+asnU2cova):
    if not (line[0] in id_list):
        id_list.append(line[0])

#function to draw the backsplice junctions up from the coverage plot
def draw_sashimi (start,end,height,sense,U2):
    xpoints=[]
    ypoints=[]
    if sense > 0:
        height=height+height/50
    else:
        height=height+height/50
    for i in range(end-start):
        xpoints.append(start+i)
        if sense > 0:
            if (i-(end-start)/2) < 0:
                height=height-((i-(end-start)/2)/((end-start)))
            elif ((i-(end-start)/2) > 0):
                height=height-((i-(end-start)/2)/((end-start)))
            else:
                height=height   
                     
        else:
            if (i-(end-start)/2) < 0:
                height=height+((i-(end-start)/2)/((end-start)))
            elif ((i-(end-start)/2) > 0):
                height=height+((i-(end-start)/2)/((end-start)))
            else:
                height=height
            
        ypoints.append(height)
    if float(U2)>0.6:
        plt.plot(xpoints,ypoints,color='black',linewidth=0.8)
    else:
        plt.plot(xpoints,ypoints,color='darkorange',linewidth=0.8)

#drawing the coverage plot with associated sashimis for each viral genomic segment
chrnb=len(id_list)
for i in range(chrnb):
    a=0
    b=0
    plt.figure(figsize=(12,4))
    beg=0
    xpoints=[]
    ypoints=[]

    seqU2ss=[]
    seqnU2ss=[]
    seqU2as=[]
    seqnU2as=[]
    for line in ssU2cova:
        if line[0]==id_list[i]:
            seqU2ss.append(line)
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='black',linewidth=0.5)
    xpoints=[]
    ypoints=[]
    for line in ssnU2cova:
        if line[0]==id_list[i]:
            seqnU2ss.append(line)
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='darkorange',linewidth=0.5)
    xpoints=[]
    ypoints=[]
    for line in asU2cova:
        if line[0]==id_list[i]:
            seqU2as.append(line)
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='black',linewidth=0.5)
    xpoints=[]
    ypoints=[]
    for line in asnU2cova:
        if line[0]==id_list[i]:
            seqnU2as.append(line)
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))
    plt.plot(xpoints, ypoints,color='darkorange',linewidth=0.5)
    en=int(totcov[-1].split("\t")[1])
    for line in first_circ:
        
        if line[3]==id_list[i]:
            if line[4]>0:
                a=0
                if line[1]>line[2]:
                    for seq in seqnU2ss[int(line[2]-(len(seqnU2ss))):int(line[1]+(len(seqnU2ss)))]:
                        if a<int(seq[2]):
                            a=int(seq[2])
                    for seq in seqU2ss[int(line[2]-(len(seqnU2ss))):int(line[1]+(len(seqnU2ss)))]:
                        if a<int(seq[2]):
                            a=int(seq[2])                       

                    draw_sashimi(line[2],line[1],(a),line[4],line[5])

                else :
                    for seq in seqnU2ss:
                        if a<int(seq[2]):
                            a=int(seq[2])
                    for seq in seqU2ss:
                        if a<int(seq[2]):
                            a=int(seq[2])   

                    draw_sashimi(line[1],line[2],(a),line[4],line[5])

            else:
                b=0
                if line[1]>line[2]:
                    
                    for seq in seqnU2as:
                        
                        if b>int(seq[2]):
                            b=int(seq[2])
                            
                    for seq in seqU2as:
                        
                        if b>int(seq[2]):
                            b=int(seq[2])

                    draw_sashimi(line[2],line[1],(b),line[4],line[5])

                else :
                    for seq in seqnU2as:
                        
                        if b>int(seq[2]):
                            b=int(seq[2])
                    for seq in seqU2as:
                        if b>int(seq[2]):
                            b=int(seq[2])

                    draw_sashimi(line[1],line[2],(b),line[4],line[5])
    for opt, arg in optlist:
        if opt in "-b":
            beg=int(arg)
        if opt in "-e":
            en =int(arg)
        plt.xlim(beg,en)


    
    plt.xlabel("Position on the genome")
    plt.ylabel("Coverage (RPB)")
    plt.hlines(0,beg,en,colors="black",label="Genome",linewidth=0.8)
    if a!=0:
        plt.hlines(a,beg,en,colors="lightgrey",label="Genome",linewidth=0.3)
    if b!=0:
        plt.hlines(b,beg,en,colors="lightgrey",label="Genome",linewidth=0.3)


    for opt, arg in optlist:
        if opt in "-o":
            plt.savefig(str(arg)+"/Circ_Coverage_{name}.png".format(name=str(id_list[i]).replace(".","_").replace("/","_")),dpi=1000,bbox_inches='tight')
    plt.clf()


