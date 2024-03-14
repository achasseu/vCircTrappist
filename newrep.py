import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import (BiopythonTranslator,GraphicFeature,GraphicRecord)
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"g:b:e:")
for opt, arg in optlist:
    if opt in "-g":
        gfffile = open (arg)
totcov=open("./aln_virus_coverage.csv")
totcov=totcov.readlines()
sites_sorting=pd.read_csv("./sites_sorting.csv",sep="\t",header=0)
sites_sorting.sort_values("Count", axis=0, ascending=False,inplace=True, na_position='first')
first_circ=sites_sorting[:102]
first_circ=first_circ.loc[:,["Count","5'_Position_(Sense_Strand)","3'_Position_(Sense_Strand)","Chromosome","Sense_Score","U2_Score","Max_CircRNA_Size"]]
circs_list=first_circ.values.tolist()
totcov=open("./aln_virus_coverage.csv")
totcov=totcov.readlines()


gff=gfffile.readlines()
features=[]
for line in gff:
    if "#" not in line[0]:
        line=line.split("\t")
        if (int(line[4])-int(line[3]))>(int(totcov[-1].split("\t")[1])/150):
            if line[2]=="CDS" and line[6] == "-":
                features.append(GraphicFeature(start=int(line[3]),end=int(line[4]),strand=-1,color="black",label=None))
            if line[2]=="CDS" and line[6] == "+":
                features.append(GraphicFeature(start=int(line[3]),end=int(line[4]),strand=+1,color="black",label=None))


id_list=[]
for line in (circs_list):
    if not (line[3] in id_list):
        id_list.append(line[3])

record=GraphicRecord(sequence_length=int(totcov[-1].split("\t")[1]),features=features)
ax1, _ = record.plot(figure_width=20)
ax1.figure.savefig("./genomeview.png",dpi=1000,bbox_inches='tight')
plt.figure(figsize=(12,2))

def draw_sashimi (start,end,height,sense,U2,max):
    xpoints=[]
    ypoints=[]
    if sense > 0:
        height=height+height/50
    else:
        height=height+height/50
    c=max
    xa=(start-end)/2
    xb=(end-start)/2
    a=c/(-xa*xb)
    for i in range(end+1-start):
        xpoints.append(start+i)
        j=i
        j=j-((end-start)/2)

        if sense > 0:
            y=-a*(j-xa)*(j-xb)   
             
        else:
            y=a*(j-xa)*(j-xb)

              

        ypoints.append(y)
        
    if float(U2)>0.6:
        plt.plot(xpoints,ypoints,color='black',linewidth=0.8)
    else:
        plt.plot(xpoints,ypoints,color='darkorange',linewidth=0.8)

total=0
for line in totcov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]


chrnb=len(id_list)
for i in range(chrnb):
    covloc=[]
    for line in totcov:
        line=line.replace("\n","")
        line=line.split("\t")
        if line[0] == id_list[i]:
            covloc.append(line)
    beg=0
    for line in circs_list:
        if line[3]==id_list[i]:
            draw_sashimi(int(line[1]),int(line[2]),0,int(line[4]),line[5],(((line[0])/total)*1000000000))
    en=int(covloc[-1][1])
    for opt, arg in optlist:
        if opt in "-b":
            beg=int(arg)
        if opt in "-e":
            en =int(arg)
        plt.xlim(beg,en)

    plt.hlines(0,0,int(covloc[-1][1]), colors="black")
    plt.xlabel("Position on the genome")
    plt.ylabel("CircRNA abundance\n(RPB)")
    plt.savefig("./Circ_List_{name}.png".format(name=str(id_list[i]).replace(".","_").replace("/","_")),dpi=1000,bbox_inches='tight')
    plt.clf()