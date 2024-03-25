import matplotlib.pyplot as plt
import pandas as pd
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"b:e:o:")
for opt, arg in optlist:
    if opt in "-o":
        totcov=open(str(arg)+"/aln_virus_coverage.csv")
        totcov=totcov.readlines()
        sites_sorting=pd.read_csv(str(arg)+"/sites_sorting.csv",sep="\t",header=0)
        sites_sorting.sort_values("Count", axis=0, ascending=False,inplace=True, na_position='first')
        first_circ=sites_sorting[:102]
        first_circ=first_circ.loc[:,["Count","5'_Position_(Sense_Strand)","3'_Position_(Sense_Strand)","Chromosome","Sense_Score","U2_Score","Max_CircRNA_Size"]]
        circs_list=first_circ.values.tolist()
        totcov=open(str(arg)+"/aln_virus_coverage.csv")
        totcov=totcov.readlines()


#determining the number and identities of the viral genomic segments
id_list=[]
for line in (circs_list):
    if not (line[3] in id_list):
        id_list.append(line[3])


#function to draw the arks representing the backsplice junctions and their abundance in the database
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

#Determining the mean coverage on the viral genome --> normalization
total=0
for line in totcov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]

#Drawing the sashimi plot. The abundance of the junction-mapping reads is represented by the height of the arks.
chrnb=len(id_list)
for i in range(chrnb):
    plt.figure(figsize=(12,4))
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

    plt.hlines(0,0,int(covloc[-1][1]),colors="black",label="Genome",linewidth=0.8)
    plt.xlabel("Position on the genome")
    plt.ylabel("CircRNA abundance\n(RPB)")
    for opt, arg in optlist:
        if opt in "-o":
            plt.savefig(str(arg)+"/Circ_List_{name}.png".format(name=str(id_list[i]).replace(".","_").replace("/","_")),dpi=1000,bbox_inches='tight')
    plt.clf()