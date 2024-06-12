import matplotlib.pyplot as plt
import pandas as pd

totcov=open("./aln_virus_coverage.csv")
totcov=totcov.readlines()
sites_sorting=pd.read_csv("./sites_sorting.csv",sep="\t",header=0)
sites_sorting.sort_values("Count", axis=0, ascending=False,inplace=True, na_position='first')
first_circ=sites_sorting[:51]
first_circ=first_circ.loc[:,["Count","5'_Position_(Sense_Strand)","3'_Position_(Sense_Strand)","Chromosome","Sense_Score","U2_Score","Max_CircRNA_Size"]]
circs_list=first_circ.values.tolist()



plt.figure(figsize=(12,4))
for circ in circs_list:
    
    size=int(circ[-1])/2
    position=(int(circ[2])+int(circ[1]))/2
    abundance=int(circ[0])
    if circ[4] > 0 :
        if circ[5] > 0.6:
            plt.scatter(position,size,abundance,c='black',linewidth=0,alpha=0.5)
            
        else :
            plt.scatter(position,size,abundance,c="darkorange",linewidth=0,alpha=0.5)
    else :
        if circ[5] > 0.6:
            plt.scatter(position,size,abundance,c='black',linewidth=0,alpha=0.5)
        else :
            plt.scatter(position,size,abundance,c='darkorange',linewidth=0,alpha=0.5)


for line in totcov[1:]:
    line=line.split("\t")
    plt.hlines(0,0,int(line[1]), colors="black")

plt.savefig("./Circ_List_{name}.png".format(name=str(totcov[-1].split("\t")[0])),dpi=1000,bbox_inches='tight')
    