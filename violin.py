import matplotlib.pyplot as plt
import numpy as np

covcirirb1b=open("ciricovrb1b.csv")
covtraprb1b=open("trapcovrb1b.csv")
covciriadeno=open("ciricovadeno.csv")
covtrapadeno=open("trapcovadeno.csv")
covciriblv=open("ciricovblv.csv") #done
covtrapblv=open("trapcovblv.csv") #done
covcirihtlv=open("ciricovhtlv.csv") #done
covtraphtlv=open("trapcovhtlv.csv") #done
covcirialhv=open("ciricovalhv.csv") #done
covtrapalhv=open("trapcovalhv.csv") #done

covvirrb1btrap=open("./covvirrb1btrap.csv")
covvirrb1bciri=open("./covvirrb1bciri.csv")

covcirirb1b=covcirirb1b.readlines()
covtraprb1b=covtraprb1b.readlines()
covvirrb1btrap=covvirrb1btrap.readlines()
covvirrb1bciri=covvirrb1bciri.readlines()

tottrap=0
totciri=0
for line in covvirrb1btrap:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    tottrap=tottrap+line[2]

for line in covvirrb1bciri:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    totciri=totciri+line[2]


xpoints=[]
ypoints=[]
plt.figure(figsize=(12,1))
for a in covcirirb1b:
    a=a.split("\t")
    xpoints.append(int(a[1]))
    ypoints.append(-((int(a[2]))/totciri)*1000000000)
plt.plot(xpoints,ypoints,color='red',linewidth=0.5)
xpoints=[]
ypoints=[]
for a in covtraprb1b:
    a=a.split("\t")
    xpoints.append(int(a[1]))
    ypoints.append((int(a[2])))

plt.plot(xpoints,ypoints,color='blue',linewidth=0.5)
plt.hlines(0,0,int(covvirrb1bciri[-1].split("\t")[1]), colors="black")
plt.xlim((int(covcirirb1b[-1].split("\t")[1])/20),(int(covcirirb1b[-1].split("\t")[1])/20)*3)
plt.ylim(-200,600)
plt.vlines((((int(covcirirb1b[-1].split("\t")[1])/20)+((int(covcirirb1b[-1].split("\t")[1])/20)*3))/2),-200,600,colors="black")
plt.savefig("./Circ_Coverage_{name}_plotrb1bdiv20.png".format(name=str("rb1btest").replace(".","_").replace("/","_")),dpi=1000,bbox_inches='tight')
