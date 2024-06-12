import matplotlib.pyplot as plt
import numpy as np

ssU2cov=open("./circ_sense_U2_coverage.csv")
ssnU2cov=open("./circ_sense_nonU2_coverage.csv")
asU2cov=open("./circ_antisense_U2_coverage.csv")
asnU2cov=open("./circ_antisense_nonU2_coverage.csv")

ssU2cov=ssU2cov.readlines()
ssnU2cov=ssnU2cov.readlines()
asU2cov=asU2cov.readlines()
asnU2cov=asnU2cov.readlines()

total=0
for line in ssU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]
for line in ssnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]
for line in asU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]
for line in asnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    total=total+line[2]

ssU2cova=[]
ssnU2cova=[]
asU2cova=[]
asnU2cova=[]
for line in ssU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*1000000
    ssU2cova.append(line)
for line in ssnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*1000000
    ssnU2cova.append(line)
for line in asU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*-1000000
    asU2cova.append(line)
for line in asnU2cov:
    line=line.replace("\n","")
    line=line.split("\t")
    line[2]=int(line[2])
    line[2]=(line[2]/total)*-1000000
    asnU2cova.append(line)


id_list=[]
for line in (ssU2cova+ssnU2cova+asU2cova+asnU2cova):
    if not (line[0] in id_list):
        id_list.append(line[0])

chrnb=len(id_list)
for i in range(chrnb):
    xpoints=[]
    ypoints=[]
    for line in ssU2cova:
        if line[0]==id_list[i]:
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='k')
    xpoints=[]
    ypoints=[]
    for line in ssnU2cova:
        if line[0]==id_list[i]:
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='g')
    xpoints=[]
    ypoints=[]
    for line in asU2cova:
        if line[0]==id_list[i]:
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))            
    plt.plot(xpoints, ypoints,color='k')
    xpoints=[]
    ypoints=[]
    for line in asnU2cova:
        if line[0]==id_list[i]:
            xpoints.append(int(line[1]))
            ypoints.append(int(line[2]))
    plt.plot(xpoints, ypoints,color='g')
    plt.savefig("./supergraph_{name}.png".format(name=str(id_list[i])))
    plt.clf()
