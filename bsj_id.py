#Splicing sites identificator and comparator
# count_file = open("./lat_bsj_count.csv")
sites_file = open("./aln_bsj_sites.txt") #this file contains all the backsplice junctions and backsplice sites and the identifiers of the associated reads
out_file = open("./aln_new_count.csv",'w') #this file will save the count of junctions identified

# count_lines=count_file.readlines()
sites_line=sites_file.readlines()

from difflib import *


bsj_list=[]
one_more_list=[]
index = -1
for line in sites_line: #for every read spanning a backsplice junction
    id_list=[]
    index=index+1
    if index%50==0:
        print(str((index/len(sites_line))*100)+"% " + "of the sequences were analyzed for counting")
    line=line.split("\t")
    if (len(line[1]) > 25): #sometimes the sequences are trimmed for unknown reasons, don't keep them
        bsj_count=1
        if not (line[2:4] in bsj_list):
            for line2 in sites_line[index+1:]: #compare the backsplice sites and junctions with every other read to count the ones that are similar
                line3=line2
                line2=line2.split("\t")
                diffa = SequenceMatcher(None, line[2], line2[2])
                diffb = SequenceMatcher(None, line[3], line2[3])
                diff2a=diffa.ratio()
                diff2b=diffb.ratio()
                diffc = SequenceMatcher(None, line[1], line2[1])
                diff2c=diffc.ratio()
                if diff2c > 0.95 and (line2[0] != line[0]) and (not line2[0]+"1" in id_list) and (not line2[0]+"2" in id_list):
                    bsj_count=bsj_count+1
                    one_more_list.append(line2[0][:-1])
                    del sites_line[sites_line.index(line3)]
                    id_list.append(line2[0])
                elif ((diff2a > 0.97) or (line[2][4:-4] in line2[2])) and ((diff2b > 0.97) or (line[3][4:-4] in line2[3])) and (len(line2[1])>25) and (line2[0] != line[0]) and (not line2[0]+"1" in id_list) and (not line2[0]+"2" in id_list):
                    bsj_count=bsj_count+1
                    one_more_list.append(line2[0])
                    del sites_line[sites_line.index(line3)]
                    id_list.append(line2[0])
            bsj_list.append(line[2:4])
            if (not line[0] in id_list) and (not line[0]+"1" in id_list) and (not line[0]+"2" in id_list):
                id_list.append(line[0])

            line[3] = line[3].replace("\n","")
            if not ("GGTTAGGG" in (line[1] + line[2] + line[3])) and not ("TTAGGGTT" in (line[1] + line[2] + line[3])):
                out_file.write(line[1]+"\t"+line[2] + "\t" + line[3] + "\t" + str(bsj_count)+"\t"+str(id_list)+"\n")