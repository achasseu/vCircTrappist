#Correction_SAM
import sys
import getopt

argv=sys.argv[1:]
optlist, args= getopt.getopt(argv,"I:S:O:")
samnc=None
samc =None

for opt, arg in optlist:
    if opt in "-I":
        ref = open(arg)
    if opt in "-S":
        samnc = open(arg)
    if opt in "-O":
        samc = open (arg,"w")

if samnc == None or samc==None:
    print("You forgot your arguments".upper())


refa=ref.readlines()
samnc=samnc.readlines()

for line in refa[:50]:
    if "@" in line[0]:
        samc.write(line)
ref.close()
for line in samnc:
    samc.write(line)

    


