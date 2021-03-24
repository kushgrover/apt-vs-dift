import random
# import os
# import re
import sys

# Input arguments
# 1 : Output file name
# 2 : Size

f = open(sys.argv[1],"w")

# NUMBER OF STATES
n = int(sys.argv[2]) - 1

# GLOBAL VARIABLES ETC
f.write("smg\n\n")
f.write("const int n = " + str(n) + ";\n")
f.write("global s: [0 .. (n+3)];    \t// order: d, tau, phi\n")
f.write("global p: [0 .. 1];\n\n")

# DEFENDER
f.write("player defender\n")
f.write("\t[trap], [notrap]\n")
f.write("endplayer\n\n")

# ATTACKER
f.write("player attacker\n")
f.write("\tgame, [dropout], [target]")
# for i in range(1,n):
#     f.write(", [a" + str(i) + "]")
f.write("\nendplayer\n\n")

# MODULE GAME 
f.write("module game\n")

# defender actions
f.write("\t// defender actions\n")
f.write("\t[notrap] s<=n-1 & p=0 -> (p'=1);\n")

for i in range(n):
    r = random.uniform(0.5,1) # trapping prob
    f.write("\t[trap] s=" + str(i) + " & p=0 -> " + str(r) + " : (s'=n+1) + " + str(1-r) + " : (p'=1);\n")

# sink transitions
f.write("\n\t// sink transitions\n")
f.write("\t[] s=n -> (s'=n+3);\n")
f.write("\t[] s=n+1 -> (s'=n+3);\n")
f.write("\t[] s=n+2 -> (s'=n+3);\n\n")

# attacker actions
f.write("\t// attacker actions\n")
f.write("\t[dropout] s<=n-1 & p=1 -> (s'=n+2);\n")

for i in range(1,n):
    try:
        if(sys.argv[3] == "-acyclic"):
            numofpre = random.randint(1,i)
            pre = random.sample(range(i), k = numofpre)
        else:
            print("wrong switch, continuining without switch")
    except: 
        numofpre = random.randint(1,n)
        pre = random.sample(range(n), k = numofpre)
    f.write("\t[] (")
    for j in pre[:-1]:
        f.write("s=" + str(j) + " | ")
    f.write("s=" + str(pre[-1]))
    f.write(") & p=1 -> (s'=" + str(i) + ") & (p'=0);\n")

# reaching target action
numofpre = random.randint(1,n)
pre = random.sample(range(n), k = numofpre)
f.write("\t[target] (")
for j in pre[:-1]:
    f.write("s=" + str(j) + " | ")
f.write("s=" + str(pre[-1]))
f.write(") & p=1 -> (s'=n);\n")

f.write("endmodule\n\n")
# END MODULE MAIN


# REWARDS
f.write("rewards \"targetreached\"\n")
f.write("\ts=n : 1;\n")
f.write("endrewards\n\n")

f.write("rewards \"trapped\"\n")
f.write("\ts=n+1 : 1;\n")
f.write("endrewards\n\n")

# COST
f.write("rewards \"cost\"\n")
f.write("\t[trap] true: 1;\n")
f.write("\t[notrap] true: 0;\n")
f.write("endrewards\n\n")

f.close()

