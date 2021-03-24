import random
# import os
# import re
import sys
import pandas as pd

csv_reader = pd.read_csv(sys.argv[1])
x1 = csv_reader['0']

output = open(sys.argv[2], "w")

# NUMBER OF STATES
n = x1.size-1

# GLOBAL VARIABLES ETC
output.write("smg\n\n") 
output.write("const int n = " + str(n) + ";\n")
output.write("global s: [0 .. (n+3)];    \t// order: d, tau, phi\n")
output.write("global p: [0 .. 1];\n\n")

# DEFENDER
output.write("player defender\n")
output.write("\t[trap], [notrap]\n")
output.write("endplayer\n\n")

# ATTACKER
output.write("player attacker\n")
output.write("\tgame, [dropout], [target]")
output.write("\nendplayer\n\n")

# MODULE GAME 
output.write("module game\n")

# defender actions
output.write("\t// defender actions\n")
output.write("\t[notrap] s<=n-1 & p=0 -> (p'=1);\n")

for i in range(n):
    r = random.uniform(0.8,1) # trapping prob
    output.write("\t[trap] s=" + str(i) + " & p=0 -> " + str(r) + " : (s'=n+1) + " + str(1-r) + " : (p'=1);\n")

# sink transitions
output.write("\n\t// sink transitions\n")
output.write("\t[] s=n -> (s'=n+3);\n")
output.write("\t[] s=n+1 -> (s'=n+3);\n")
output.write("\t[] s=n+2 -> (s'=n+3);\n\n")

# attacker actions
output.write("\t// attacker actions\n")
output.write("\t[dropout] s<=n-1 & p=1 -> (s'=n+2);\n")

def compute_pre(x):
    pre = []
    for i in range(n-1):
        if(x[i] == 1):
            pre.append(i)
    return pre

for i in range(1,n-1):
    x = csv_reader[str(i)]
    pre = compute_pre(x)
    if(pre != []):
        output.write("\t[] (")
        for j in pre[:-1]:
            output.write("s=" + str(j) + " | ")
        output.write("s=" + str(pre[-1]))
        output.write(") & p=1 -> (s'=" + str(i) + ") & (p'=0);\n")

# reaching target action
x = csv_reader[str(n)]
pre = compute_pre(x)
output.write("\t[target] (")
for j in pre[:-1]:
    output.write("s=" + str(j) + " | ")
output.write("s=" + str(pre[-1]))
output.write(") & p=1 -> (s'=n);\n")

output.write("endmodule\n\n")
# END MODULE MAIN

# REWARDS
output.write("rewards \"targetreached\"\n")
output.write("\ts=n : 1;\n")
output.write("endrewards\n\n")

output.write("rewards \"trapped\"\n")
output.write("\ts=n+1 : 1;\n")
output.write("endrewards\n\n")

# COST
output.write("rewards \"cost\"\n")
output.write("\t[trap] true: 1;\n")
output.write("\t[notrap] true: 0;\n")
output.write("endrewards\n\n")

output.close()

