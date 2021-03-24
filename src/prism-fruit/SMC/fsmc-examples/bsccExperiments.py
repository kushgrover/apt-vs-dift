#!/usr/bin/python
import subprocess
import re
import os

paternRes=re.compile("^result:\t\t(.*)$", re.MULTILINE)
paternBSCC=re.compile("^sim\. BSCCs:\t(.*)$", re.MULTILINE)

files= ["/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/brp/brp(5,5)_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/brp/brp(20,20)_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/leader/synchronous/leader4_2_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/leader/synchronous/leader6_2_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/leader/synchronous/leader4_6_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/leader/synchronous/leader5_6_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/self-stabilisation/herman/herman9_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/self-stabilisation/herman/herman13_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/contract_eql/contract.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/crowds/crowds_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/nand/nand(5,5)_more.pm",
"/home/przemek/IST/projects/StatEquiv/prism-stat-eq/examples/nand/nand(10,10)_more.pm"]


#confs= [(500,100),(500,1000),(500,300),(5000,300),(50000,300)]
confs= [(500,300)]



for con in confs:
    (sampleNo,visitNo) = con
    print("Samples: "+str(sampleNo)+", visitNo: "+str(visitNo)+"\n")

    for f in files:
        name = os.path.basename(f)
        print("processing "+name)
        cmd = ["/home/przemek/IST/projects/StatEquiv/prism-stat-eq/bin/prism","-simbscc", "-simsamples",str(sampleNo),"-simvisitno",str(visitNo),f]
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        output = p1.communicate()[0]
        matchRes = paternRes.search(output, re.MULTILINE)
        matchBSCC = paternBSCC.search(output, re.MULTILINE)

        if matchRes:
            print(matchRes.group(1))
        if matchBSCC:
            print(matchBSCC.group(1))

        if not (matchRes or matchBSCC):
            print("ERROR: Couldn't parse output!");
            print(output) 
            exit(1)
            


    print("\n\n")
    
