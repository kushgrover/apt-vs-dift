#!/usr/bin/python3
# Usage:
# ./run_experiments ... > result.log
# ./result_parser result.log

import sys
import csv

log_file = sys.argv[1]
f = open(log_file)

keyphrases = ["Command line",
            "Number of cnd. (max, avg.)", 
            "Size of largest candidate", 
            "Result", 
            "User time (seconds)", 
            "Percent of CPU this job got",
            "Maximum resident set size (kbytes)",
            "Simulation result details:"]

results = []

for line in f:
    if line.startswith("PRISM"):
        results.append(dict())
    if line.startswith("Error"):
        results[-1]["error"] = True
    for phrase in keyphrases:
        if line.lstrip().startswith(phrase):
            if phrase == "Command line":
                command = line.split(":")[1].strip().split()
                results[-1]["model"] = command[1].split("/")[-1]
                results[-1]["property"] = command[2].split("/")[-1]
                results[-1]["const"] = command[-1]
                results[-1]["pmin"] = command[command.index("-simpmin")+1]                      if "-simwhitebox" in command:
                    results[-1]["whitebox"] = True
                else:
                    results[-1]["whitebox"] = False
            else:
                results[-1][phrase] = line.split(":")[1].strip()

names = ["model", "property", "const", "pmin", "whitebox", "error"] + keyphrases[1:]
writer = csv.DictWriter(sys.stdout, fieldnames=names)
writer.writeheader()
writer.writerows(results)
