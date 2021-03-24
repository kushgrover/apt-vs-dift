#!/usr/bin/python3

import sys
import csv
import os

log_folder = sys.argv[1]

print("Model;Approach;Repetition;Resets;stepsUntilLastReset;lastCandidatePrint")

for file in os.listdir(log_folder):
	if not ".log" in file:
		continue
	model = file.split("-")[1]
	approach = file.split("-")[3]
	rep = file.split("-")[-1].split(".")[0]
	resets="TO"
	stepsUntilLastReset="TO"
	lastline="init"
	with open(log_folder+"/"+file,"r") as f:
		for line in f:
			if line.lstrip().startswith("Simulation result details:"):
				resets = ("".join(line.split(":")[1:]).strip()).split(" ")[1]
				stepsUntilLastReset = ("".join(line.split(":")[1:]).strip()).split(" ")[3]
			if line.lstrip().startswith("Resets:"):
				lastline = line
	print(";".join([model,approach,rep,resets,stepsUntilLastReset,lastline]))
