#!/usr/bin/python3

import sys
import csv
import os

target_file = sys.argv[1]

#dict mapping model cross approach tuples to total number of resets over all repetitions resp. the total number of steps
mxa2result = dict()


with open(target_file,"r") as file:
	for line in file:
		if line.startswith("Model") or "o" not in line:
			continue
		model,app,rep,resets,steps,lastline = line.split(";")

		if model+";"+app not in mxa2result.keys():
			result = dict()
			result["resets"]=0
			result["steps"]=0
			result["samples"]=0
			result["TOs"]=0
			result["lines"]=""
			mxa2result[model+";"+app]=result

		result = mxa2result[model+";"+app]
		if resets == "TO":
			result["TOs"]+=1
		else:
			result["samples"]+=1
			result["resets"] += int(resets)
			result["steps"] += int(steps)
		result["lines"] += lastline.strip()+";"
		mxa2result[model+";"+app]=result


print("Model;Approach;Avg Resets;Avg StepsUntilReset;TOs;Lines...")
for mxa in mxa2result.keys():
	res = mxa2result[mxa]
	tos = res["TOs"]
	if res["samples"]==0:
		print(mxa+";-1;-1;"+str(tos)+";"+res["lines"])
	else:
		denominator=res["samples"]
		print(mxa+";"+str(res["resets"]/denominator)+";"+str(res["steps"]/denominator)+";"+str(tos)+";"+res["lines"])



