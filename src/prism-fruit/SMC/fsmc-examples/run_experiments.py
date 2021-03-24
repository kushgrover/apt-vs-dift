#!/usr/bin/python3
import subprocess
import re
import os, sys
import time
import pipes
import argparse

TIMEOUT="15m"

patternP=re.compile("^(P.+)$")
CPUS=[15,16,17,18,19,30,31,32,33,34,35,36,37,38,39]
NEXTCPU=0

def run(model, prop, const, minprob, arguments, num_arg, scriptname, logname):
    global NEXTCPU
    global CPUS
    cmd = "/usr/bin/time -v timeout -s 9 "+TIMEOUT+" ../prism/bin/prism "+ pipes.quote(model)+" "+prop+" -simminprob "+minprob+" "+arguments+" "+num_arg;
    if len(const) > 0:
        cmd = cmd + " -const "+const
    #print("="*30)
    #print("Running "+cmd)
    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print(p.stdout.read().decode('UTF-8'))
    CPUID = CPUS[NEXTCPU]
    print(f"pueue add -g cpu{CPUID} -- \"taskset -c {CPUID} {cmd} > '{logname}.log' 2>&1\"")
    NEXTCPU = (NEXTCPU+1)%len(CPUS)
                    
    

def generate(benchname, expname,repeat):
    benchfile = open(benchname,'r')
    benchmarks = []
    for line in benchfile.readlines():
        elems = line.split(';')
        if len(elems)>=4 and elems[0][0] != "#":
            tpl = (elems[0], elems[1], elems[2], elems[3], elems[4].rstrip())
            benchmarks.append(tpl)

    # read experiments
    expfile = open(expname,'r')
    experiments = []
    for line in expfile.readlines():
        elems = line.split(';')
        if len(elems)==3 and elems[0][0] != "#":
            tpl = (elems[0], elems[1],elems[2].rstrip())
            experiments.append(tpl)

    # do experiments
    for experiment in experiments:
        (short,name,arguments) = experiment
        for benchmark in benchmarks:
            (model,prop,const,minprob,num_arg) = benchmark
            modelName = os.path.basename(model)
            pfile = open(prop, 'r')
            plines = pfile.read()
            matchP = patternP.search(plines)

            if (matchP):
                propName = matchP.group(1)
            else:
                propName = str(os.path.basename(prop))

            for r in range(repeat):
                logname = "logs/log-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short+"-"+str(r)
                # resname = "result/res-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short
                scriptname = "scripts/script-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short+"-"+str(r)
                run(model, prop, const, minprob, arguments, num_arg, scriptname, logname)
                # print(scriptname)

    

def main(genscripts, readscripts):
    if len(sys.argv) < 4:
        print(USAGE, file=sys.stderr)
        exit(1)

    outdir = sys.argv[3]    
        
    repeat = 1
    if len(sys.argv) == 5:
            repeat = int(sys.argv[4])

    # read benchmarks
    benchfile = open(sys.argv[1],'r')
    benchmarks = []
    for line in benchfile.readlines():
        elems = line.split(';')
        if len(elems)>=4 and elems[0][0] != "#":
            tpl = (elems[0], elems[1], elems[2], elems[3], elems[4].rstrip())
            benchmarks.append(tpl)

    # read experiments
    expfile = open(sys.argv[2],'r')
    experiments = []
    for line in expfile.readlines():
        elems = line.split(';')
        if len(elems)==3 and elems[0][0] != "#":
            tpl = (elems[0], elems[1],elems[2].rstrip())
            experiments.append(tpl)

    # do experiments
    for experiment in experiments:
        (short,name,arguments) = experiment
        
        for benchmark in benchmarks:
            (model,prop,const,minprob,num_arg) = benchmark
            modelName = os.path.basename(model)
            pfile = open(prop, 'r')
            plines = pfile.read()
            matchP = patternP.search(plines)

            if (matchP):
                propName = matchP.group(1)
            else:
                propName = str(os.path.basename(prop))

            logname = outdir+"/log-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short
            resname = outdir+"/res-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short
            scriptname = outdir+"/script-"+modelName+"("+const+")-"+os.path.basename(prop)+"-"+short
            logfile = open(logname,'w')
            resfile = open(resname,'w')
            original = sys.stdout
            sys.stdout = Tee(sys.stdout, resfile)

            print("---------------------------------------")
            print(name)
            
            if len(const)>0:
                print("Checking "+str(modelName)+"("+const+") with property "+propName)
                
            else:
                print("Checking "+str(modelName)+" with property "+propName)
            print()

            result = {}

            for i in range(repeat):
                

                if genscripts:
                    generateScript(model, prop, const, minprob, arguments, num_arg, scriptname, logname)
                    continue

                print("[try="+str(i)+"] ", end=' ')
                pres = doSim(model, prop, const, minprob, arguments, num_arg, logfile)
                
                if 'result' not in pres:
                    print("ERROR: Couldn't compute statistics", file=sys.stderr)
                    exit(1)

                if pres['result'] == 'TO':
                    result['result'] = 'TO'
                    print()
                    break

                if 'result' not in result:
                    result['result'] = pres['result']
                elif result['result'] != pres['result']:
                    result['result'] = 'MISMATCH'

                for key in ["time", "path_len_avg", "reach_anal", "inclusion", "iterations","e2b","trans"]:
                    if (key in pres):
                        if key not in result:
                            result[key] = pres[key]/repeat
                        else:
                            result[key] += pres[key]/repeat

                print()

            if genscripts:
                return
            # print stats
            print("[Summary]", end=' ')
            if result['result'] == 'TO':
                print("result=TO")
            else:
                for key in list(result.keys()):
                    print(key+"="+str(result[key])+"\t", end=' ')
            
            print()
            print()

            
            logfile.close()
            sys.stdout = original
            resfile.close()



            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulations for SMC.')
    parser.add_argument('-e','--experiments',nargs=1,help="experiments file")
    parser.add_argument('-b','--benchmarks',nargs=1,help="benchmarks file")
    parser.add_argument('-r','--repeat',type=int,default=1,help="how many repetitions")
    args = parser.parse_args()           

    if (args.benchmarks == None or args.experiments == None):
        parser.print_help()
        exit(1)
    generate(args.benchmarks[0], args.experiments[0],args.repeat)
