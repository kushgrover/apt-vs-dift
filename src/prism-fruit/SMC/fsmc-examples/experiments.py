#!/usr/bin/python
import subprocess
import re
import os, sys
import time
import pipes
import argparse
USAGE = "Usage: "+sys.argv[0] + " BENCHMARK_FILE EXPERIMENTS_FILE OUTPUTDIR [REPEAT]"
TIMEOUT="15m"
SERVER=True

paternRes=re.compile("^Result:\s+(\S*)\s.*$", re.MULTILINE)
patternP=re.compile("^(P.+)$")
paternSimParams=re.compile("^Simulation method parameters:.*$", re.MULTILINE)
paternParams=re.compile("min. prob=.+, iterations=.+$", re.MULTILINE)
paternIters=re.compile("^Sampling complete: (.+) iterations in (.+) seconds", re.MULTILINE)
paternFalseCnd=re.compile("^No of false candidates:\s*(.+)$", re.MULTILINE)
paternError=re.compile("^Error:.+$", re.MULTILINE);
paternVN=re.compile("Visit no. \(base, increment\):\s*(.+)$", re.MULTILINE)
paternCND=re.compile("Number of cnd. \(max, avg.\):\s*(.+)$", re.MULTILINE)
paternReachTime=re.compile("^Time for reachability analysis:\s*Time=(.+) secs$", re.MULTILINE)
paternInclTime=re.compile("^Time for inclusion check:\s*Time=(.+) secs$", re.MULTILINE)
paternPathlen=re.compile("Path length statistics: average (.+), min (.+), max (.+)" ,re.MULTILINE)
paternMem = re.compile('Maximum resident set size \(kbytes\): (.*)')
paternTrans = re.compile("Avg. \(transient,BSCC\) path len.:\s*\((.+), (.+)\)")
paternETB = re.compile("Avg. entry to birth time:\s*(.+)")

paternTime = re.compile("Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\d+):(\d+)\.(\d+)")
#paternTime = re.compile("User time \(seconds\):\s*(.+)")
paternLog = re.compile(".*log-(.+\(.*\))-(.+)-(.+)-([^-]+)")

model_data = {
'bluetooth_oneinit.pm(mrec=4)' : ('bluetooth(4)', '149K', '$7.8\cdot10^{-3}$', '3K, 1'),
'bluetooth_oneinit.pm(mrec=7)' : ('bluetooth(7)', '569K', '$7.8\cdot10^{-3}$', '5.8K, 1'),
'bluetooth_oneinit.pm(mrec=10)' : ('bluetooth(10)', '>569K', '$7.8\cdot10^{-3}$', '>5.8K, 1'),
'brp_nodl.pm(N=10000-MAX=10000)' : ('brp(10K,10K)', '>40M ', '$0.01$', '>4.5K, 1'),
'brp_nodl.pm(N=1500-MAX=1500)' : ('brp(1.5K,1.5K)', '? ', '$0.01$', '4.5K, 1'),
'brp_nodl.pm(N=2000-MAX=2000)' : ('brp(2K,2K)', '40M ', '$0.01$', '4.5K, 1'),
'brp_nodl.pm(N=500-MAX=500)' : ('brp(500,500)', '4.5M ', '$0.01$', '1.5K, 1'),
'crowds_nodl.pm(TotalRuns=6-CrowdSize=15)' : ('crowds(6,15)', '7.3M', '$0.066$', '>3K, 1'),
'crowds_nodl.pm(TotalRuns=7-CrowdSize=20)' : ('crowds(7,20)', '17M', '$0.05$', '>3K, 1'),
'crowds_nodl.pm(TotalRuns=8-CrowdSize=20)' : ('crowds(8,20)', '68M', '$0.05$', '>3K, 1'),
'egl.pm(N=15-L=10)' : ('eql(15,10)', '616G ', '$0.5$', '1, 1'),
'egl.pm(N=20-L=15)' : ('eql(20,15)', '1279T', '$0.5$', '1, 1'),
'egl.pm(N=20-L=20)' : ('eql(20,20)', '1719T', '$0.5$', '1, 1'),
'extreme.pm(N=18)' : ('Fig.\\ref{fig:self-looping states}(18)', '37', '$0.5$', '1, 1,'),
'extreme.pm(N=19)' : ('Fig.\\ref{fig:self-looping states}(19)', '39', '$0.5$', '1, 1'),
'extreme.pm(N=20)' : ('Fig.\\ref{fig:self-looping states}(20)', '41', '$0.5$', '1, 1'),
'fms.sm(n=10)' : ('fms(10)', '234M', '$4\cdot10^{-4}$', '1, 25M'),
'fms.sm(n=20)' : ('fms(20)', '>234M', '$2\cdot10^{-4}$', '1, >25M'),
'fms.sm(n=40)' : ('fms(40)', '>234M', '$1\cdot10^{-4}$', '1, >25M'),
'gridworld.pm(N=300-p=0.999)' : ('gridworld(300)', '162M', '$1\cdot10^{-3}$', '598, 89K'),
'gridworld.pm(N=400-p=0.999)' : ('gridworld(400)', '384M', '$1\cdot10^{-3}$', '798, 160K'),
'gridworld.pm(N=500-p=0.999)' : ('gridworld(500)', '750M', '$1\cdot10^{-3}$', '998, 250K'),
'herman17_oneinit.pm()' : ('herman(17) ', '129M', '$7.6\cdot10^{-6}$', '1, 34'),
'herman19_oneinit.pm()' : ('herman(19) ', '1162M', '$1.9\cdot10^{-6}$', '1, 38'),
'herman21_oneinit.pm()' : ('herman(21) ', '10G', '$4.7\cdot10^{-7}$','1, 42'),
'leader6_11.pm()' : ('leader(6,11)', '>280K', '$5.6\cdot10^{-7}$', '1, 1'),
'leader6_5.pm()' : ('leader(6,5)', '?', '?', '?'),
'leader6_6.pm()' : ('leader(6,6)', '280K ', '$2.1\cdot10^{-5}$', '1, 1'),
'leader6_8.pm()' : ('leader(6,8)', '>280K', '$3.8\cdot10^{-6}$', '1, 1'),
'nand.pm(N=50-K=3)' : ('nand(50,3)', '11M', '$0.02$', '51, 1'),
'nand.pm(N=60-K=4)' : ('nand(60,4)', '29M', '$0.02$', '61, 1'),
'nand.pm(N=70-K=5)' : ('nand(70,5)', '67M', '$0.02$', '71, 1'),
'scale.pm(N=1000-M=5)' : ('Fig.\\ref{fig:add}(1K,5)', '4022', '$0.5$', '2, 5'),
'scale.pm(N=1000-M=50)' : ('Fig.\\ref{fig:add}(1K,50)', '4202', '$0.5$', '2, 50'),
'scale.pm(N=1000-M=500)' : ('Fig.\\ref{fig:add}(1K,500)', '6002', '$0.5$', '2, 500,'),
'scale.pm(N=10000-M=5)' : ('Fig.\\ref{fig:add}(10K,5)', '40K', '$0.5$', '2, 5'),
'scale.pm(N=100000-M=5)' : ('Fig.\\ref{fig:add}(100K,5)', '400K', '$0.5$', '2, 5'),
'tandem.sm(c=1000)' : ('tandem(1K)', '1.7M', '$9.9\cdot10^{-5}$', '1, 501K'),
'tandem.sm(c=2000)' : ('tandem(2K)', '>1.7M', '$4.9\cdot10^{-5}$ ', '1, >501K'),
'tandem.sm(c=500)' : ('tandem(500)', '>1.7M', '$2.4\cdot10^{-5}$ ', '1, >501K'),
'triangle(12).pm(M=100-p=0.8)' : ('triangle(12)', '?', '?', '?'),
'triangle(13).pm(M=100-p=0.8)' : ('triangle(13)', '?', '?', '?'),
'triangle(14).pm(M=100-p=0.8)' : ('triangle(14)','?', '?', '?')
}

correct_result = {
'bluetooth_oneinit.pm(mrec=10)' : {'time_qual.pctl' : True },
'bluetooth_oneinit.pm(mrec=4)' : {'time_qual.pctl' : True },
'bluetooth_oneinit.pm(mrec=7)' : {'time_qual.pctl' : True },
'brp_nodl.pm(N=10000-MAX=10000)' : {'p1_qual.pctl' : False },
'brp_nodl.pm(N=1500-MAX=1500)' : {'p1_qual.pctl' : False },
'brp_nodl.pm(N=2000-MAX=2000)' : {'p1_qual.pctl' : False },
'brp_nodl.pm(N=500-MAX=500)' : {'p1_qual.pctl' : False },
'crowds_nodl.pm(TotalRuns=6-CrowdSize=15)' : {'positive_qual.pctl' : True },
'crowds_nodl.pm(TotalRuns=7-CrowdSize=20)' : {'positive_qual.pctl' : True },
'crowds_nodl.pm(TotalRuns=8-CrowdSize=20)' : {'positive_qual.pctl' : True },
'egl.pm(N=15-L=10)' : {'unfairA_qual.pctl' : True },
'egl.pm(N=20-L=15)' : {'unfairA_qual.pctl' : True },
'egl.pm(N=20-L=20)' : {'unfairA_qual.pctl' : True },
'extreme.pm(N=18)' : {'prop1_qual.pctl' : True },
'extreme.pm(N=19)' : {'prop1_qual.pctl' : True },
'extreme.pm(N=20)' : {'prop1_qual.pctl' : True },
'fms.sm(n=10)' : {'reach_qual.pctl' : True },
'fms.sm(n=20)' : {'reach_qual.pctl' : True },
'fms.sm(n=40)' : {'reach_qual.pctl' : True },
'gridworld.pm(N=300-p=0.999)' : {'prop_qual.pctl' : True },
'gridworld.pm(N=400-p=0.999)' : {'prop_qual.pctl' : False },
'gridworld.pm(N=500-p=0.999)' : {'prop_qual.pctl' : False },
'herman17_oneinit.pm()' : {'reach.pctl' : True },
'herman19_oneinit.pm()' : {'reach.pctl' : True },
'herman21_oneinit.pm()' : {'reach.pctl' : True },
'leader6_11.pm()' : {'elected_qual.pctl' : True },
'leader6_5.pm()' : {'elected_qual.pctl' : True },
'leader6_6.pm()' : {'elected_qual.pctl' : True },
'leader6_8.pm()' : {'elected_qual.pctl' : True },
'nand.pm(N=50-K=3)' : {'reliable_qual.pctl' : True },
'nand.pm(N=60-K=4)' : {'reliable_qual.pctl' : True },
'nand.pm(N=70-K=5)' : {'reliable_qual.pctl' : True },
'scale.pm(N=1000-M=5)' : {'prop1_qual.pctl' : False },
'scale.pm(N=1000-M=50)' : {'prop1_qual.pctl' : False },
'scale.pm(N=1000-M=500)' : {'prop1_qual.pctl' : False },
'scale.pm(N=10000-M=5)' : {'prop1_qual.pctl' : False },
'scale.pm(N=100000-M=5)' : {'prop1_qual.pctl' : False },
'tandem.sm(c=1000)' : {'reach_qual.csl' : True },
'tandem.sm(c=2000)' : {'reach_qual.csl' : True },
'tandem.sm(c=500)' : {'reach_qual.csl' : True },
'triangle(12).pm(M=100-p=0.8)' : {'reach_qual.pctl' : True },
'triangle(13).pm(M=100-p=0.8)' : {'reach_qual.pctl' : True },
'triangle(14).pm(M=100-p=0.8)' : {'reach_qual.pctl' : True }
}



class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)



def generateScript(model, prop, const, minprob, arguments, num_arg, scriptname, logname):

    cmd = "/usr/bin/time -v timeout -s 9 "+TIMEOUT+" bin/prism "+ pipes.quote(model)+" "+prop+" -simminprob "+minprob+" "+arguments+" "+num_arg;
    if len(const) > 0:
        cmd = cmd + " -const "+const

    scriptfile = open(scriptname, 'w')
    scriptfile.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -M przemek@ist.ac.at
#$ -q main.q
#$ -N PRISM
#$ -j y
#$ -l h_rt=0:20:00
#$ -l h_vmem=10G
#$ -l mf=10G
#$ -wd /cluster/home/przemek/prism-stat-eq
""")
    logname = logname.replace(',','-')
    scriptfile.write("#$ -o \""+logname+"\"\n")

    modelName = os.path.basename(model)

    pfile = open(prop, 'r')
    plines = pfile.read()
    matchP = patternP.search(plines)

    if (matchP):
        propName = matchP.group(1)
    else:
        propName = str(os.path.basename(prop))
    pfile.close()
    
    # scriptfile.write("rm \""+logname+"\"\n")

    if len(const)>0:
        scriptfile.write("echo \"Checking "+str(modelName)+"("+const+") with property "+propName+"\"\n")        
    else:
        scriptfile.write("echo \"Checking "+str(modelName)+" with property "+propName+"\"\n")

    scriptfile.write(cmd)
    scriptfile.close()

                      

def parselog(log):
    logfile = open(log)
    txt = logfile.read()
    logfile.close()
    result = {}

    matchTime = paternTime.search(txt, re.MULTILINE)
    if matchTime:
        simTime = float(matchTime.group(1))*60+float(matchTime.group(2))
    else:
        print("Could not get time in "+log, file=sys.stderr)
        exit(1)
                       
          
    matchRes = paternRes.findall(txt, re.MULTILINE)
    if len(matchRes)>0:
        simRes = matchRes[-1]
        # print "res="+str(simRes)+"\t",
        # print "time="+("{:.2f}".format(simTime))+"s\t",
        result['result'] = simRes
        result['time'] = simTime

        memMatch = paternMem.search(txt, re.MULTILINE)
        if memMatch:
            size = int(memMatch.group(1))/1024
            # print "mem="+str(size)+"MB\t",

        matchFalseCnd = paternFalseCnd.search(txt, re.MULTILINE)
        if matchFalseCnd:
            fc = matchFalseCnd.group(1)
            # print "false cnd="+fc+"\t",

#        matchVN=paternVN.search(txt, re.MULTILINE)
#        if matchVN:
            # print "visit no. (base, inc)="+matchVN.group(1)+"\t",
#        matchCND=paternCND.search(txt, re.MULTILINE)
#        if matchCND:
            # print "no of cand. (max, avg)="+matchCND.group(1)+"\t",
        matchPathlen=paternPathlen.search(txt, re.MULTILINE)
        if matchPathlen:
            # print "path len (max,avg)=("+matchPathlen.group(3)+","+matchPathlen.group(1)+")\t",
            result['path_len_avg'] = float(matchPathlen.group(1))

        matchReachTime=paternReachTime.search(txt, re.MULTILINE)
        if matchReachTime:
            # print "reachability analysis="+matchReachTime.group(1)+"\t",
            result['reach_anal'] = float(matchReachTime.group(1))

        matchInclTime=paternInclTime.search(txt, re.MULTILINE)
        if matchInclTime:
            # print "inclusion check="+matchInclTime.group(1)+"\t",
            result['inclusion'] = float(matchInclTime.group(1))

        matchPaternIters = paternIters.search(txt, re.MULTILINE)
        if matchPaternIters:
            # print "iterations="+matchPaternIters.group(1)+"\t",
            result['iterations'] = float(matchPaternIters.group(1))

        matchPaternTrans = paternTrans.search(txt, re.MULTILINE)
        if matchPaternTrans:
            # print "trans.="+matchPaternTrans.group(1)+"\t",
            result['trans'] = float(matchPaternTrans.group(1))

        matchPaternETB = paternETB.search(txt, re.MULTILINE)
        if matchPaternETB:
            # print "e2b="+matchPaternETB.group(1)+"\t",
            result['e2b'] = float(matchPaternETB.group(1))


        matchSimParams = paternSimParams.findall(txt)
#        if len(matchSimParams) >=2:
            # print matchSimParams[1]
 #       else:
            # print ""


            
        # matchParams = paternParams.findall(txt)
        # if len(matchParams) >=2:
        #     print matchParams[1]
        # else:
        #     print ""

    if len(matchRes) == 0:
        matchError = paternError.search(txt, re.MULTILINE)
        if matchError:
            print("Error in "+log+":\t"+matchError.group(), file=sys.stderr)
            result['result'] = 'ERROR'
        elif txt.find("Command terminated by signal 9") >= 0:
            result['result'] = 'TO'
        else:
            print("Could not parse outpout in "+log+" (memout?)", file=sys.stderr)
            result['result'] = 'MO'
    return result

# gives modelname -> property -> experiment -> list of results
def getResults(logs):
    logfile = open(logs)
    alldict = {}


    for log in logfile:
        match = paternLog.search(log)
        if match:
            modelName = match.group(1)
            propName = match.group(2)
            expName = match.group(3)
            no = int(match.group(4))
            if modelName in alldict:
                mdict = alldict[modelName]
            else:
                mdict = {}
                alldict[modelName] = mdict

            if propName in mdict:
                pdict = mdict[propName]
            else:
                pdict = {}
                mdict[propName] = pdict

            if expName in pdict:
                edict = pdict[expName]
            else:
                edict = []
                pdict[expName] = edict
            
            res = parselog(log.strip())
            # print(res)
            edict.append(res)

        else:
            print("Could not match log "+log)
            exit(1)
        

    logfile.close()
    return alldict

def getSummary(results):
    summary = {}
    for model in results:
        mres = results[model]
        msum = {}
        summary[model] = msum

        for prop in mres:
            pres = mres[prop]
            psum = {}
            msum[prop] = psum

            for exp in pres:
                eres = pres[exp]
                esum = {}
                psum[exp] = esum

                for res in eres:
                    if 'result' not in res:
                        print("ERROR: Couldn't compute summary", file=sys.stderr)
                        exit(1)

                    if res['result'] == 'TO' or res['result'] == 'ERROR' or res['result'] == 'MO':
                        esum['result'] = res['result']
                        break

                    if 'result' not in esum:
                        esum['result'] = res['result']
                    elif esum['result'] != res['result']:
                        esum['result'] = 'MISMATCH'

                    for key in ["time", "path_len_avg", "reach_anal", "inclusion", "iterations","e2b","trans"]:
                        if (key in res):
                            if key not in esum:
                                esum[key] = res[key]/len(eres)
                            else:
                                esum[key] += res[key]/len(eres)

                # print model+", "+prop+", "+exp+" -> ",
                # print esum
    return summary
                



def printTable(summary, experiments,metadata=True,iterations=True,path_len=False):
    for model in sorted(summary.keys()):
        if model in model_data:
            print("{0:40s} ".format(model_data[model][0]), end=' ') 
        else:
            print("{0:40s} ".format(model), end=' ') 

        if metadata:
            if model in model_data:
                print("& {0:20s} ".format(model_data[model][1]), end=' ') 
                print("& {0:20s} ".format(model_data[model][2]), end=' ') 
                print("& {0:20s} ".format(model_data[model][3]), end=' ') 
            else:
                print("& {0:20s} ".format('?'), end=' ') 
                print("& {0:20s} ".format('?'), end=' ') 
                print("& {0:20s} ".format('?'), end=' ') 


        properties = list(summary[model].keys())
        if len(properties) != 1:
            print("Don't have an unique property for "+model, file=sys.stderr)
            exit(1)
        else:
            prop = properties[0]

        psummary = summary[model][prop]

        # find the best experimennt for the model
        best = ""
        best_time = 10000.0
        allbest = ""
        allbest_time = 10000.0
        for exp in experiments:
            if exp in psummary:
                esummary = psummary[exp]
                res = esummary['result']
                if exp.find("simcnd")>=0 or exp.find("simrnd") >= 0:
                    if res == "true" or res == "false":
                        if esummary['time'] < best_time:
                            best_time = esummary['time']
                            best = exp

                if res == "true" or res == "false":
                    if esummary['time'] < allbest_time:
                        allbest_time = esummary['time']
                        allbest = exp
                

        for exp in experiments:
            if exp not in psummary:
                print("&"+" "*19, end=' ')
            else:
                esummary = psummary[exp]
                res = esummary['result']
                if res == "ERROR":
                    print("& $\\error$   ", end=' ')
                elif res == "TO":
                    print("& $\\tout$    ", end=' ')
                elif res == "MO":
                    print("& $\\memout$  ", end=' ')
                else:
                    if model in correct_result and prop in correct_result[model]:
                        cres = correct_result[model][prop]
                        if res == 'true' and not cres:
                            print("& $\\wr$      ", end=' ')
                        elif res == 'false' and cres:
                            print("& $\\wr$      ", end=' ')
                        elif exp == best:
                            print("& \\ba{"+"{0:5.1f}s".format(esummary['time'])+"}", end=' ')
                        elif exp == allbest:
                            print("& \\bb{"+"{0:5.1f}s".format(esummary['time'])+"}", end=' ')
                        else:
                            print("& {0:7.1f}s   ".format(esummary['time']), end=' ')
                    else:
                        print("Could not check correctnes of "+model+", "+prop, file=sys.stderr)
                        print("& {0:7.1f}s   ".format(esummary['time']), end=' ')


                if iterations and exp.find("sim") >= 0:
                    if 'iterations' in esummary:
                        print("& {0:5.0f}    ".format(esummary['iterations']), end=' ')
                    else:
                        print("&   -      ", end=' ')

                if path_len and exp.find("sim") >= 0:
                    if 'path_len_avg' in esummary:
                        print("& {0:5.0f}    ".format(esummary['path_len_avg']), end=' ')
                    else:
                        print("&   -      ", end=' ')

                
                if exp.find("simreach") >= 0:
                    if 'reach_anal' in esummary:
                        print("& {0:5.1f}s   ".format(esummary['reach_anal']), end=' ')
                    else:
                        print("&   -      ", end=' ')

        print("  \\\\")


def printTableCSV(summary, experiments,metadata=True,iterations=True,path_len=False):
    for model in sorted(summary.keys()):
        if model in model_data:
            print("{0:40s}".format(model_data[model][0]), end=' ')
        else:
            print("{0:40s}".format(model), end=' ')

        if metadata:
            if model in model_data:
                print(";{0:20s}".format(model_data[model][1]), end=' ') 
                print(";{0:20s}".format(model_data[model][2]), end=' ') 
                print(";{0:20s}".format(model_data[model][3]), end=' ') 
            else:
                print(";{0:20s}".format('?'), end=' ') 
                print(";{0:20s}".format('?'), end=' ') 
                print(";{0:20s}".format('?'), end=' ') 


        properties = list(summary[model].keys())
        if len(properties) != 1:
            print("Don't have an unique property for "+model, file=sys.stderr)
            exit(1)
        else:
            prop = properties[0]

        psummary = summary[model][prop]

        # find the best experimennt for the model
        best = ""
        best_time = 10000.0
        allbest = ""
        allbest_time = 10000.0
        for exp in experiments:
            if exp in psummary:
                esummary = psummary[exp]
                res = esummary['result']
                if exp.find("simcnd")>=0 or exp.find("simrnd") >= 0:
                    if res == "true" or res == "false":
                        if esummary['time'] < best_time:
                            best_time = esummary['time']
                            best = exp

                if res == "true" or res == "false":
                    if esummary['time'] < allbest_time:
                        allbest_time = esummary['time']
                        allbest = exp
                

        for exp in experiments:
            if exp not in psummary:
                print(";"+" "*19, end=' ')
            else:
                esummary = psummary[exp]
                res = esummary['result']
                if res == "ERROR":
                    print(";error", end=' ')
                elif res == "TO":
                    print(";TO", end=' ')
                elif res == "MO":
                    print(";MO", end=' ')
                else:
                    if model in correct_result and prop in correct_result[model]:
                        cres = correct_result[model][prop]
                        if res == 'true' and not cres:
                            print(";wr", end=' ')
                        elif res == 'false' and cres:
                            print(";wr", end=' ')
                        elif exp == best:
                            print(";?{0:5.1f}".format(esummary['time']), end=' ')
                        elif exp == allbest:
                            print(";??{0:5.1f}".format(esummary['time']), end=' ')
                        else:
                            print(";{0:7.1f}".format(esummary['time']), end=' ')
                    else:
                        print("Could not check correctnes of "+model+", "+prop, file=sys.stderr)
                        print("& {0:7.1f}s   ".format(esummary['time']), end=' ')


                if iterations and exp.find("sim") >= 0:
                    if 'iterations' in esummary:
                        print(";{0:5.0f}".format(esummary['iterations']), end=' ')
                    else:
                        print(";", end=' ')

                if path_len and exp.find("sim") >= 0:
                    if 'path_len_avg' in esummary:
                        print(";{0:5.0f}".format(esummary['path_len_avg']), end=' ')
                    else:
                        print(";-", end=' ')

                
                if exp.find("simreach") >= 0:
                    if 'reach_anal' in esummary:
                        print(";{0:5.1f}".format(esummary['reach_anal']), end=' ')
                    else:
                        print(";-", end=' ')

        print("")

def parse(logs):
    results = getResults(logs)
    summary = getSummary(results)
    print(summary)
    #printTable(summary,["simcnd","simrnd_1e_4","simreach","mc"])
    printTableCSV(summary,["simcnd","gsimcnd"],False,True,True) 
    #printTable(summary,["simcnd_cb1"],False,False,True)
    # printTable(summary,["simcnd_ltl","mc_ltl","simcnd_rew","mc_rew"],False,False)


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
                generateScript(model, prop, const, minprob, arguments, num_arg, scriptname, logname)
                print(scriptname)

    



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
    parser.add_argument('-g','--generate',action='store_true',help='generate test scripts')
    parser.add_argument('-p','--parse',action='store_true',help='parse test logs')
    parser.add_argument('-e','--experiments',nargs=1,help="experiments file")
    parser.add_argument('-b','--benchmarks',nargs=1,help="benchmarks file")
    parser.add_argument('-l','--logs',nargs=1,help="logs file")
    parser.add_argument('-r','--repeat',type=int,default=1,help="how many repeats")
    args = parser.parse_args()           

    if ((not args.generate) and (not args.parse)):
        parser.print_help()
        exit(1)
    if (args.generate):
        if (args.benchmarks == None or args.experiments == None):
            parser.print_help()
            exit(1)
        generate(args.benchmarks[0], args.experiments[0],args.repeat)
    else:
        if (args.logs == None):
            parser.print_help()
            exit(1)
        parse(args.logs[0])
