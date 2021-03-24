#!/bin/bash

#Solver names
PRISM="/root/CECsForCAV/games-brtdp/code/prism_4_4/prism-4.4.beta-linux64/bin/prism"
HEUR="/root/CECsForCAV/games-brtdp/code/maxi/bin/prism"
STORM="/root/CECsForCAV/games-brtdp/code/storm/storm-stable/build/bin/storm"
PRISMG="/root/CECsForCAV/games-brtdp/code/prism-games-2.0.beta3-linux64/bin/prism"
ATVA="/root/CECsForCAV/games-brtdp/code/atva/prism-mdp-learning/bin/prism"

mainLog="/root/CECsForCAV/games-brtdp/code/maxi/benchmarks/logs/data$1"
mkdir $mainLog

for solverName in PGex PRISM-games BVIA; do #PRISM_default PRISM_explicit PRISM_interval PRISM_large PRISM-games PGex Storm BVIA RTDP_0 RTDP_1 RTDP_5 RTDP_65 RTDP_69 ATVA; do 
    echo $(date +%d.%m.%y---%H:%M:%S) $solverName
    logfile="$mainLog/$solverName"
    mkdir $logfile

    #solver defs (solver, 2stormopts, logfile and options, const(=-const or -constants), canDoGames)
    if [[ $solverName == "PRISM_default" ]] ; then
        solver=$PRISM
        options="-javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=false
        reps=0
    elif [[ $solverName == "PRISM_explicit" ]] ; then
        solver=$PRISM
        options="-ex -javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=false
        reps=0
    elif [[ $solverName == "PRISM_interval" ]] ; then
        solver=$PRISM
        options="-ii -ex -ddextraactionvars 100 -javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=false
        reps=0
    elif [[ $solverName == "PRISM_large" ]] ; then
        solver=$PRISM
        options="-m -jor -javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=false
        reps=0
    elif [[ $solverName == "PRISM-games" ]] ; then
        solver=$PRISMG
        options="-javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=true
        reps=0
    elif [[ $solverName == "PGex" ]] ; then
        solver=$PRISMG
        options="-ex -javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=true
        reps=0
    elif [[ $solverName == "Storm" ]] ; then
        solver=$STORM
        options="--verbose --timemem --sylvan:maxmem 8192 --cudd:maxmem 8192"
        #Storm stuff
        StormPrism="--prism"
        StormProp="--prop"
        const="--constants"
        canDoGames=false
        reps=0
    elif [[ $solverName == "BVIA" ]] ; then
        solver=$HEUR
        options="-BVI_A -javamaxmem 16g"
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=true
        reps=0
    elif [[ $solverName == "RTDP"* ]] ; then
        solver=$HEUR
        options="-heuristic RTDP_ADJ -RTDP_ADJ_OPTS $(echo $solverName | cut -d '_' -f 2) -javamaxmem 16g -heuristic_verbose" #CAREFUL: NOT USED; OVERWRITTEN LATER FOR MDP OPT
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=true
        reps=1
    elif [[ $solverName == "ATVA" ]] ; then
        solver=$ATVA
        options="-heuristic RTDP_UNBOUNDED -next_state MAX_DIFF -heuristic_verbose -heuristic_epsilon 1e-6 -ex" #unable to set javamaxmem
        #Storm stuff
        StormPrism=""
        StormProp=""
        const="-const"
        canDoGames=false
        reps=2
    else
        echo "I don't know this solver: $solverName Exiting."
        exit
    fi

  

        




    for m in mdsm teamform cdmsn cloud; do # mdsm teamform cdmsn cloud; do

        echo $(date +%d.%m.%y---%H:%M:%S) $m

        #model defs (model, props, params,isMDP)
        
        if [[ $m == "firewire" ]] ; then
            par1=( 220 240 260 280 )
            parName1="deadline"
            par2=( 36 )
            parName2="delay"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/firewire/impl/deadline.nm"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/firewire/impl/deadline.pctl"
            isMDP=true
        elif [[ $m == "wlan" ]] ; then
            par1=( 2 4 6) #( 2 3 4 5 6 )
            parName1="k"
            par2=(2 6) #( 2 3 4 5 6 )
            parName2="COL"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/wlan/benchmarks/wlan"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/wlan/benchmarks/wlan.pctl"
            isMDP=true
        elif [[ $m == "zeroconf" ]] ; then
            par1=( 2 10 )
            parName1="K"
            par2=( 20 500 1000)
            parName2="N"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/zeroconf/zeroconf.nm"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/zeroconf/zeroconf.pctl"
            isMDP=true
        elif [[ $m == "csma" ]] ; then
            par1=( 2 3 )
            parName1="N"
            par2=( 2 4 6 )
            parName2="K"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/csma/csma"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/csma/csma.pctl"
            isMDP=true
        elif [[ $m == "leader" ]] ; then
            par1=(3 4 5 6 )
            parName1=N
            par2=( 1 )
            parName2="none"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/leader/asynchronous/leader"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/leader/asynchronous/leader.pctl"
            isMDP=true
        elif [[ $m == "mer" ]]  ; then
            par1=( 0.0001 0.1 0.5 )
            parName1="x"
            par2=( 1500 3000 ) #( 1500 3000 4500 )
            parName2="n"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/mer/mer2_2_full.nm"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/mer/mer.pctl"
            isMDP=true
        elif [[ $m == "mdsm" ]] ; then
            par1=( 1 2 )
            parName1="prop"
            par2=( 1 )
            parName2="none"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/mdsm/mdsm.prism"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/mdsm/mdsm.props"
            isMDP=false
        elif [[ $m == "teamform" ]] ; then
            par1=( 3 4 )
            parName1="N"
            par2=( 1 )
            parName2="none"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/team-form/team-form-"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/team-form/team-form.props"
            isMDP=false
        elif [[ $m == "cdmsn" ]] ; then
            par1=( 1 )
            parName1="A"
            par2=( 1 )
            parName2="none"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/cdmsn/cdmsn.prism"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/cdmsn/cdmsn.props"
            isMDP=false
        elif [[ $m == "cloud" ]] ; then
            par1=( 5 6 )
            parName1="N"
            par2=( 1 )
            parName2="none"
            model="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/cloud/cloud_"
            props="/root/CECsForCAV/games-brtdp/code/maxi/examples/models/cloud/cloud.props"
            isMDP=false
        else
            echo "I don't know this model: $m Exiting."
            exit
        fi

        if ! $isMDP && ! $canDoGames ; then
            continue;
        fi

        mkdir $logfile/$m

        for p1 in ${par1[@]}; do
        for p2 in ${par2[@]}; do

	    echo $(date +%d.%m.%y---%H:%M:%S) $p1,$p2

            if [[ $parName2 == "none" ]] ; then
                outputfile=$logfile/$m/$parName1\_$p1
            else
                outputfile=$logfile/$m/$parName1\_$p1\_$parName2\_$p2
            fi

            #modelStr and params for each model 

            if [[ $m == "firewire" ]] ; then
                params="$const $parName1=$p1,$parName2=$p2,fast=0.5"
                modelStr=""
            elif [[ $m == "wlan" ]] ; then
                params="$const $parName2=$p2"
                modelStr="$p1.nm"
            elif [[ $m == "zeroconf" ]] ; then
                params="$const $parName1=$p1,$parName2=$p2,reset=false,err=0.1"
                modelStr=""
            elif [[ $m == "csma" ]] ; then
                params="$const k=1"
                modelStr=$p1\_$p2.nm
            elif [[ $m == "leader" ]] ; then
                params=""
                modelStr="$p1.nm"
            elif [[ $m == "mer" ]]  ; then
                params="$const $parName1=$p1,$parName2=$p2,K=30"
                modelStr=""
            elif [[ $m == "mdsm" ]] ; then
                params="-prop $p1"
                modelStr=""
            elif [[ $m == "teamform" ]] ; then
                params=""
                modelStr="$p1.prism"
            elif [[ $m == "cdmsn" ]] ; then
                params=""
                modelStr=""
            elif [[ $m == "cloud" ]] ; then
                params=""
                modelStr="$p1.prism"
            fi
            
            #optimize RTDP for MDPs
            if [[ $solverName == "RTDP"* ]] ; then
                if $isMDP ; then
                    options="-heuristic RTDP_ADJ -RTDP_ADJ_OPTS $(echo "$(echo $solverName | cut -d '_' -f 2)+128" | bc ) -heuristic_verbose -javamaxmem 16g"
                else
                    options="-heuristic RTDP_ADJ -RTDP_ADJ_OPTS $(echo $solverName | cut -d '_' -f 2) -heuristic_verbose -javamaxmem 16g"
                fi
            fi

            if [[ $reps == 0 ]] ; then
                taskset -c 3 /usr/bin/time -v -o $outputfile.stat timeout 15m $solver $StormPrism $model$modelStr $StormProp $props $params $options > $outputfile.log
            else
                for ((r=1;r <= $reps; r++)); do
                    outputfileReps=$outputfile\-$r
                    taskset -c 4 /usr/bin/time -v -o $outputfileReps.stat timeout 15m $solver $StormPrism $model$modelStr $StormProp $props $params $options > $outputfileReps.log
                done
            fi
        done
        done
    done
done

