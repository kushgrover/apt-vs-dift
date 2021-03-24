
PRISM="/home/maxi/masterarbeit/games-brtdp/code/prism_4_4/prism-4.4.beta-linux64/bin/prism"
HEUR="/home/maxi/masterarbeit/games-brtdp/code/maxi/bin/prism"
STORM="/home/maxi/masterarbeit/games-brtdp/code/storm/storm-stable/build/bin/storm"
PRISMG="/home/maxi/masterarbeit/games-brtdp/code/prism-games-2.0.beta3-src/bin/prism"

##################### MDPs ######################

#Firewire
firewireM="../../examples/models/firewire/impl/deadline.nm"
firewireP="../../examples/models/firewire/impl/deadline.pctl"
firewireLog="../logs/data$1/firewire"
mkdir $firewireLog

#impl
deadline.nm deadline.pctl -const deadline=$deadline,delay=3,fast=0.5 -s


wlanM="../../examples/models/wlan/benchmarks/wlan"
wlanP="../../examples/models/wlan/benchmarks/wlan.pctl"
wlanLog="../logs/data$1/wlan"
mkdir $wlanLog

$PRISM wlan0_time_bounded.nm wlan_time_bounded.pctl -const TRANS_TIME_MAX=10,DEADLINE=100 -nopre -m -prop 1
prism examples/models/wlan/benchmarks/wlan2.nm examples/models/wlan/benchmarks/wlan.pctl -const COL=2


zeroconfM="../../examples/models/zeroconf/zeroconf.nm"
zeroconfP="../../examples/models/zeroconf/zeroconf.pctl"
zeroconfLog="../logs/data$1/zeroconf"
mkdir $zeroconfLog

merM="../../examples/models/mer/mer2_2_full.nm"
merP="../../examples/models/mer/mer.pctl"
merLog="../logs/data$1/mer"
mkdir $merLog



#prism csma2_2.nm csma.pctl -const k=1 -s (bis 3_6; prop one of the 10, 5 seems good)

#leader/asynchronous/leader6.nm leader/asynchrous/leader.pctl -prop 1 (up til 6, that wont work for BRTDP)

#SMGs
mdsmM="../../examples/models/mdsm/mdsm.prism"
mdsmP="../../examples/models/mdsm/mdsm.props"
mdsmLog="../logs/data$1/mdsm"
mkdir $mdsmLog

teamformM="../../examples/models/team-form/team-form-"
teamformP="../../examples/models/team-form/team-form.props" 
teamformLog="../logs/data$1/teamform"
mkdir $teamformLog

cdmsnM="../../examples/models/cdmsn/cdmsn.prism"
cdmsnP="../../examples/models/cdmsn/cdmsn.props"
cdmsnLog="../logs/data$1/cdmsn"
mkdir $cdmsnLog prop 1

cloudM="../../examples/models/cloud/cloud_"
cloudP="../../examples/models/cloud/cloud.props" 
cloudLog="../logs/data$1/cloud"
mkdir $cloudLog

for modelLog in $firewireLog $wlanLog $zeroconfLog $merLog; do
    for solver in "PRISM_default" "PRISM_interval" "PRISM_large" "Storm" "RTDP_5" "RTDP_6" "VIA" "VIC"; do
        mkdir $modelLog/$solver
    done
done

for modelLog in $mdsmLog $teamformLog $cdmsnLog $cloudLog; do
    for solver in "PRISM-games" "RTDP_5" "RTDP_6" "RTDP_69" "RTDP_70" "VIA" "VIC"; do
        mkdir $modelLog/$solver
    done
done


function runBench { #outputfile,timeout,solver,model,props
    #echo "perf stat -o $1.stat timeout $2 $3 $4 $5 > $1.log"
    perf stat -o $1.stat timeout $2 $3 $4 $5 > $1.log
}



echo "firewire"
LOG=$firewireLog
for d in ; do #220 240 260 280; do
    echo $d
    ####CAREFUL: DEADLINE AND DELAY ARE SWITCHED IN OUTPUT
    #PRISM and Storm are fine
    echo "Storm"
    HEURLOG="$LOG/Storm"
    #$(runBench $HEURLOG/\_deadline\_36\_delay\_$d 15m "$STORM --prism" $firewireM "--prop $firewireP --constants deadline=$d,delay=36,fast=0.5 --timemem --cudd:maxmem 7000" )


    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for rep in ; do #((rep=1; rep<21; rep++)); do
        #echo runBench $HEURLOG/deadline\_36\_delay\_$d\_rep\_$rep 10m $HEUR $firewireM "$firewireP -ex -const deadline=$d_d,delay=36,fast=0.5 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose"
        $(runBench $HEURLOG/deadline\_36\_delay\_$d\_rep\_$rep 10m $HEUR $firewireM "$firewireP -ex -const deadline=$d,delay=36,fast=0.5 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for rep in ; do #((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/deadline\_36\_delay\_$d\_rep\_$rep 10m $HEUR $firewireM "$firewireP -ex -const deadline=$d,delay=36,fast=0.5 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done

    #VIA and VIC are fine
done


echo "wlan"
LOG=$wlanLog
for k in ; do # 2 3 4 5 6; do
for COL in 2 3 4 5 6; do
    echo "$k,$COL"
    
    echo "PRISM_default"
    HEURLOG="$LOG/PRISM_default"
	$(runBench $HEURLOG/k\_$k\_COL\_$COL 10m $PRISM "$wlanM$k.nm" "$wlanP -ex -const COL=$COL -javamaxmem 4g")

    echo "PRISM_interval"    
    HEURLOG="$LOG/PRISM_interval"
	$(runBench $HEURLOG/k\_$k\_COL\_$COL 10m $PRISM "$wlanM$k.nm" "$wlanP -ex -const COL=$COL -ii -javamaxmem 4g")

    echo "PRISM_large"    
    HEURLOG="$LOG/PRISM_large"
	$(runBench $HEURLOG/k\_$k\_COL\_$COL 10m $PRISM "$wlanM$k.nm" "$wlanP -ex -const COL=$COL -m -jor -javamaxmem 4g")

	echo "Storm"
    HEURLOG="$LOG/Storm"
    $(runBench $HEURLOG/k\_$k\_COL\_$COL 10m "$STORM --prism" "$wlanM$k.nm" "--prop $wlanP --constants COL=$COL --timemem --cudd:maxmem 7000" )

    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/k\_$k\_COL\_$COL\_rep\_$rep 10m $HEUR "$wlanM$k.nm" "$wlanP -ex -const COL=$COL -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/k\_$k\_COL\_$COL\_rep\_$rep 10m $HEUR "$wlanM$k.nm" "$wlanP -ex -const COL=$COL -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done

    echo "VIA"
    HEURLOG="$LOG/VIA"
    perf stat -o $HEURLOG/k\_$k\_COL\_$COL.stat timeout 10m $HEUR $wlanM$k.nm $wlanP -ex -const COL=$COL -prop 1 -BVI_A > $HEURLOG/k\_$k\_COL\_$COL.log

    echo "VIC"
    HEURLOG="$LOG/VIC"
    perf stat -o $HEURLOG/k\_$k\_COL\_$COL.stat timeout 10m $HEUR $wlanM$k.nm $wlanP -ex -const COL=$COL -prop 1 -BVI_C > $HEURLOG/k\_$k\_COL\_$COL.log
done
done

echo "zeroconf"
LOG=$zeroconfLog
for K in ; do # 2 10; do #4
for N in 20 500 1000; do #1000
    echo "$K,$N"
    
    #PRISM and Storm are fine
    echo "Storm"
    HEURLOG="$LOG/Storm"
    #$(runBench $HEURLOG/K\_$K\_N\_$N 10m "$STORM --prism" $zeroconfM "--prop $zeroconfP --constants reset=false,K=$K,N=$N,err=0.1 --timemem --cudd:maxmem 7000" )


    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/K\_$K\_N\_$N\_rep\_$rep 10m $HEUR $zeroconfM "$zeroconfP -ex -const reset=false,K=$K,N=$N,err=0.1 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/K\_$K\_N\_$N\_rep\_$rep 10m $HEUR $zeroconfM "$zeroconfP -ex -const reset=false,K=$K,N=$N,err=0.1 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done


    echo "VIA"
    HEURLOG="$LOG/VIA"
    $(runBench $HEURLOG/K\_$K\_N\_$N 10m $HEUR $zeroconfM "$zeroconfP -ex -const reset=false,K=$K,N=$N,err=0.1 -BVI_A" )

    echo "VIC"
    HEURLOG="$LOG/VIC"
    $(runBench $HEURLOG/K\_$K\_N\_$N 10m $HEUR $zeroconfM "$zeroconfP -ex -const reset=false,K=$K,N=$N,err=0.1 -BVI_C" )
done
done





echo "mer moved down"



echo "SMGs"


echo "mdsm"
LOG=$mdsmLog
for prop in 1 2; do
    echo $prop

    #Prism is done

    #echo "RTDP_5"
    #HEURLOG="$LOG/RTDP_5"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/prop\_$prop\_rep\_$rep 10m $HEUR $mdsmM "$mdsmP -ex -prop $prop -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    #done

    #echo "RTDP_6"
    #HEURLOG="$LOG/RTDP_6"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/prop\_$prop\_rep\_$rep 10m $HEUR $mdsmM "$mdsmP -ex -prop $prop -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    #done

    echo "RTDP_69"
    HEURLOG="$LOG/RTDP_69"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/thresh\_prop\_$prop\_rep\_$rep 10m $HEUR $mdsmM "$mdsmP -ex -prop $prop -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 69 -heuristic_verbose" )
    done

    echo "RTDP_70"
    HEURLOG="$LOG/RTDP_70"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/thresh\_prop\_$prop\_rep\_$rep 10m $HEUR $mdsmM "$mdsmP -ex -prop $prop -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 70 -heuristic_verbose" )
    done

    #VIA and VIC are done
done


echo "teamform"
LOG=$teamformLog
for N in ; do #4; do #3 4; do
    echo $N

    #Prism is done

    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$teamformM$N.prism" "$teamformP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$teamformM$N.prism" "$teamformP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done

    echo "RTDP_69"
    HEURLOG="$LOG/RTDP_69"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$teamformM$N.prism" "$teamformP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 69 -heuristic_verbose" )
    done

    echo "RTDP_70"
    HEURLOG="$LOG/RTDP_70"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$teamformM$N.prism" "$teamformP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 70 -heuristic_verbose" )
    done

    #VIA and VIC are done
done


echo "cdmsn"
LOG=$cdmsnLog
for n in ; do 
    echo "PRISM-games"
    HEURLOG="$LOG/PRISM-games"
    $(runBench $HEURLOG/A\_1 10m $PRISMG $cdmsnM "$cdmsnP -ex -javamaxmem 4g" )

    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/A\_1\_rep\_$rep 10m $HEUR $cdmsnM "$cdmsnP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/A\_1\_rep\_$rep 10m $HEUR $cdmsnM "$cdmsnP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done

    echo "RTDP_69"
    HEURLOG="$LOG/RTDP_69"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/A\_1\_rep\_$rep 10m $HEUR $cdmsnM "$cdmsnP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 69 -heuristic_verbose" )
    done

    echo "RTDP_70"
    HEURLOG="$LOG/RTDP_70"
    for ((rep=1; rep<21; rep++)); do
        $(runBench $HEURLOG/A\_1\_rep\_$rep 10m $HEUR $cdmsnM "$cdmsnP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 70 -heuristic_verbose" )
    done

    echo "VIA"
    HEURLOG="$LOG/VIA"
    $(runBench $HEURLOG/A\_1 10m $HEUR $cdmsnM "$cdmsnP -ex -BVI_A" )

    echo "VIC"
    HEURLOG="$LOG/VIC"
    $(runBench $HEURLOG/A\_1 10m $HEUR $cdmsnM "$cdmsnP -ex -BVI_C" )
done

echo "cloud"
LOG=$cloudLog
for N in ; do #5 6; do
    echo $N

    echo "PRISM-games"
    HEURLOG="$LOG/PRISM-games"
    #$(runBench $HEURLOG/N\_$N 10m $PRISMG "$cloudM$N.prism" "$cloudP -ex -javamaxmem 4g" )

    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    #done

    #echo "RTDP_6"
    #HEURLOG="$LOG/RTDP_6"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    #done

    #echo "RTDP_69"
    #HEURLOG="$LOG/RTDP_69"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 69 -heuristic_verbose" )
    #done

    #echo "RTDP_70"
    #HEURLOG="$LOG/RTDP_70"
    #for ((rep=1; rep<21; rep++)); do
    #    $(runBench $HEURLOG/N\_$N\_rep\_$rep 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 70 -heuristic_verbose" )
    #done

    echo "VIA"
    HEURLOG="$LOG/VIA"
    $(runBench $HEURLOG/N\_$N 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -BVI_A" )

    echo "VIC"
    HEURLOG="$LOG/VIC"
    $(runBench $HEURLOG/N\_$N 10m $HEUR "$cloudM$N.prism" "$cloudP -ex -BVI_C" )
done



echo "mer"
LOG=$merLog
for x in ; do # 0.9999 0.0001; do
for N in 4500; do #1500 3000 4500; do
    echo "$x,$N"
    if [[ x == 0.9999 ]] && [[ N == 1500 ]]
        then
        continue;
    fi
    
    echo "PRISM_default"
    HEURLOG="$LOG/PRISM_default"
#    $(runBench $HEURLOG/x\_$x\_n\_$N 15m $PRISM $merM "$merP -ex -const n=$N,x=$x,K=30 -javamaxmem 7g" )

    echo "PRISM_interval" 
    HEURLOG="$LOG/PRISM_interval"
#    $(runBench $HEURLOG/x\_$x\_n\_$N 15m $PRISM $merM "$merP -ex -const n=$N,x=$x,K=30 -ii -javamaxmem 7g" )

    echo "PRISM_large" 
    HEURLOG="$LOG/PRISM_large"
#    $(runBench $HEURLOG/x\_$x\_n\_$N 15m $PRISM $merM "$merP -ex -const n=$N,x=$x,K=30 -m -jor -javamaxmem 7g" )

    echo "Storm"
    HEURLOG="$LOG/Storm"
    $(runBench $HEURLOG/x\_$x\_n\_$N 15m "$STORM --prism" $merM "--prop $merP --constants n=$N,x=$x,K=30 --timemem --cudd:maxmem 7000" )

    echo "RTDP_5"
    HEURLOG="$LOG/RTDP_5"
    for ((rep=1; rep<2; rep++)); do
        $(runBench $HEURLOG/x\_$x\_n\_$N\_rep\_$rep 15m $HEUR $merMS "$merPS -ex -const n=$N,x=$x,K=30 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 5 -heuristic_verbose" )
    done

    echo "RTDP_6"
    HEURLOG="$LOG/RTDP_6"
    for ((rep=1; rep<2; rep++)); do
        $(runBench $HEURLOG/x\_$x\_n\_$N\_rep\_$rep 15m $HEUR $merMS "$merPS -ex -const n=$N,x=$x,K=30 -heuristic RTDP_ADJ -next_state MAX_DIFF -RTDP_ADJ_OPTS 6 -heuristic_verbose" )
    done

    #VIA and VIC are done

done
done





