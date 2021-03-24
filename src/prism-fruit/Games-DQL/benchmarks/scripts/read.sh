#!/bin/bash


#replace with Elapsed (wall clock) time; then use another field
function read_time {
    if [[ $4 == "none" ]] && [[ $6 == "none" ]]
    then
        line=$(grep -s "User time (seconds):" $1/$2\_$3.stat);
        result=$(grep -s "Result" $1/$2\_$3.log);
    elif [[ $4 == "none" ]]
    then 
        line=$(grep -s "User time (seconds):" $1/$2\_$3-$6.stat);
        result=$(grep -s "Result" $1/$2\_$3-$6.log);
    elif [[ $6 == "none" ]]
    then
        line=$(grep -s "User time (seconds):" $1/$2\_$3\_$4\_$5.stat);
        result=$(grep -s "Result" $1/$2\_$3\_$4\_$5.log);
    else 
        line=$(grep -s "User time (seconds):" $1/$2\_$3\_$4\_$5-$6.stat);
        result=$(grep -s "Result" $1/$2\_$3\_$4\_$5-$6.log);
    fi

    if [ $? -eq 1 ]
    then
        echo "X"; #no result
    elif [[ $line == "" ]]
    then   
        echo "I"; #inexistant
    else
        
        echo $(echo "scale=2;$(echo $line | cut -d ':' -f 2)/1" | bc)
    fi
}  



##################################### START

BASEPATH="../logs/data0"

OUTPUT="../../../../paper/tables/time.tex"

solvers=( PRISM_default PRISM_interval PRISM_large PRISM-games Storm BVIA RTDP_1 RTDP_65 ATVA )
models=( firewire wlan zeroconf csma leader mer mdsm teamform cdmsn cloud )

columns="|c|c|"
header="Model & Parameters"
num=2
for s in ${solvers[@]}; do
    columns=$columns"c|"
    header="$header & $s"
    num=$(echo "$num+1" | bc)
done

str="\begin{table}[]
\centering
\small
\caption{Experimental results Time}
\label{table:time}
\makebox[\linewidth]{
\begin{tabular}{$columns}
\hline
$header \\\\ \hline \hline"
echo $str > $OUTPUT


for m in ${models[@]}; do
    echo $m
    
    if [[ $m == "firewire" ]] ; then
        par1=( 220 240 260 280 )
        parName1="deadline"
        par2=( 36 )
        parName2="delay"
        rows=4
    elif [[ $m == "wlan" ]] ; then
        par1=( 2 4 6) #( 2 3 4 5 6 )
        parName1="k"
        par2=(2 6) #( 2 3 4 5 6 )
        parName2="COL"
        rows=6
    elif [[ $m == "zeroconf" ]] ; then
        par1=( 2 10 )
        parName1="K"
        par2=( 20 500 1000)
        parName2="N"
        rows=6
    elif [[ $m == "csma" ]] ; then
        par1=( 2 3 )
        parName1="N"
        par2=( 2 4 6 )
        parName2="K"
        rows=6
    elif [[ $m == "leader" ]] ; then
        par1=(3 4 5 6 )
        parName1=N
        par2=( 1 )
        parName2="none"
        rows=4
    elif [[ $m == "mer" ]]  ; then
        par1=( 0.0001 0.1 0.5 )
        parName1="x"
        par2=( 1500 3000 ) #( 1500 3000 4500 )
        parName2="n"
        rows=6
    elif [[ $m == "mdsm" ]] ; then
        par1=( 1 2 )
        parName1="prop"
        par2=( 1 )
        parName2="none"
        rows=2
    elif [[ $m == "teamform" ]] ; then
        par1=( 3 4 )
        parName1="N"
        par2=( 1 )
        parName2="none"
        rows=2
    elif [[ $m == "cdmsn" ]] ; then
        par1=( 1 )
        parName1="A"
        par2=( 1 )
        parName2="none"
        rows=1
    elif [[ $m == "cloud" ]] ; then
        par1=( 5 6 )
        parName1="N"
        par2=( 1 )
        parName2="none"
        rows=2
    else
        echo "I don't know this model: $m Exiting."
        exit
    fi
     
    echo "\multirow{$rows}{*}{$m} " >> $OUTPUT


    for p1 in ${par1[@]}; do
    for p2 in ${par2[@]}; do
            if [[ $parName2 == "none" ]]
            then
                echo "& $parName1=$p1 " >> $OUTPUT
            else
                echo "& $parName1=$p1 $parName2=$p2" >> $OUTPUT
            fi
            echo "& $parName1=$p1 $parName2=$p2"

            for solver in ${solvers[@]}; do
                CURRPATH="$BASEPATH/$solver/$m"
    
                if [[ $solver == "RTDP"* ]] || [[ $solver == "ATVA" ]]; then
                        sumTime=0
                        timeouts=0
                        samples=0
                        for ((rep=1; rep<21; rep++)); do
                            time=$(read_time $CURRPATH $parName1 $p1 $parName2 $p2 $rep)
                            if [[ $time == "X" ]]
                            then
                                timeouts=$(echo "$timeouts+1" | bc)
                            elif [[ $time == "I" ]]
                            then
                                time="File inexistant, do nothing"
                            else
                                samples=$(echo "$samples+1" | bc)
                                sumTime=$(echo "$sumTime+$time" | bc)
                            fi
                        done
                        if [[ $samples == 0 ]]
                        then
                            echo "& X ($samples,$timeouts) " >> $OUTPUT
                        else
                            time=$(echo "scale=2;$sumTime/$samples" | bc)
                            echo "& $time ($samples,$timeouts) " >> $OUTPUT
                        fi
                else
                    time=$(read_time $CURRPATH $parName1 $p1 $parName2 $p2 "none")
                    echo "& $time " >> $OUTPUT
                fi
            done 

            echo "\\\\ \cline{2-$num}" >> $OUTPUT
    done
    done
    echo "\hline" >> $OUTPUT
done



########################################## END
str="\end{tabular}
}
\end{table}"
echo $str >> $OUTPUT
