#!/bin/bash

if [[ -z "$1" || -z "$2" ]]; then
	echo "$0 <GENERATED_CONFIG_FILE> <NUM_REPETITIONS>";
	exit 1
else
	if [ -f $1 ]; then
	   RUNCONFIG=$1
	else
	   echo "File $1 does not exist"
	   echo "Run $0 <GENERATED_CONFIG_FILE> <NUM_REPETITIONS> with a valid file containing generated run configurations";
	   exit 1;
	fi
	
fi

PRISM="../bin/prism"
MAINLOG="logs"
REPS=$2
now="$( date +"%Y-%m-%d-%H-%M-%S" )"
LOGDIR="$MAINLOG/$now"

mkdir -p $LOGDIR

printf "Saving logs to %s\n\n" $LOGDIR

# For each line in file
grep -v '^#' $RUNCONFIG | while IFS= read -r line; do
    model=$(echo $line | gawk 'match($0, /models\/([^\ ]+)/, a) { print a[1] }')
    property=$(echo $line | gawk 'match($0, /\-prop(\ )+([^\ ]+)/, a) { print a[2] }')
    consts=$(echo $line | gawk 'match($0, /\-const(\ )+([^\ ]+)/, a) { print a[2] }')
    colour=$(echo $line | gawk 'match($0, /\-RTDP_ADJ_OPTS(\ )+([0-9])/, a) { print a[2] }')
    case $colour in 
    	0)
			colour="white"
			;;
		1)
			colour="grey"
			;;
		2)
			colour="black"
			;;
	esac
    fileprefix="$model-$property-$consts-$colour"
    for ((r=1;r<=$REPS; r++)); do
		filename="$fileprefix-$r"
		printf 'Running (%s/%s) %s\nConfig: %s\n' "$r" "$REPS" "$filename" "$line"
		pueue add -- /usr/bin/time -v -o $LOGDIR/$filename.stat timeout 60m $PRISM $line '2>&1' '|' tee $LOGDIR/$filename.log
		printf 'Output and stats in\n\t- %s.stat\n\t- %s.log\n\n' "$LOGDIR/$filename" "$LOGDIR/$filename"

	done;
done;

printf "Run 'pueue status' to see progress. Logs would become available in %s.\n\n" $LOGDIR

