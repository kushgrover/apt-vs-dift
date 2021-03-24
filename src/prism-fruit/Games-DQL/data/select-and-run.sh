#!/bin/bash

# Use this variable to change the timeout (currently set to 30 minutes)
# Make sure you empty logs/evaluation folder if you change the timeout
# Also make sure no other jobs are writing to logs/evaluation (run 'pueue reset')
TIMEOUT=30m

PRISM="../bin/prism"
MAINLOG="logs"
# now="$( date +"%Y-%m-%d-%H-%M-%S" )"
# LOGDIR="$MAINLOG/$now"
LOGDIR="$MAINLOG/evaluation"

# Start the parallel queueing daemon
printf "Checking if pueue is running\n"
pueue status
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
	pueue --daemon
	if [ $? -ne 0 ]; then
		printf "\nCheck if pueue is installed (pip install pueue)."
		exit 1
	fi
fi


declare -A grey_models
grey_models=( 
	["consensus"]="models/consensus.2.prism properties/consensus.props -prop disagree -const K=2 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:272;Av:2;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["csma"]="models/csma.2-2.prism properties/csma.props -prop some_before -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:1038;Av:2;e:1e-8;d:0.01;p:0.25;post:4\"" \
	["firewire"]="models/firewire.true.prism properties/firewire.true.props -prop deadline -const delay=3,deadline=200 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:83153;Av:3;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["ij-3"]="models/ij.3.v1.prism properties/ij.3.v1.props -prop stable -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:7;Av:3;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["ij-10"]="models/ij.10.v1.prism properties/ij.10.v1.props -prop stable -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:1023;Av:10;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["pacman"]="models/pacman.nm properties/pacman.props -prop crsh -const MAXSTEPS=5 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:498;Av:5;e:1e-8;d:0.01;p:0.08;post:3\"" \
	["philosophers-3"]="models/philosophers-mdp.3.v1.prism properties/philosophers-mdp.3.v1.props -prop eat -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:956;Av:6;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["pnueli-zuck-3"]="models/pnueli-zuck.3.v1.prism properties/pnueli-zuck.v1.props -prop live -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:2701;Av:6;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["rabin-3"]="models/rabin.3.prism properties/rabin.3.props -prop live -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:27766;Av:3;e:1e-8;d:0.01;p:0.03125;post:6\"" \
	["wlan-0"]="models/wlan.0.prism properties/wlan.props -prop sent -const COL=0 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:2954;Av:3;e:1e-8;d:0.01;p:0.0625;post:16\"" \
	["zeroconf"]="models/zeroconf.prism properties/zeroconf.props -prop correct_max -const N=20,K=2,reset=true -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:670;Av:3;e:1e-8;d:0.01;p:1.025262467191601E-4;post:6\"" \
	["cdmsn"]="models/cdmsn.prism properties/cdmsn.props -prop all_prefer_one -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:1240;Av:2;e:1e-8;d:0.01;p:0.059071881973375796;post:5\"" \
	["cloud-5"]="models/cloud_5.prism properties/cloud.props -prop eventually_deploy -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:8842;Av:11;e:1e-8;d:0.01;p:0.001;post:2\"" \
	["mdsm-1"]="models/mdsm.prism properties/mdsm.props -prop player_1_deviate -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:62245;Av:2;e:1e-8;d:0.01;p:0.0684;post:5\"" \
	["mdsm-2"]="models/mdsm.prism properties/mdsm.props -prop player_2_deviate -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:62245;Av:2;e:1e-8;d:0.01;p:0.0684;post:5\"" \
	["team-form-3"]="models/team-form-3.prism properties/team-form.props -prop completed -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams \"S:12476;Av:3;e:1e-8;d:0.01;p:0.02040816326530612;post:49\""	
	)

declare -A black_models
black_models=(
	["consensus"]="models/consensus.2.prism properties/consensus.props -prop disagree -const K=2 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:272;Av:2;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["csma"]="models/csma.2-2.prism properties/csma.props -prop some_before -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:1038;Av:2;e:1e-8;d:0.01;p:0.25;post:4\"" \
	["firewire"]="models/firewire.true.prism properties/firewire.true.props -prop deadline -const delay=3,deadline=200 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:83153;Av:3;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["ij-3"]="models/ij.3.v1.prism properties/ij.3.v1.props -prop stable -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:7;Av:3;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["ij-10"]="models/ij.10.v1.prism properties/ij.10.v1.props -prop stable -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:1023;Av:10;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["pacman"]="models/pacman.nm properties/pacman.props -prop crsh -const MAXSTEPS=5 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:498;Av:5;e:1e-8;d:0.01;p:0.08;post:3\"" \
	["philosophers-3"]="models/philosophers-mdp.3.v1.prism properties/philosophers-mdp.3.v1.props -prop eat -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:956;Av:6;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["pnueli-zuck-3"]="models/pnueli-zuck.3.v1.prism properties/pnueli-zuck.v1.props -prop live -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:2701;Av:6;e:1e-8;d:0.01;p:0.5;post:2\"" \
	["rabin-3"]="models/rabin.3.prism properties/rabin.3.props -prop live -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:27766;Av:3;e:1e-8;d:0.01;p:0.03125;post:6\"" \
	["wlan-0"]="models/wlan.0.prism properties/wlan.props -prop sent -const COL=0 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:2954;Av:3;e:1e-8;d:0.01;p:0.0625;post:16\"" \
	["zeroconf"]="models/zeroconf.prism properties/zeroconf.props -prop correct_max -const N=20,K=2,reset=true -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:670;Av:3;e:1e-8;d:0.01;p:1.025262467191601E-4;post:6\"" \
	["cdmsn"]="models/cdmsn.prism properties/cdmsn.props -prop all_prefer_one -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:1240;Av:2;e:1e-8;d:0.01;p:0.059071881973375796;post:5\"" \
	["cloud-5"]="models/cloud_5.prism properties/cloud.props -prop eventually_deploy -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:8842;Av:11;e:1e-8;d:0.01;p:0.001;post:2\"" \
	["mdsm-1"]="models/mdsm.prism properties/mdsm.props -prop player_1_deviate -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:62245;Av:2;e:1e-8;d:0.01;p:0.0684;post:5\"" \
	["mdsm-2"]="models/mdsm.prism properties/mdsm.props -prop player_2_deviate -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:62245;Av:2;e:1e-8;d:0.01;p:0.0684;post:5\"" \
	["team-form-3"]="models/team-form-3.prism properties/team-form.props -prop completed -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 2 -colourParams \"S:12476;Av:3;e:1e-8;d:0.01;p:0.02040816326530612;post:49\""
)

SELECTION=$(whiptail --title "Select Benchmarks" --checklist \
	"Choose benchmarks to run using the Up, Down and Space keys" 25 35 16 \
"consensus" "MDP" OFF \
"csma" "MDP" OFF \
"firewire" "MDP" OFF \
"ij-3" "MDP" OFF \
"ij-10" "MDP" OFF \
"pacman" "MDP" OFF \
"philosophers-3" "MDP" OFF \
"pnueli-zuck-3" "MDP" OFF \
"rabin-3" "MDP" OFF \
"wlan-0" "MDP" OFF \
"zeroconf" "MDP" OFF \
"cdmsn" "Game" OFF \
"cloud-5" "Game" OFF \
"mdsm-1" "Game" OFF \
"mdsm-2" "Game" OFF \
"team-form-3" "Game" OFF 3>&1 1>&2 2>&3)

exitstatus=$?
if [ $exitstatus -ne 0 ]; then
	echo "Cancelled."
	exit 1
fi

REPS=$(whiptail --inputbox "Number of times each benchmark is to be run" 8 78 1 --title "Repetitions" 3>&1 1>&2 2>&3)

exitstatus=$?
if [ $exitstatus -ne 0 ]; then
    echo "Cancelled."
    exit 1
fi

echo "Selected benchmarks: " $SELECTION
echo "Repetitions: " $REPS

RUNCONFIG=()

for ITEM in $SELECTION
do
	RUNCONFIG+=("${grey_models[`eval echo $ITEM`]}")
	RUNCONFIG+=("${black_models[`eval echo $ITEM`]}")
done

if [ -d "$LOGDIR" ] 
then
    printf "$LOGDIR already exists.\n\n"
    read -p "In order to ensure consistency of the results, it is recommended to delete $LOGDIR. If you have not changed any setting (eg. timeout), you may wish to retain the logs from the previous runs. Delete $LOGDIR? (y/n) " -n 1 -r
	echo #
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
	    rm -rf "$LOGDIR"
	fi 
fi

mkdir -p $LOGDIR
printf "Saving logs to %s\n\n" $LOGDIR

# For each line in file
for line in "${RUNCONFIG[@]}"
do
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
		pueue add -- /usr/bin/time -v -o $LOGDIR/$filename.stat timeout $TIMEOUT $PRISM $line '2>&1' '|' tee $LOGDIR/$filename.log
		printf 'Output and stats in\n\t- %s.stat\n\t- %s.log\n\n' "$LOGDIR/$filename" "$LOGDIR/$filename"

	done;
done;

printf "Run 'pueue status' to see progress.\n\n\nLogs would become available in %s.\n\n" $LOGDIR