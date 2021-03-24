#!/usr/bin/env bash

if [[ -z "$1" ]]; then
  echo "Need to pass a configuration file";
else
  grep -v '^#' $1 | while IFS= read -r line; do
    echo ${line} | tr -d '( )' | awk -F';' '{printf "models/%s properties/%s -prop %s",
     $1, $2, $3; 
     if (length($4) != 0) CONST=" -const " $4;
     printf CONST;
     printf " -heuristic RTDP_ADJ";
     printf " -RTDP_ADJ_OPTS %s", $5;
     while ( ( "../bin/prism models/"$1" properties/"$2" "CONST" -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 | tail -1" | getline result ) > 0 ) { 
     	if (length(result) != 0) {printf " -colourParams \""result"\"" }
     }; 
     print ""}'
  done
fi
