#!/bin/bash

# Finds minimum transition probabilities in the model

prism $@ -exporttrans tr
./extractMin.py tr
# cat tr | awk '{print $3'} | sort -n > probs
# head -2 probs

