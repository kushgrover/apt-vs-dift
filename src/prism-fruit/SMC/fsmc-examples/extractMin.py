#!/usr/bin/python
import subprocess
import re
import os, sys
import time
import pipes


def main():

    fl = open(sys.argv[1],'r')
    max_sum_rate = 0.0
    max_state = 0
    sum_rate = 0.0
    prev_state = 0
    min_rate = 10000.0;
    lines = fl.readlines()
    
    for line in lines[1:]:
        elems = line.split(' ')
        state = int(elems[0])
        rate =  float(elems[2].replace(',','.'))

        if (state != prev_state):
            prev_state = state
            if (sum_rate > max_sum_rate):
                max_sum_rate = sum_rate
                max_state = prev_state
            sum_rate = 0.0


        sum_rate += rate
        if (rate < min_rate):
            min_rate = rate

    print("max_sum_rate="+str(max_sum_rate)+", max_state="+str(max_state)+", min_rate="+str(min_rate))
            

if __name__ == "__main__":
    main()
