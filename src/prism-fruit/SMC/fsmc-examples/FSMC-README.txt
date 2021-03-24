0. ABOUT

This is PRISM-adaptive, an extension of the PRISM model checker to statistical
model checking of full linear temporal logic and mean payoff.

Author:
Przemyslaw Daca (przemek.daca@gmail.com)
----------------------------
1. LICENSE

The tool is distributed under the GNU General Public License (GPL).

----------------------------
2. INSTALL

Follow instructions from README_PRISM.tex. In short
 - enter the PRISM directory and type "make"
 - to run, execute bin/xprism or bin/prism

----------------------------
3. RUN BENCHMARKS

To run the benchmark Python 2.7 and the Linux "time" utility are required.
Execute the following commands:
./run_experiments.py -e experimentsReach.set -b benchmarksReach.set
./run_experiments.py -e experimentsLTL.set -b benchmarksLTL.set
./run_experiments.py -e experimentsRew.set -b benchmarksRew.set

----------------------------
4. RUN INSTANCE

Example:
bin/prism examples/crowds/crowds_nodl.pm examples/crowds/positive_qual.pctl -simminprob 0.0666 -sim -simcandidate -simmethod sprte -simfalsecnd 0.001 -simratio 0.5 -simconf 0.01 -simwidth 0.01 -simpathlen 10000000 -simcheckbound 1000  -const TotalRuns=6,CrowdSize=15


