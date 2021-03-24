# Instructions for benchmarking

1. Create a configuration file, say `configurations.txt`, each line containing `(MODELFILE; PROPERTYFILE; PROPERTYNAME; CONSTS; COLOR; EXPECTEDRESULT)`; e.g. `(crowds.prism; crowds.props; positive; TotalRuns=5,CrowdSize=10; 1; 0.10478678887151971)`
2. Run `generate.sh config/all.txt > all-generated.txt` to obtain a file containing the run configurations with the colour parameters computed automatically
3. Use `run_experiments.sh all-generated.txt N` to run each experiment in `all-generated.txt` N times.

### Notes
1. This file is assumed to be located in a `data` directory within the PRISM directory containing `bin`, `src` and so on.
2. The `data` directory must contain folder `models` and `properties` with a flat structure containing all the models and the properties respectively.
2. The log directory defaults to `data/logs` and the path to PRISM defualts to `../bin/prism`
