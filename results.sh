mkdir -p results
rm -rf results/*

mkdir -p results/random_10
mkdir -p results/random_10/output
mkdir -p results/random_10/plots

echo "Random acyclic graph of size 10:"
echo "      Generating random model of size 10"
python3 src/random_model_gen.py results/random_10/model.prism 10 -acyclic
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/random_10/model.prism props/init.props -outputDir results/random_10 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/random_10/output/step1.txt
echo "      done"

colourParams=$(<results/random_10/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/random_10/model.prism props/init.props -outputDir results/random_10 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/random_10/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/random_10/model.prism props/main.props -pareto > results/random_10/model_pareto.txt
src/prism-games/prism/bin/prism results/random_10/best.prism props/main.props -pareto > results/random_10/best_pareto.txt
src/prism-games/prism/bin/prism results/random_10/worst.prism props/main.props -pareto > results/random_10/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/random_10/ results/random_10/plots/ > results/random_10/output/plot.txt
echo "      done"



mkdir -p results/random_10_cyc
mkdir -p results/random_10_cyc/output
mkdir -p results/random_10_cyc/plots

echo "Random graph of size 10:"
echo "      Generating random model of size 10"
python3 src/random_model_gen.py results/random_10_cyc/model.prism 10
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/random_10_cyc/model.prism props/init.props -outputDir results/random_10_cyc -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/random_10_cyc/output/step1.txt
echo "      done"

colourParams=$(<results/random_10_cyc/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/random_10_cyc/model.prism props/init.props -outputDir results/random_10_cyc -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/random_10_cyc/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/random_10_cyc/model.prism props/main.props -pareto > results/random_10_cyc/model_pareto.txt
src/prism-games/prism/bin/prism results/random_10_cyc/best.prism props/main.props -pareto > results/random_10_cyc/best_pareto.txt
src/prism-games/prism/bin/prism results/random_10_cyc/worst.prism props/main.props -pareto > results/random_10_cyc/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/random_10_cyc/ results/random_10_cyc/plots/ > results/random_10_cyc/output/plot.txt
echo "      done"




mkdir -p results/random_100
mkdir -p results/random_100/output
mkdir -p results/random_100/plots

echo "Random acyclic graph of size 100:"
echo "      Generating random model of size 100"
python3 src/random_model_gen.py results/random_100/model.prism 100 -acyclic
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/random_100/model.prism props/init.props -outputDir results/random_100 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/random_100/output/step1.txt
echo "      done"

colourParams=$(<results/random_100/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/random_100/model.prism props/init.props -outputDir results/random_100 -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/random_100/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/random_100/model.prism props/main.props -pareto > results/random_100/model_pareto.txt
src/prism-games/prism/bin/prism results/random_100/best.prism props/main.props -pareto > results/random_100/best_pareto.txt
src/prism-games/prism/bin/prism results/random_100/worst.prism props/main.props -pareto > results/random_100/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/random_100/ results/random_100/plots/ > results/random_100/output/plot.txt
echo "      done"



mkdir -p results/random_100_cyc
mkdir -p results/random_100_cyc/output
mkdir -p results/random_100_cyc/plots

echo "Random graph of size 100:"
echo "      Generating random model of size 100"
python3 src/random_model_gen.py results/random_100_cyc/model.prism 100
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/random_100_cyc/model.prism props/init.props -outputDir results/random_100_cyc -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/random_100_cyc/output/step1.txt
echo "      done"

colourParams=$(<results/random_100_cyc/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/random_100_cyc/model.prism props/init.props -outputDir results/random_100_cyc -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/random_100_cyc/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/random_100_cyc/model.prism props/main.props -pareto > results/random_100_cyc/model_pareto.txt
src/prism-games/prism/bin/prism results/random_100_cyc/best.prism props/main.props -pareto > results/random_100_cyc/best_pareto.txt
src/prism-games/prism/bin/prism results/random_100_cyc/worst.prism props/main.props -pareto > results/random_100_cyc/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/random_100_cyc/ results/random_100_cyc/plots/ > results/random_100_cyc/output/plot.txt
echo "      done"





mkdir -p results/screen_grab
mkdir -p results/screen_grab/output
mkdir -p results/screen_grab/plots

echo "Screen Grab:"
echo "      Generating PRISM model for Screen Grab graph"
python3 src/matrix_model_gen.py examples/screen_grab/screen_grab.csv results/screen_grab/model.prism
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/screen_grab/model.prism props/init.props -outputDir results/screen_grab -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/screen_grab/output/step1.txt
echo "      done"

colourParams=$(<results/screen_grab/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/screen_grab/model.prism props/init.props -outputDir results/screen_grab -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/screen_grab/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/screen_grab/model.prism props/main.props -pareto > results/screen_grab/model_pareto.txt
src/prism-games/prism/bin/prism results/screen_grab/best.prism props/main.props -pareto > results/screen_grab/best_pareto.txt
src/prism-games/prism/bin/prism results/screen_grab/worst.prism props/main.props -pareto > results/screen_grab/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/screen_grab/ results/screen_grab/plots/ > results/screen_grab/output/plot.txt
echo "      done"


mkdir -p results/nation_attack
mkdir -p results/nation_attack/output
mkdir -p results/nation_attack/plots

echo "Nation Attack:"
echo "      Generating PRISM model for Nation Attack graph"
python3 src/matrix_model_gen.py examples/nation_attack/nation_attack.csv results/nation_attack/model.prism
echo "      done"

echo "      Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism results/nation_attack/model.prism props/init.props -outputDir results/nation_attack -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > results/nation_attack/output/step1.txt
echo "      done"

colourParams=$(<results/nation_attack/colourParams.txt)

echo "      Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism results/nation_attack/model.prism props/init.props -outputDir results/nation_attack -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > results/nation_attack/output/step2.txt
echo "      done"

echo "      Generating Pareto frontiers"
src/prism-games/prism/bin/prism results/nation_attack/model.prism props/main.props -pareto > results/nation_attack/model_pareto.txt
src/prism-games/prism/bin/prism results/nation_attack/best.prism props/main.props -pareto > results/nation_attack/best_pareto.txt
src/prism-games/prism/bin/prism results/nation_attack/worst.prism props/main.props -pareto > results/nation_attack/worst_pareto.txt
echo "      done"

echo "      Generating plots"
python3 src/plot_pareto.py results/nation_attack/ results/nation_attack/plots/ > results/nation_attack/output/plot.txt
echo "      done"

