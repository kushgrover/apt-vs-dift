mkdir -p temp
rm -rf temp/*

mkdir -p temp/output
mkdir -p temp/plots

echo "Generating random model of size $1"
python3 src/random_model_gen.py temp/model.prism $1
echo "done"

# echo "Generating PRISM model from adjacency matrix"
# python3 src/matrix_model_gen.py examples/screen_grab/screen_grab.csv temp/model.prism
# # python3 src/matrix_model_gen.py examples/nation_attack/nation_attack.csv temp/model.prism
# echo "done"

echo "Generating colour parameters"
src/prism-fruit/Games-DQL/bin/prism temp/model.prism props/init.props -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 3 > temp/output/step1.txt
echo "done"

colourParams=$(<temp/colourParams.txt)

echo "Generating best and worst case games"
src/prism-fruit/Games-DQL/bin/prism temp/model.prism props/init.props -outputDir temp -heuristic RTDP_ADJ -RTDP_ADJ_OPTS 1 -colourParams $colourParams > temp/output/step2.txt
echo "done"

echo "Generating Pareto frontiers"
src/prism-games/prism/bin/prism temp/model.prism props/main.props -pareto > temp/model_pareto.txt
src/prism-games/prism/bin/prism temp/best.prism props/main.props -pareto > temp/best_pareto.txt
src/prism-games/prism/bin/prism temp/worst.prism props/main.props -pareto > temp/worst_pareto.txt
echo "done"

echo "Generating plots"
python3 src/plot_pareto.py temp/ temp/plots/ > temp/output/plot.txt
echo "done"


