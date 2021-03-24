sed -E 's/OPTS 0/OPTS 1/g' white-generated.txt > grey-generated.txt
sed -E 's/OPTS 0/OPTS 2/g' white-generated.txt > black-generated.txt