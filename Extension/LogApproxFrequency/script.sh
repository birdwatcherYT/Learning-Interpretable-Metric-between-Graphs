#!/bin/bash

datasets=("AIDS" "BZR" "DD" "DHFR" "FRANKENSTEIN" "Mutagenicity" "NCI1" "COX2" "DBLP_v1")
maxpats=(30 15 30 15 15 10 15 15 30)

# datasets=("AIDS" "BZR" "DD" "FRANKENSTEIN")
# maxpats=(30 15 30 15)

datasets=("DBLP_v1")
maxpats=(30)

datasets=("AIDS")
maxpats=(30)

for i in `seq 0 $((${#datasets[@]}-1))`;do
	(for j in `seq 1 10`;do
		echo $j ${datasets[$i]} ${maxpats[$i]}
	done) | xargs -L 1 -P 5 ./run.sh
done

# for i in `seq 0 $((${#datasets[@]}-1))`;do
# 	(for j in `seq 1 10`;do
# 		echo $j ${datasets[$i]} ${maxpats[$i]}
# 	done) | xargs -L 1 -P 5 ./run.sh
# done
