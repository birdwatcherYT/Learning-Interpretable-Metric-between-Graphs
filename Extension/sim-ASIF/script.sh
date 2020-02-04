#!/bin/bash

datasets=("AIDS" "BZR" "COX2")
maxpats=(8 8 8)

# datasets=("DHFR" "Mutagenicity" "NCI1")
# maxpats=(8 8 8)

datasets=("AIDS" "COX2" "BZR" "DHFR" "Mutagenicity" "NCI1")
maxpats=(10 10 10 10 10 10)


for i in `seq 0 $((${#datasets[@]}-1))`;do
	(for j in `seq 1 10`;do
		echo $j ${datasets[$i]} ${maxpats[$i]}
	done) | xargs -L 1 -P 5 ./run.sh
done

