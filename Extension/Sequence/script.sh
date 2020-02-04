#!/bin/bash

datasets=("promoters" "splice")
maxpats=(30 30)

for i in `seq 0 $((${#datasets[@]}-1))`;do
	(for j in `seq 1 10`;do
		echo $j ${datasets[$i]} ${maxpats[$i]}
	done) | xargs -L 1 -P 5 ./run.sh
done
