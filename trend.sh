#!/bin/bash

for fold in 0.01 0.05 0.15 0.25 0.35 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 
do
	cd $fold
	for i in 1 2 3 4 5
	do
		awk 'NR==n{$1=a}1' n=$i a=$i < trend.dat >> int2.txt
		rm trend.dat
		mv int2.txt trend.dat
	done
	#rm pq*
	cd ../
done


