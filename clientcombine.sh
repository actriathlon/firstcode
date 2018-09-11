#!/bin/bash

starttime=200
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13
do
	cd new"$j"_1
	echo new"$j"_1
	echo $j 
	python ../../../clientcombine.py $j
	cd ../
done

