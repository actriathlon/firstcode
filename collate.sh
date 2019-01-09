#!/bin/bash


for fold in 0.01 0.05 0.15 0.25 0.35 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
do
	cd $fold
	cat pqclientresults1.dat > trend.dat
	cat pqclientresults2.dat >> trend.dat
	cat pqclientresults3.dat >> trend.dat
	cat pqclientresults4.dat >> trend.dat
	cat pqclientresults5.dat >> trend.dat
	cd ../
done


