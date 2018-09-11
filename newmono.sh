#!/bin/bash

starttime=200
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13
do
	echo $j 

	echo $starttime
	cd new"$j"_1
	for i in 0 #10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190
	do
		declare -i a	
		a=$i+$starttime
		echo $a
		#cp -r 200monoclient "$a"monoclient
		#cd "$a"monoclient 	
		cd "$a"diclient
		#sed '/freeze/s/200/'$a'/g' setupclient.txt > setupclient2.txt 
                sed '/film/s/TRUE/FALSE/g' setupclient.txt > setupclient2.txt
                rm setupclient.txt
                mv setupclient2.txt setupclient.txt
		cd ../"$a"triclient
		sed '/film/s/TRUE/FALSE/g' setupclient.txt > setupclient2.txt
                rm setupclient.txt
                mv setupclient2.txt setupclient.txt


		#sbatch cl.slurm
		cd ../

	done
	cd ../
done

