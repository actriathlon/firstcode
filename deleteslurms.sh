#!/bin/bash
j=5
#starttime=400

for fold in 0.01 0.05 0.15 0.25 0.35 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 
do

	cd $fold
	for starttime in 100 300 #500 700 #190 290 #2 3 4 5 6 7 8 9 10 11 12 13
	do
		echo $j 

		echo $starttime
		#cd new"$j"_1
		
		for i in 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 #200
		do
			declare -i a	
			a=$i+$starttime
			declare -i b
			b=$a-1
			echo $a
		
		
			cd 3linker"$a"monoclient
                	rm slurm*
			cd ../3linker"$a"diclient
			rm slurm*

                	cd ../3linker"$a"triclient
			rm slurm*
                	cd ../


			cd 3linker"$a"tetraclient
			rm slurm*	
			cd ../3linker"$a"pentaclient 	
			rm slurm*
			cd ../
		done
		#cd ../
	done




	python ../clientcombine.py 1 $j pqbindingaverages1.dat pqclientresults1.dat
	python ../clientcombine.py 2 $j pqbindingaverages2.dat pqclientresults2.dat
	python ../clientcombine.py 3 $j pqbindingaverages3.dat pqclientresults3.dat
	python ../clientcombine.py 4 $j pqbindingaverages4.dat pqclientresults4.dat
	python ../clientcombine.py 5 $j pqbindingaverages5.dat pqclientresults5.dat

	cd ../

done

