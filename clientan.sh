#!/bin/bash

starttime=200
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13
do
	echo $j 

	echo $starttime
	cd new"$j"_1
	#cp movetagged.vtf m.vtf
	#cp energy.dat e.dat
	#cp coms.dat c.dat
	for i in 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190
	do
		declare -i a	
		a=$i+$starttime
		declare -i b
		b=$a-1
		echo $a
		../../../comver 1 $b 1 1 0
		
		#for file in inputp.dat;
		#do
		#	while IFS='|' read -r w x y z
		#done
		cat 'inputp.dat' | while read w x y z
		do
			echo $w $x $y $z			
		
		
			#echo $w $x $y $z
			#cd "$a"monoclient
			#echo "$a"monoclient
			#sbatch cl.slurm

			cd "$a"monoclient
                        #sbatch cl.slurm
                        python ../../../../clientanalysis.py $w $x 1 $y $z
			#cd "$a"diclient 	
			#sbatch cl.slurm
			#python ../../../../clientanalysis.py $w $x 2 $y $z
			#cd ../"$a"triclient
                	#sbatch cl.slurm
			#python ../../../../clientanalysis.py $w $x 3 $y $z
			cd ../
		done
	done
	cd ../
done

