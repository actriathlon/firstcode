#!/bin/bash

l=3

for i in mono di tri tetra penta
do
	echo $i
	nm=$i
	echo $nm
	cp -r 100"$nm"client "$l"linker100"$nm"client

	cd "$l"linker100"$nm"client
	rm slurm*
	
	awk 'NR==n{$4=a}1' n=7 a=$l < 2systemsetupclient.txt > int2.txt

	rm 2systemsetupclient.txt
	mv int2.txt 2systemsetupclient.txt

	cd ../

done

