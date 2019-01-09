#!/bin/bash

for fold in 0.01 0.05 0.15 0.1 0.2 0.25 0.3 0.35 0.4 0.5 0.6 0.7 0.8

do
        cd $fold

	rm pq*
	rm trend.dat
	rm trendl.dat
	cd ../

done
