#!/bin/bash
#will take data in the ALL directories and select the .CSV file
#and feed it into the python file for processing
for file in ALL00??/*.CSV
do
#	echo $file
	python dataextract-t1.py $file
done
