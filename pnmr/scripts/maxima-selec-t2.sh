#!/bin/bash
#same as the *-t1.sh script except it feeds into a different python file
for file in ALL00??/*.CSV
do
	python dataextract-t2.py $file
done
