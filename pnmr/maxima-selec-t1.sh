#!/bin/bash

for file in ALL00??/*.CSV
do
	python dataextract-t1.py $file
done
