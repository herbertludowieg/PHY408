#!/bin/bash

for file in ALL00??/*.CSV
do
	python dataextract-t2.py $file
done
