#!/bin/bash
echo "Rb87-------------------------"
for file in ~/Documents/PHY408/optical-pumping/data-iv/Rb87/ALL00??/*CH2.CSV
do
	echo $file
	python ~/Documents/PHY408/optical-pumping/scripts/data_extract.py $file
done
echo ""
echo "Rb85------------------------"
for file in ~/Documents/PHY408/optical-pumping/data-iv/Rb85/ALL00??/*CH2.CSV
do
	echo $file
        python ~/Documents/PHY408/optical-pumping/scripts/data_extract.py $file
done

