#!/bin/bash

gmsh_folder=$(ls | grep gmsh)
if [ ! -d "$gmsh_folder" ]; then
	./requirements.sh
fi

export LD_LIBRARY_PATH="./$(ls | grep gmsh)/lib;$LD_LIBRARY_PATH"

make
if [ $? -eq 0 ]; then
	echo -e "Error in compilation, exiting"
	exit 1
fi

rm -f test/sphere*.dat

./eikonal test/sphere.geo -p -s 100

if [ $? -eq 0 ]; then
	python ./plotBurnArea.py test/sphere_0.dat Sphere 5
fi