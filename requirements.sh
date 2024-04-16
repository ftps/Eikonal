#!/bin/bash

# Download latest version of gmsh
rm -rf gmsh-*
wget https://gmsh.info/bin/Linux/gmsh-stable-Linux64-sdk.tgz
tar -zxvf gmsh-stable-Linux64-sdk.tgz
rm gmsh-stable-Linux64-sdk.tgz

# Get correct GMSH version into make file
x="GMSH = ./$(ls | grep gmsh)" # | cut -d'-' -f 2
y=$(sed '6q;d' makefile)
sed -i -e "s%$y%$x%g" "makefile"