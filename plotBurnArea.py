#!/bin/python -B

import matplotlib.pyplot as plt
import sys
from math import pi

def readBurnArea(filename):
    wb = []
    A = []
    with open(filename, "r") as fp:
        for line in fp:
            line = line.split()
            wb.append(float(line[0]))
            A.append(float(line[1]))
            
    return (wb, A)

def plotBurnArea(wb, A, filename):
    plt.figure()
    plt.plot(wb, A, label="Eikonal")
    plt.xlabel("Web burn - units of length")
    plt.ylabel("Burn area - units of area")
    plt.grid()
    plt.title(f"Burn area plot of {filename}")

def BATESfull(wb, R, r, L):
    if 2*wb >= (R - r) or 2*wb >= L:
        return 0
    
    aux1 = r + wb
    aux2 = R - wb
    A = 2*pi*aux1*(L - 2*wb)
    A += 2*pi*aux2*(L - 2*wb)
    A += 2*pi*(aux2*aux2 - aux1*aux1)
    return A

def BATES(wb, R, r, L):
    if wb >= (R - r) or 2*wb >= L:
        return 0
    
    aux = r + wb
    A = 2*pi*aux*(L - 2*wb)
    A += 2*pi*(R*R - aux*aux)
    return A

def getBATESDefault(wb, tag):
    A = []
    tag = tag.split("-")
    R = float(tag[0])
    r = float(tag[1])
    L = float(tag[2])
    for w in wb:
        A.append(BATES(w, R, r, L))
    return A

def getSphereDefault(wb, R):
    R = float(R)
    A = []
    for w in wb:
        A.append(0 if w > R else 4*pi*((R-w)**2))
    return A

def plotDefault(wb, A):
    plt.plot(wb, A, label="Analytical")
    plt.legend()
    

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("ERROR USING PLOT FILE")
		print("Usage of test file is:\n\n./plotBurnArea.py path_to_data.txt")
    
	wb, A = readBurnArea(sys.argv[1])
	plotBurnArea(wb, A, sys.argv[1])

	if len(sys.argv) == 4:
		if sys.argv[2] == "BATES":
			A_def = getBATESDefault(wb, sys.argv[3])
		elif sys.argv[2] == "Sphere":
			A_def = getSphereDefault(wb, sys.argv[3])
		plotDefault(wb, A_def)
 
	plt.show()