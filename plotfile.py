import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab as pl
from matplotlib import collections  as mc
import sys
from os import listdir
from os.path import isfile, join
import os
import subprocess
import itertools

allMarkers = list(('x', '+', '*', 'o', ','))

def plotRuntimeYDividedByLogN(allData):
	#allData is a list of tuples with form ((n, value), legend)
	marker = itertools.cycle(allMarkers)
	
	plt.figure()
	
	for ((inputsize, value), legend) in allData:
		plt.plot(inputsize, value , marker=next(marker), linestyle ='--', label=legend)
		highestY = value[len(value)-1]
	
	plt.legend(loc=2,borderaxespad=0.)
	plt.xlabel('Inputsize')
	plt.ylabel('Runtime (ms)')
	plt.xscale("log", basex=2)
	plt.axvline(x=2**20, color='black', linestyle=':', label="L3 cache size")
	ymin, ymax = plt.gca().get_ylim();
	plt.text(2**20, highestY*1.4, "L3 cache size",  rotation=270)
	
def plotValues(allData, ylabel):
	#allData is a list of tuples with form ((n, value), legend)
	plt.figure()
	
	marker = itertools.cycle(allMarkers)
	
	print(allData)
	for ((inputsize, value), legend) in allData:
		plt.plot(inputsize, value , marker=next(marker), linestyle ='--', label=legend)
		highestY = value[len(value)-1]
	
	plt.legend(loc=2, borderaxespad=0.)
	plt.xlabel('Inputsize')
	plt.ylabel(ylabel)
	plt.xscale("log", basex=2)
	plt.show()
	#plt.axvline(x=2**13, color='black', linestyle='--', label="L1 cache size")
	#plt.text(2**13, highestY*1.4	, "L1 cache size",  rotation=270)

	

def readFile(filename):	
	inputSizes = []
	values = []
	with open(filename) as f:
		for line in f:
			thisLine = line.split(',')
			inputSizes.append(float(thisLine[0]))
			values.append(float(thisLine[1]))
			
	return inputSizes, values
	
	
def plot(folder):
	#Compiling and running c++ implementation
	#subprocess.call(["g++",  "HelloTest.cpp", "-o", "hello.exe"])
	#subprocess.call(["hello.exe"])
	
	# List elements are tuples with layout ((N, value), legend)
	#Plotting runtime divided by log(n)
	allRuntimes = []
	allRuntimes.append((readFile("./" + folder + "/data/Simpletimes.txt"), "Sorted array"))
	allRuntimes.append((readFile("./" + folder + "/data/BFStimes.txt"), "BFS Structure"))
	allRuntimes.append((readFile("./" + folder + "/data/BFSLItimes.txt"), "BFSLI Structure"))
	plotRuntimeYDividedByLogN(allRuntimes)
	
	
	#Plotting branch mispredictions
	allBranchMissPreds = []
	allBranchMissPreds.append((readFile("./" + folder + "/data/SimpleBM.txt"), "Sorted array"))
	allBranchMissPreds.append((readFile("./" + folder + "/data/BFSBM.txt"), "BFS Structure"))
	allBranchMissPreds.append((readFile("./" + folder + "/data/BFSLIBM.txt"), "BFSLI Structure"))
	plotValues(allBranchMissPreds, 'Branch Mispredictions')
	
	#Plotting L1 cache misses
	allCachemisses = []
	allCachemisses.append((readFile("./" + folder + "/data/SimpleL1.txt"), "Sorted array"))
	allCachemisses.append((readFile("./" + folder + "/data/BFSL1.txt"), "BFS Structure"))
	allCachemisses.append((readFile("./" + folder + "/data/BFSLIL1.txt"), "BFSLI Structure"))
	plotValues(allCachemisses, 'L1 cache misses')
	
	#Plotting L2 cache misses
	allCachemisses = []
	allCachemisses.append((readFile("./" + folder + "/data/SimpleL2.txt"), "Sorted array"))
	allCachemisses.append((readFile("./" + folder + "/data/BFSL2.txt"), "BFS Structure"))
	allCachemisses.append((readFile("./" + folder + "/data/BFSLIL2.txt"), "BFSLI Structure"))
	plotValues(allCachemisses, 'Instructions')
	
	plt.show()
	
if __name__ == "__main__":
	plot(sys.argv[1])
