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
	ymin, ymax = plt.gca().get_ylim()
	plt.text(2**20, highestY*1.4, "L3 cache size",  rotation=270)

def create_normalized_list(inputsize, values):
	res = []
	for (n, t) in zip(inputsize, values):
		res.append(t / (n * n))

	return res

def plotValues(allData, ylabel):
	#allData is a list of tuples with form ((n, value), legend)
	plt.figure()
	
	marker = itertools.cycle(allMarkers)
	
	print(allData)
	for ((inputsize, value), legend) in allData:
		plt.plot(inputsize, value, marker=next(marker), linestyle ='--', label=legend)
		highestY = value[len(value)-1]
	
	plt.axhline(y=(4/3), color='r', linestyle='-')
	plt.axhline(y=1, color='b', linestyle='-')
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
	allRuntimes.append((readFile("./" + folder + "/results.txt"), "Global Affine"))
	plotValues(allRuntimes, 'Approx/Exact ratio')
	
	plt.show()
	
if __name__ == "__main__":
	plot(sys.argv[1])
