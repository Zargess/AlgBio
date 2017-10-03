import numpy as np
import os, sys
import os.path
import time
import global_linear as linear
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util

score = {}
gap_cost = 5
gap_symbol = "-"

T = None
A = ""
B = ""
C = ""

def sub(a, b):
    return score[(a.upper(), b.upper())]

	
def find_center(sequences):
	bestSum = float('inf')
	bestSeq = -1 
	for i in range(0, len(sequences)):
		sum = 0
		for j in range (0, len(sequences)):
			if(i != j):
				s = linear.runAlgo(sequences[i], sequences[j], s_mat=score, gc=gap_cost)
				
				sum +=  s
		print("I: {}, SUM: {}".format(i,sum))		
		if(sum < bestSum):
			bestSum = sum
			bestSeq = i
	return bestSeq
	
def construct_alignment(sequences, bestSeqIdx):
	center = sequences[bestSeqIdx]
	oldM = np.empty((1, len(center)), dtype='int')
	
	# Convert to ASCII codes
	for i in range(0, len(sequences[bestSeqIdx])):
		oldM[0, i] = ord(center[i])
	
	
	for i in range(0, len(sequences)):
		if(i != bestSeqIdx): 
			A = linear.runAlgoWithBacktrack(sequences[bestSeqIdx], sequences[i], s_mat=score, gc=gap_cost)
			newM = extend_alignment(oldM, A)
			oldM = newM
	return oldM

def extend_alignment(oldM, A):
	m = max([len(w) for w in A])
	newM = np.resize(oldM, (oldM.shape[0]+1, oldM.shape[1] + len(A[1])))
	diff = 0
	j = 0
	for c in range (0, m):
		while(j <= n):
			if A[0][c] != '-':
				if chr(oldM[0, j]) == '-':
					insertOldColAndGap(newM, oldM, j, j+diff)
					j += 1
				else:
					insertOldColAndSym(newM, oldM, j, j+diff, A[1][c])
					j += 1
					break			
			else:
				if chr(oldM[0,j]) != '-':
					insertGapColAndSym(newM, j+diff, A[1][c])
					diff += 1
					break
				else:
					insertOldColAndSym(newM, oldM, j, j+diff, A[1][c])
					j += 1
					break
		if(j > n):
			insertGapColAndSym(newM, c, A[1][c])
	
	while(j <= n):
		insertOldColAndGap(newM, j+diff)
		j += 1
	return newM

def insertOldColAndGap(newM, oldM, j, jdiff):
	newM[0:oldM.shape[0],jdiff] = oldM[:,j]
	newM[newM.shape[0]-1, jdiff] = ord('-')

def insertOldColAndSym(newM, oldM, j, jdiff, sym):
	newM[0:oldM.shape[0],jdiff] = oldM[:,j]
	newM[newM.shape[0]-1, jdiff] = ord(sym)
	
def insertGapColAndSym(newM, jdiff, sym):
	newM[0:newM.shape[0]-1,jdiff] = ord('-')
	newM[newM.shape[0]-1, jdiff] = ord(sym)
	
	
def decode_matrix(mat):
	f = lambda x : chr(x)
	func = np.vectorize(f)
	return func(mat)
if __name__ == "__main__":
    args = sys.argv
    
    
    score, _ = util.read_score_matrix_and_alphabet(args[1])
    gap_cost = args[2]
    fastaDictionary = util.read_fasta_file(args[3])
    sequences = [s.replace(" ", "") for s in fastaDictionary.values()]

    #A = next(sequences)
    #B = next(sequences)
    #C = next(sequences)
    # Length of the strings plus 1
    n = len(A)
    m = len(B)
    l = len(C)

    center = find_center(sequences)
    m = construct_alignment(sequences, center)
    print(decode_matrix(m))

    