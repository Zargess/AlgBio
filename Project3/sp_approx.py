import numpy as np
import os, sys
import os.path
import itertools
import time
import global_linear as linear
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util
import msa_sp_score_3k

score = {}
gap_cost = 5
gap_symbol = "-"

T = None
A = ""
B = ""
C = ""

def sub(a, b):
    return score[(a.upper(), b.upper())]

def pp_matrix(M):
	for i in range(0, M.shape[0]):
			print (''.join(decode_matrix(M)[i,:]))

	
def find_center(sequences):
	bestSum = float('inf')
	bestSeq = -1 
	for i in range(0, len(sequences)):
		sum = 0
		for j in range (0, len(sequences)):
			if(i != j):
				#print(sequences[i])
				#print(sequences[j])
				s = linear.runAlgo(sequences[i], sequences[j], s_mat=score, gc=gap_cost)
				#print(s)
				sum +=  s
		#print(sum)
		if(sum < bestSum):
			bestSum = sum
			bestSeq = i
	
	return bestSeq

def construct_alignment_fac(sequences, bestSeqIdx):
	center = sequences[bestSeqIdx]
#	oldM = np.empty((1, len(center)), dtype='int')
	
	# Convert to ASCII codes
#	for i in range(0, len(sequences[bestSeqIdx])):
#		oldM[0, i] = ord(center[i])
	
	indicies = [0,1,2,3,4]
	indicies.remove(bestSeqIdx)
	allPerms = list(itertools.permutations(indicies))
	c = 0
	for perm in allPerms: 		
		oldM = np.empty((1, len(center)), dtype='int')
		
		# Convert to ASCII codes
		for i in range(0, len(sequences[bestSeqIdx])):
			oldM[0, i] = ord(center[i])
		c += 1
		for i in perm:
			A = linear.runAlgoWithBacktrack(sequences[bestSeqIdx], sequences[i], s_mat=score, gc=gap_cost)
			#print(A)
			#print(decode_matrix(oldM))
			newM = extend_alignment(oldM, A)
			#print(decode_matrix(newM))
			oldM = newM
		dm = decode_matrix(oldM)
		content = []
		for i in range(0, dm.shape[0]):
			content.append(str(''.join(dm[i,:])))
		util.write_fasta_file("testfile", content)
		print("Score for perm: {} is: {}".format(perm, msa_sp_score_3k.compute_sp_score("testfile.fa")))
		
		
	return oldM

	
def construct_alignment(sequences, bestSeqIdx):
	center = sequences[bestSeqIdx]
	oldM = np.empty((1, len(center)), dtype='int')
	
	# Convert to ASCII codes
	for i in range(0, len(sequences[bestSeqIdx])):
		oldM[0, i] = ord(center[i])
	
	
	for i in range(0, len(sequences)):
		if(i != bestSeqIdx): 
			A = linear.runAlgoWithBacktrack(sequences[bestSeqIdx], sequences[i], s_mat=score, gc=gap_cost)
			#print(A)
			#print(decode_matrix(oldM))
			newM = extend_alignment(oldM, A)
			#print(decode_matrix(newM))
			oldM = newM
	return oldM

def extend_alignment(oldM, A):
	m = max([len(w) for w in A])
	n = oldM.shape[1]
	#print("Alignment: ")
	#print(A)
	newM = np.zeros((oldM.shape[0]+1, oldM.shape[1] + len(A[1])), dtype='int')
	diff = 0
	j = 0
	for c in range (0, m):
		while(j < n):
			if A[0][c] != '-':
				if chr(oldM[0, j]) == '-':
					insertOldColAndGap(newM, oldM, j, j+diff)
					j += 1
					#print("First case")
				else:
					insertOldColAndSym(newM, oldM, j, j+diff, A[1][c])
					j += 1
					#print("second case")
					break
			else:
				if chr(oldM[0,j]) != '-':
					insertGapColAndSym(newM, j+diff, A[1][c])
					diff += 1
					#print("third case")
					break
				else:
					insertOldColAndSym(newM, oldM, j, j+diff, A[1][c])
					j += 1
					#print("fourth case")
					break
		if(j > n):
			insertGapColAndSym(newM, c, A[1][c])
			#print("is it here?")
	#if(j+diff >= c): 
	#	newM = np.resize(newM, (newM.shape[0], newM.shape[1]-len(A[1])))
	#if(j+diff <= newM.shape[1]):
	#	newM = np.resize(newM, (newM.shape[0], j+diff))
	#print(j)
	#print(diff)
	#print(newM.shape)
	#print("After extension")
	#pp_matrix(newM)	
	#while(j < n):
	#	insertOldColAndGap(newM, oldM, j, j+diff)
	#	j += 1
		
	
	return newM

def insertOldColAndGap(newM, oldM, j, jdiff):
	newM[0:oldM.shape[0],jdiff] = oldM[:,j]
	#print(decode_matrix(oldM[:,j]))
	#print(decode_matrix(newM[:,jdiff]))
	newM[newM.shape[0]-1, jdiff] = ord('-')
	#print(decode_matrix(newM[:,jdiff]))

def insertOldColAndSym(newM, oldM, j, jdiff, sym):
	#print("NewM before inserting old col and sym")
	#pp_matrix(newM)
	newM[0:oldM.shape[0],jdiff] = oldM[:,j]
	
	newM[newM.shape[0]-1, jdiff] = ord(sym)
	#pp_matrix(newM)
	#print("NewM after inserting old col and sym")
	
	
def insertGapColAndSym(newM, jdiff, sym):
	newM[0:newM.shape[0]-1,jdiff] = ord('-')
	newM[newM.shape[0]-1, jdiff] = ord(sym)
	
	
def decode_matrix(mat):
	f = lambda x : chr(x)
	func = np.vectorize(f)
	return func(mat)
	
def compute_score(s_mat, gc, sequences):
	global score
	score = s_mat
	global gap_cost
	gap_cost = gc
	
	center = find_center(sequences)
	m = construct_alignment(sequences, center)
	dm = decode_matrix(m)
	content = []
	for i in range(0, dm.shape[0]):
		content.append(str(''.join(dm[i,:])))
		
	util.write_fasta_file("testfile", content)
	return msa_sp_score_3k.compute_sp_score("testfile.fa")
	
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
    #sequences = ["AGTACC", "AGATCC", "TTATG", "ACTACG", "ACTTGG"] #, "AAACTA", "AGCTAA"]
    center = find_center(sequences)
    
    m = construct_alignment(sequences, center)
    dm = decode_matrix(m)
	
    content = []
    for i in range(0, dm.shape[0]):
        content.append(str(''.join(dm[i,:])))
    util.write_fasta_file("testfile", content)
    print("This is center: {}".format(sequences[center]))
    print(msa_sp_score_3k.compute_sp_score("testfile.fa"))
	
    