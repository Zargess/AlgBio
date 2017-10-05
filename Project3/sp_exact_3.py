import numpy as np
import os, sys
import os.path
import time
import global_linear as linear
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util

score = {}
gap_cost = 0
gap_symbol = "-"

T = None
A = ""
B = ""
C = ""

def sub(a, b):
    return score[(a.upper(), b.upper())]

def SP(a, b, c):
    if a != gap_symbol and b != gap_symbol and c != gap_symbol:
        return sub(a, b) + sub(b, c) + sub(a, c)
    if a != gap_symbol and b != gap_symbol and c == gap_symbol:
        return sub(a, b) + gap_cost + gap_cost
    if a != gap_symbol and b == gap_symbol and c != gap_symbol:
        return gap_cost + gap_cost + sub(a, c)
    if a == gap_symbol and b != gap_symbol and c != gap_symbol:
        return gap_cost + sub(b, c) + gap_cost
    return 2 * gap_cost

	
def calc_T(n, m, l):
    for i in range(0, n+1):
        for j in range(0, m+1):
            for k in range(0, l+1):
                v0 = float("inf")
                v1 = float("inf")
                v2 = float("inf")
                v3 = float("inf")
                v4 = float("inf")
                v5 = float("inf")
                v6 = float("inf")
                v7 = float("inf")

                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    v1 = T[i-1, j-1, k-1] + SP(A[i-1], B[j-1], C[k-1])
                if i > 0 and j > 0 and k >= 0:
                    v2 = T[i-1, j-1, k] + SP(A[i-1], B[j-1], gap_symbol)
                if i > 0 and j >= 0 and k > 0:
                    v3 = T[i-1, j, k-1] + SP(A[i-1], gap_symbol, C[k-1])
                if i >= 0 and j > 0 and k > 0:
                    v4 = T[i, j-1, k-1] + SP(gap_symbol, B[j-1], C[k-1])
                if i > 0 and j >= 0 and k >= 0:
                    v5 = T[i-1, j, k] + SP(A[i-1], gap_symbol, gap_symbol)
                if i >= 0 and j > 0 and k >= 0:
                    v6 = T[i, j-1, k] + SP(gap_symbol, B[j-1], gap_symbol)
                if i >= 0 and j >= 0 and k > 0:
                    v7 = T[i, j, k-1] + SP(gap_symbol, gap_symbol, C[k-1])

                T[i, j, k] = min(v0, v1, v2, v3, v4, v5, v6, v7)

	
def backtrack(i, j, k, output1, output2, output3):
	if (i > 0) and (j > 0) and (k > 0) and (T[i,j, k] == (T[i-1, j-1, k-1] + SP(A[i-1], B[j-1], C[k-1]))):
		return backtrack(i-1, j-1, k-1, A[i-1] + output1, B[j-1] + output2, C[k-1] + output3)		
	if (i > 0) and (j > 0) and (k >= 0) and (T[i,j,k] == (T[i-1, j-1, k] + SP(A[i-1], B[j-1], gap_symbol))):
		return backtrack(i-1, j-1, k, A[i-1] + output1, B[j-1] + output2, gap_symbol + output3)
	if (i > 0) and (j >= 0) and (k > 0) and (T[i,j,k] == (T[i-1, j, k-1] + SP(A[i-1], gap_symbol, C[k-1]))):
		return backtrack(i-1, j, k-1, A[i-1] + output1, gap_symbol + output2, C[k-1] + output3)
	if (i >= 0) and (j > 0) and (k > 0) and (T[i,j, k] == (T[i, j-1, k-1] + SP(gap_symbol, B[j-1], C[k-1]))):
		return backtrack(i, j-1, k-1, gap_symbol + output1, B[j-1] + output2, C[k-1] + output3)			
	if (i > 0) and (j >= 0) and (k >= 0) and (T[i,j, k] == (T[i-1, j, k] + SP(A[i-1], gap_symbol, gap_symbol))):
		return backtrack(i-1, j, k, A[i-1] + output1, gap_symbol + output2, gap_symbol + output3)		
	if (i >= 0) and (j > 0) and (k >= 0) and (T[i,j, k] == (T[i, j-1, k] + SP(gap_symbol, B[j-1], gap_symbol))):
		return backtrack(i, j-1, k, gap_symbol + output1, B[j-1] + output2, gap_symbol + output3)		
	if (i >= 0) and (j >= 0) and (k > 0) and (T[i,j, k] == (T[i, j, k-1] + SP(gap_symbol, gap_symbol, C[k-1]))):
		return backtrack(i, j, k-1, gap_symbol + output1, gap_symbol + output2, C[k-1] + output3)		
	return output1 + "\n" + output2 + "\n" + output3

def compute_score(s_mat, gc, sequences):
	assert len(sequences)==3
	global A
	global B
	global C
	global score
	global gap_cost
	global T
	score = s_mat
	gap_cost = int(gc)
	A = sequences[0].replace(" ", "")
	B = sequences[1].replace(" ", "")
	C = sequences[2].replace(" ", "")
	n = len(A)
	m = len(B)
	l = len(C)
	T = np.empty([n+1, m+1, l+1])

	T[:] = float("inf")
	calc_T(n, m, l)
	return T[m,n,l]
	
if __name__ == "__main__":
    args = sys.argv

    score, _ = util.read_score_matrix_and_alphabet(args[1])
    gap_cost = int(args[2])
    fastaDictionary = util.read_fasta_file(args[3])
    sequences = iter(fastaDictionary.values())

    A = next(sequences).replace(" ", "")
    B = next(sequences).replace(" ", "")
    C = next(sequences).replace(" ", "")
    # Length of the strings plus 1
    n = len(A)
    m = len(B)
    l = len(C)
	

    T = np.empty([n+1, m+1, l+1])
    T[:] = float("inf")
    calc_T(n, m, l)
	
    print("Score of optimal multiple alignment: \n{}".format(T[n,m,l]))
    print("Optimal multiple alignment: \n{}".format(backtrack(n,m,l,"","","")))