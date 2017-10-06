import numpy as np
import os, sys
import os.path
import time
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util



A = ""
B = ""
score_matrix = {}
gap_cost = 0
alphabet = []
T = None

def score_func(key):
	if key in score_matrix:
		return score_matrix[key]
	else:
		return 0

def cost_iter(n,m):
	for i in range(0, n+1):
		for j in range(0, m+1):
			v1 = float("inf")
			v2 = float("inf")
			v3 = float("inf")
			v4 = float("inf")

			if i > 0 and j > 0:
				v1 = T[i-1, j-1] + score_func((A[i-1], B[j-1]))
			if i > 0 and j >= 0:
				v2 = T[i-1, j] + int(gap_cost)
			if i >= 0 and j > 0:
				v3 = T[i, j-1] + int(gap_cost)
			if i == 0 and j == 0:
				v4 = 0

			T[i, j] = min(v1, v2, v3, v4)
	return T[i,j]

def cost(i, j):
    if T[i,j] == float("inf"):
        v1 = float("inf")
        v2 = float("inf")
        v3 = float("inf")
        v4 = float("inf")
        
        if i > 0 and j > 0:
            v1 = cost(i-1, j-1) + score_matrix[(A[i-1], B[j-1])]
        if i > 0 and j >= 0:
            v2 = cost(i-1, j) + int(gap_cost)
        if i >= 0 and j > 0:
            v3 = cost(i, j-1) + int(gap_cost)
        if i == 0 and j == 0:
            v4 = 0

        T[i,j] = min(v1, v2, v3, v4)
    return T[i,j]

def backtrack(i, j, output1, output2):
    if (i > 0) and (j > 0) and (T[i,j] == (T[i-1, j-1] + score_func((A[i-1], B[j-1])))):
        return backtrack(i-1, j-1, A[i-1] + output1, B[j-1] + output2)
    if (i > 0) and (j >= 0) and (T[i,j] == T[i-1, j] + int(gap_cost)):
        return backtrack(i-1, j, A[i-1] + output1, "-" + output2)
    if (i >= 0) and (j > 0) and (T[i,j] == T[i,j-1] + int(gap_cost)):
        return backtrack(i, j-1, "-" + output1, B[j-1] + output2)

    return output1 + "\n" + output2
	
def backtrack_iter(i,j):
    output1 = ""
    output2 = ""
    while i > 0 or j > 0:
        if (i > 0 and j > 0) and (T[i,j] == T[i-1, j-1] + score_func((A[i-1], B[j-1]))):
            output1 = A[i-1] + output1
            output2 = B[j-1] + output2
            i -= 1
            j -= 1
        elif (i > 0) and (j >= 0) and (T[i,j] == T[i-1, j] + int(gap_cost)):
            output1 = A[i-1] + output1
            output2 = "-" + output2
            i -= 1
        elif (i >= 0) and (j > 0) and (T[i,j] == T[i, j-1] + int(gap_cost)):
            output1 = "-" + output1
            output2 = B[j-1] + output2
            j -= 1
    return output1 + "\n" + output2

def runAlgo(string1, string2, s_mat=score_matrix, gc=gap_cost):
    global A
    global B
    global T
    global score_matrix
    global gap_cost    

    A = string1
    B = string2
    
    n = len(A)
    m = len(B)
    score_matrix = s_mat
    gap_cost = gc
    
    
	
    T = np.empty([n+1, m+1])
    T[:] = float("inf")
    s = cost_iter(n,m)
    return s

def runAlgoWithBacktrack(string1, string2, s_mat=score_matrix, gc=gap_cost):
	runAlgo(string1, string2, s_mat=score_matrix, gc=gap_cost)
	return backtrack_iter(len(string1), len(string2)).split("\n")
	
def run_experiment(startN, iterations):
    length = startN
    set1, set2 = util.generate_data_equal_length(startN, iterations)
    lengths = []
    values = []
    counter = 0
    print("Starting experiment")
    for i in range(0, iterations):

        start = time.time()
        for j in range(0, 5):
            runAlgo(set1[counter], set2[counter])
            counter += 1
        end = time.time()
        print("Iteration: " + str(i) + " is complete")
        lengths.append(length)
        values.append((end - start) / 5)
        length = int(length * 1.3)
    return lengths, values

if __name__ == "__main__":
    args = sys.argv
    score_matrix, gap_cost, should_output_allignment, alphabet, fastaSeq1, fastaSeq2 = util.parse_arguments(args)
    
    A = next(iter(fastaSeq1.values())).replace(" ", "")
    B = next(iter(fastaSeq2.values())).replace(" ", "")
    print(A)
    print(B)
    	
    n = len(A)
    m = len(B)
	
    print(runAlgo(A, B, s_mat=score_matrix, gc=gap_cost))
    #T = np.empty([n+1, m+1])
    #T[:] = float("inf")
    #print(np.dtype(T[0,0]))
    #res = cost(n,m)
    # TODO: If should_output_allignment == 1 -> print an optimal allignment
    if should_output_allignment == 1:
        allignment = backtrack_iter(n, m).split("\n")
        print(">seq1\n" + allignment[0] + "\n")
        print(">seq2\n" + allignment[1] + "\n")
    #print("optimal cost for the two sequences is: \n" + str(res) + "\n")
    
    #lengths, values = run_experiment(10, 15)
    #util.write_experiment_results_to_file("data/global_linear.data", lengths, values)
