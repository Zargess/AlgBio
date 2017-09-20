import numpy as np
import os, sys
import os.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util

A = ""
B = ""
score_matrix = {}
alpha = 0
beta = 0
alphabet = []

S = None
D = None 
I = None

def cost_S(i, j):
    if S[i,j] == float("inf"):
        v1 = float("inf")
        v2 = float("inf")
        v3 = float("inf")
        v4 = float("inf")

        if i > 0 and j > 0:
            v1 = cost_S(i-1, j-1) + score_matrix[(A[i-1], B[j-1])]    
        if i > 0 and j >= 0:
            v2 = cost_D(i,j)
        if i >= 0 and j > 0:
            v3 = cost_I(i,j)
        if i == 0 and j == 0:
            v4 = 0

        S[i,j] = min(v1, v2, v3, v4)
    return S[i,j]
	
def cost_D(i, j):
    if D[i,j] == float("inf"):
        v1 = float("inf")
        v2 = float("inf")

        if i > 0 and j >= 0:
            v1 = cost_S(i-1,j) + (alpha + beta)
        if i > 1 and j >= 0:
            v2 = cost_D(i-1,j) + alpha 
        D[i,j] = min(v1,v2)
    return D[i,j]

def cost_I(i, j):
    if I[i,j] == float("inf"):
        v1 = float("inf")
        v2 = float("inf")

        if i >= 0 and j > 0:
            v1 = cost_S(i,j-1) + (alpha + beta)
        if i >= 0 and j > 1:
            v2 = cost_I(i,j-1) + alpha 
        I[i,j] = min(v1,v2)
    return I[i,j]

def backtrack(i, j, output1, output2):
    if (i > 0) and (j >= 0) and (S[i,j] == D[i,j]):
        return backtrack(i-1, j, A[i-1] + output1, "-" + output2)
    if (i >= 0) and (j > 0) and (S[i,j] == I[i,j]):
        return backtrack(i, j-1, "-" + output1, B[j-1] + output2)
    if (i > 0) and (j > 0) and (S[i,j] == (S[i-1, j-1] + score_matrix[A[i-1], B[j-1]])):
        return backtrack(i-1, j-1, A[i-1] + output1, B[j-1] + output2)

    return output1 + "\n" + output2

def parse_arguments(args):
    score_matrix_file = args[1]
    a = int(args[2])
    b = int(args[3])
    file1 = args[4]
    file2 = args[5]
    should_output_allignment = int(args[6])
    score_matrix,alphabet = util.read_score_matrix_and_alphabet(score_matrix_file)
    fastaSeq1 = ""
    fastaSeq2 = ""

    if os.path.isfile(file1):
        fastaSeq1 = util.read_fasta_file(file1)
    else:
        fastaSeq1 = {"Seq1" : file1.upper()}

    if os.path.isfile(file2):
        fastaSeq2 = util.read_fasta_file(file2)
    else:
        fastaSeq2 = {"Seq2" : file2.upper()}

    return score_matrix, a, b, should_output_allignment, alphabet, fastaSeq1, fastaSeq2

if __name__ == "__main__":
    args = sys.argv
    score_matrix, alpha, beta, should_output_allignment, alphabet, fastaSeq1, fastaSeq2 = parse_arguments(args)

    A = fastaSeq1["Seq1"].replace(" ", "")
    B = fastaSeq2["Seq2"].replace(" ", "")
    n = len(A)
    m = len(B)

    S = np.empty([n+1, m+1])
    D = np.empty([n+1, m+1])
    I = np.empty([n+1, m+1])
    S[:] = float("inf")
    D[:] = float("inf")
    I[:] = float("inf")

    res = cost_S(n,m)
    if should_output_allignment == 1:
        allignment = backtrack(n, m, "", "").split("\n")
        print(">seq1\n" + allignment[0] + "\n")
        print(">seq2\n" + allignment[1] + "\n")

    print("optimal cost for the two sequences is: \n" + str(res) + "\n")