import numpy as np
import os, sys
import os.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util

A = ""
B = ""
score_matrix = {}
gap_cost = 0
alphabet = []
T = None

def cost(i, j):
    if T[i,j] == float("inf"):
        v1 = -float("inf")
        v2 = -float("inf")
        v3 = -float("inf")
        v4 = -float("inf")

        if i > 0 and j > 0:
            v1 = cost(i-1, j-1) + score_matrix[(A[i-1], B[j-1])]
        if i > 0 and j >= 0:
            v2 = cost(i-1, j) + gap_cost
        if i >= 0 and j > 0:
            v3 = cost(i, j-1) + gap_cost
        if i == 0 and j == 0:
            v4 = 0

        T[i,j] = max(v1, v2, v3, v4)
    return T[i,j]

def backtrack(i, j, output1, output2):
    if (i > 0) and (j > 0) and (T[i,j] == (T[i-1, j-1] + score_matrix[A[i-1], B[j-1]])):
        return backtrack(i-1, j-1, A[i-1] + output1, B[j-1] + output2)
    if (i > 0) and (j >= 0) and (T[i,j] == T[i-1, j] + gap_cost):
        return backtrack(i-1, j, A[i-1] + output1, "-" + output2)
    if (i >= 0) and (j > 0) and (T[i,j] == T[i,j-1] + gap_cost):
        return backtrack(i, j-1, "-" + output1, B[j-1] + output2)

    return output1 + "\n" + output2

if __name__ == "__main__":
    args = sys.argv
    score_matrix, gap_cost, should_output_allignment, alphabet, fastaSeq1, fastaSeq2 = util.parse_arguments(args)

    A = fastaSeq1["Seq1"].replace(" ", "")
    B = fastaSeq2["Seq2"].replace(" ", "")
    n = len(A)
    m = len(B)

    T = np.empty([n+1, m+1])
    T[:] = float("inf")
    res = cost(n,m)
    # TODO: If should_output_allignment == 1 -> print an optimal allignment
    if should_output_allignment == 1:
        print("an optimal alignment is: \n" + backtrack(n, m, "", "") + "\n")
    print("optimal cost for the two sequences is: \n" + str(res) + "\n")