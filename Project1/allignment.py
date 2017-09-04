import numpy as np
import sys

score = np.matrix([
    [10, 2, 5, 2],
    [2, 10, 2, 5],
    [5, 2, 10, 2],
    [2, 5, 2, 10]
])

gap_cost = -5

T = None

A = ""
B = ""

def getIndex(string, i):
    symbol = string[i-1]
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    return 3
    

def cost(i, j):
    if T[i,j] == float("inf"):
        v1 = -float("inf")
        v2 = -float("inf")
        v3 = -float("inf")
        v4 = -float("inf")

        if i > 0 and j > 0:
            v1 = cost(i-1, j-1) + score[getIndex(A,i), getIndex(B,j)]
        if i > 0 and j >= 0:
            v2 = cost(i-1, j) + gap_cost
        if i >= 0 and j > 0:
            v3 = cost(i, j-1) + gap_cost
        if i == 0 and j == 0:
            v4 = 0

        T[i,j] = max(v1, v2, v3, v4)
    return T[i,j]

if __name__ == "__main__":
    args = sys.argv
    file1 = args[1]
    file2 = args[2]

    #TODO : Read fasta files and get lengths of strings
    A = "TCCAGAGA"
    B = "TCGAT"
    n = len(A)
    m = len(B)

    T = np.empty([n+1, m+1])
    T[:] = float("inf")
    res = cost(n,m)
    print(res)
