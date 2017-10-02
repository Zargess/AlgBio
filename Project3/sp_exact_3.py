import numpy as np
import os, sys
import os.path
import time
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
    for i in range(0, n):
        for j in range(0, m):
            for k in range(0, l):
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

    print(T[n,m,l])

if __name__ == "__main__":
    args = sys.argv

    score, _ = util.read_score_matrix_and_alphabet(args[1])
    print(score)
    A = "aacgt"
    B = "ggatc"
    C = "aaaaa"
    # Length of the strings plus 1
    n = len(A)
    m = len(B)
    l = len(C)

    T = np.empty([n+1, m+1, l+1])

    T[:] = float("inf")

    calc_T(n, m, l)

    #print(T)