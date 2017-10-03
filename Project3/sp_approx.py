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


if __name__ == "__main__":
    args = sys.argv

    score, _ = util.read_score_matrix_and_alphabet(args[1])
    fastaDictionary = util.read_fasta_file(args[2])
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