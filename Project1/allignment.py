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

def read_fasta_file(filename):
    """
    Reads the given FASTA file f and returns a dictionary of sequences.

    Lines starting with ';' in the FASTA file are ignored.
    """
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences

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

def backtrack(i, j, output1, output2):
    if (i > 0) and (j > 0) and (T[i,j] == (T[i-1, j-1] + score[getIndex(A,i), getIndex(B,j)])):
        return backtrack(i-1, j-1, A[i-1] + output1, B[j-1] + output2)
    if (i > 0) and (j >= 0) and (T[i,j] == T[i-1, j] + gap_cost):
        return backtrack(i-1, j, A[i-1] + output1, "-" + output2)
    if (i >= 0) and (j > 0) and (T[i,j] == T[i,j-1] + gap_cost):
        return backtrack(i, j-1, "-" + output1, B[j-1] + output2)

    return output1 + "\n" + output2

def count(i, j, counter):
    if (i > 0) and (j > 0) and (T[i,j] == (T[i-1, j-1] + score[getIndex(A,i), getIndex(B,j)])):
        counter = count(i-1, j-1, counter)
    if (i > 0) and (j >= 0) and (T[i,j] == T[i-1, j] + gap_cost):
        counter = count(i-1, j, counter)
    if (i >= 0) and (j > 0) and (T[i,j] == T[i,j-1] + gap_cost):
        counter = count(i, j-1, counter)
    if(i == 0) and (j == 0):
        counter = counter + 1
    return counter
    
if __name__ == "__main__":
    args = sys.argv
    file1 = args[1]
    file2 = args[2]

    fastaSeq1 = read_fasta_file(file1)
    fastaSeq2 = read_fasta_file(file2)
    

    #TODO : Read fasta files and get lengths of strings
    #A = "AATAAT" 
    #B = "AAGG" 
    A = fastaSeq1["Seq1"].replace(" ", "")
    B = fastaSeq2["Seq2"].replace(" ", "")
    n = len(A)
    m = len(B)

    T = np.empty([n+1, m+1])
    T[:] = float("inf")
    res = cost(n,m)
    counter = 0
    print("optimal cost for seqs: " + str(A) + " and " + str(B) + " is: \n" + str(res) + "\n")
    print("an optimal alignment is: \n" + backtrack(n, m, "", "") + "\n")

    print("number of possible optimal alignments is: \n" + str(count(n, m, counter)))
    
