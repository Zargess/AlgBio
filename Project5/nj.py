import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import util
import newick


id_map = {}
id_map_rev = {}
S = {}

def compute_avg_row_sum(disMatrix, i):
    res = 0
    for m in range(0, disMatrix.shape[1]):
        if(disMatrix[i, m] != float("inf")):
            res += disMatrix[i, m]
    res /= (len(S) - 2)
    return res

def compute_idx_min_N_entry(disMatrix):
    N = np.empty(disMatrix.shape)
    for i in range(disMatrix.shape[0]):
        for j in range(disMatrix.shape[1]):
            N[i,j] = disMatrix[i, j] - (compute_avg_row_sum(disMatrix, i) + compute_avg_row_sum(disMatrix, j))
    np.fill_diagonal(N, float("inf"))
    return np.unravel_index(N.argmin(), N.shape)

def rename_recursive(node):
    if(node.name != "-1"):
        node.name = id_map_rev[int(node.name)]
    for child in node.descendants:
        rename_recursive(child)


if __name__ == "__main__":
    if (len(sys.argv) != 2):
        exit(-1)

    filename = sys.argv[1]
    #filename = "example_slide4.phy"
    dict_matrix, taxa = util.read_score_matrix_and_alphabet(filename)
    disMatrix = np.empty((len(taxa), len(taxa)))

    # Mapping between symbols and internal IDs.
    id = 0
    for taxon in taxa:
        id_map[taxon] = id
        id_map_rev[id] = taxon
        node = newick.Node.create(name=str(id), length=str(float('inf')))
        S[id] = node
        id += 1

    for entry in dict_matrix:
        i = id_map[entry[0]]
        j = id_map[entry[1]]
        disMatrix[i, j] = dict_matrix[entry]


    while len(S) > 3:
        (i, j) = compute_idx_min_N_entry(disMatrix)
        iRowSum = compute_avg_row_sum(disMatrix, i)
        jRowSum = compute_avg_row_sum(disMatrix, j)

        iNode = S[i]
        iNode.length = (0.5 * (disMatrix[i, j] + iRowSum - jRowSum))
        jNode = S[j]
        jNode.length = (0.5 * (disMatrix[i, j] + jRowSum - iRowSum))
        newID = disMatrix.shape[0]
        newNode = newick.Node.create(name=str(newID), descendants=[iNode, jNode])
        S.pop(i)
        S.pop(j)
        newName = str(id_map_rev[i]) + str(id_map_rev[j])
        id_map[newName] = newID
        id_map_rev[newID] = newName

        S[newID] = newNode

        disMatrix = np.insert(disMatrix, disMatrix.shape[0], 0, axis=0)
        disMatrix = np.insert(disMatrix, disMatrix.shape[1], 0, axis=1)

        for m in range(disMatrix.shape[0] - 1):
            if m != i and m != j:
                disMatrix[disMatrix.shape[0] - 1, m] = 0.5 * (disMatrix[i, m] + disMatrix[j, m] - disMatrix[i, j])
                disMatrix[m, disMatrix.shape[0] - 1] = 0.5 * (disMatrix[i, m] + disMatrix[j, m] - disMatrix[i, j])

        disMatrix[i, :] = float("inf")
        disMatrix[j, :] = float("inf")
        disMatrix[:, i] = float("inf")
        disMatrix[:, j] = float("inf")



    keys = list(S.keys())
    dim = disMatrix[keys[0], keys[2]]
    dij = disMatrix[keys[0], keys[1]]
    djm = disMatrix[keys[1], keys[2]]
    gamma_v_i = 0.5 * (dij + dim - djm)
    gamma_v_j = 0.5 * (dij + djm - dim)
    gamma_v_m = 0.5 * (dim + djm - dij)

    S[keys[0]].length = gamma_v_i
    S[keys[1]].length = gamma_v_j
    S[keys[2]].length = gamma_v_m

    S[-1] = newick.Node.create(name="-1", descendants=[S[keys[0]], S[keys[1]], S[keys[2]]])

    rename_recursive(S[-1])
    print(newick.dumps(S[-1]))








