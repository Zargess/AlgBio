import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import util
import newick


id_map = {}
id_map_rev = {}
S = {}
shape = None

def compute_avg_row_sum(disMatrix, i):
    res = np.nansum(disMatrix[i,:])
    res /= (len(S) - 2)
    return res

def compute_idx_min_N_entry(disMatrix):
    bestValue = float("inf")
    bestTuple = None
    rowSumMap = {}
    for i in range(shape[0]):
        rowSumMap[i] = compute_avg_row_sum(disMatrix, i)

    for i in range(shape[0]):
        for j in range(shape[1]):
            if(i == j):
                continue
            value = disMatrix[i, j] - (rowSumMap[i] + rowSumMap[j])
            if(value < bestValue):
                bestValue = value
                bestTuple = (i, j, rowSumMap[i], rowSumMap[j])
    return bestTuple

def rename_recursive(node):
    if(node.name != "root"):
        node.name = id_map_rev[int(node.name)]
    for child in node.descendants:
        rename_recursive(child)


if __name__ == "__main__":
    global shape
    if (len(sys.argv) != 2):
        exit(-1)

    filename = sys.argv[1]
    dict_matrix, taxa = util.read_score_matrix_and_alphabet(filename)
    disMatrix = np.full((2*len(taxa), 2*len(taxa)), np.nan)
    shape = (len(taxa) , len(taxa))


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
        #print("Starting iteration with |S| = {}".format(len(S)))
        (i, j, iRowSum, jRowSum) = compute_idx_min_N_entry(disMatrix)

        iNode = S[i]
        iNode.length = (0.5 * (disMatrix[i, j] + iRowSum - jRowSum))
        jNode = S[j]
        jNode.length = (0.5 * (disMatrix[i, j] + jRowSum - iRowSum))
        newID = shape[0]
        newNode = newick.Node.create(name=str(newID), descendants=[iNode, jNode])
        S.pop(i)
        S.pop(j)
        newName = str(id_map_rev[i]) + str(id_map_rev[j])
        id_map[newName] = newID
        id_map_rev[newID] = newName

        S[newID] = newNode


        for m in range(shape[0]):
            if m != i and m != j:
                newDist = 0.5 * (disMatrix[i, m] + disMatrix[j, m] - disMatrix[i, j])
                disMatrix[shape[0], m] = newDist
                disMatrix[m, shape[0]] = newDist
        shape = (shape[0] + 1, shape[1] + 1)

        disMatrix[i, :] = np.nan #float("inf")
        disMatrix[j, :] = np.nan #float("inf")
        disMatrix[:, i] = np.nan #float("inf")
        disMatrix[:, j] = np.nan #float("inf")



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

    S["root"] = newick.Node.create(name="root", descendants=[S[keys[0]], S[keys[1]], S[keys[2]]])

    rename_recursive(S["root"])
    print(newick.dumps(S["root"]))








