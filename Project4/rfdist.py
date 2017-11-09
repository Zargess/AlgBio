import sys
import io
import newick
import math
import random
import time
import itertools
import os

class TreeNode:
    def __init__(self, id, children):
        self.isRootOrSubRoot = False
        self.children = children
        self.id = id
        self.interval = None
        self.subtree_size = 0

    def isLeaf(self):
        return len(self.children) == 0

    def addChild(self, child):
        self.children.append(child)

counter = 0
dfsdict = {}

def pp_tree(tree : TreeNode, lvl):
    filler = '----'
    for i in range(0, lvl):
        sys.stdout.write(filler)

    #if(tree.isLeaf()):
    #    print("ID: {} DFS_NUM: {}".format(tree.id, dfsdict[tree.id]))
    if(not tree.isLeaf()):
        print("ID: {} INTERVAL: {} SUBTREESIZE: {}".format(tree.id, tree.interval, tree.subtree_size))
    #if (not tree.isLeaf()):
    #    print("ID: {} SUBTREESIZE: {}".format(tree.id, tree.subtree_size))
    for child in tree.children:
        pp_tree(child, lvl+1)

def count_leaves(tree : TreeNode):
    if(tree.isLeaf()):
        #print("Hello I am leaf with name: {}".format(tree.id))
        return 1
    else:
        counter = 0
        for child in tree.children:
            counter += count_leaves(child)
        return counter


def dfs_numbering(tree : TreeNode, currentid):
    if (tree.isLeaf()):
        dfsdict[tree.id] = currentid
        return currentid + 1

    for child in tree.children:
        currentid = dfs_numbering(child, currentid)

    return currentid

def annotate_tree_intervals(tree: TreeNode):
    if(tree.isLeaf()):
        dfs_num = dfsdict[tree.id]
        tree.interval = (dfs_num, dfs_num)
        return

    min = float('inf')
    max = -float('inf')
    for child in tree.children:
        annotate_tree_intervals(child)

        childmin = child.interval[0]
        childmax = child.interval[1]

        if(childmin < min):
            min = childmin
        if(childmax > max):
            max = childmax
    tree.interval = (min, max)

def annotate_tree_subsetsize(tree: TreeNode):
    for child in tree.children:
        annotate_tree_subsetsize(child)
        if(child.isLeaf()):
            tree.subtree_size += 1
        else:
            tree.subtree_size += child.subtree_size

def report_intervals(tree: TreeNode, reportOnlyCandidates):
    if(tree.isLeaf()):
        return []

    res = []
    isCandidateInterval = tree.interval[1] - tree.interval[0] + 1 == tree.subtree_size
    if(not tree.isRootOrSubRoot and (not reportOnlyCandidates or isCandidateInterval)):
        res.append(tree.interval)

    for child in tree.children:
        if not child.isLeaf() :
            childintervals = report_intervals(child, reportOnlyCandidates)
            res.extend(childintervals)
    return res

def radix_sort_intervals(intervals, n):
    for k in reversed(range(2)):
        buckets = [[] for _ in range(n)]
        for interval in intervals:
            buckets[interval[k]].append(interval)
        intervals = [i for bucket in buckets for i in bucket]
    return intervals

def intervals_equal(int1, int2):
    return int1[0] == int2[0] and int1[1] == int2[1]

def compare_intervals_le(int1, int2):
    if(int1[0] == int2[0]):
        return int1[1] < int2[1]
    else:
        return int1[0] < int2[0]


    return int1[0] < int2[0] and int1[1] < int2[1]

def count_shared_intervals(intervalsT1, intervalsT2):
    #This methods assumes they are sorted and intervalsT2 only contains candidate intervals.
    total = 0
    i = 0
    j = 0
    l1 =len(intervalsT1)
    l2 = len(intervalsT2)
    while(i < l1 and j < l2):
        int_t1 = intervalsT1[i]
        int_t2 = intervalsT2[j]
        if(intervals_equal(int_t1, int_t2)):
            i += 1
            j += 1
            total += 1
        elif(compare_intervals_le(int_t1, int_t2)):
            i += 1
        else:
            j += 1

    return total

def generate_test_trees():
    # T1
    t4 = TreeNode(4, [])
    t3 = TreeNode(3, [])
    t43 = TreeNode(43, [t4, t3])

    t5 = TreeNode(5, [])
    t543 = TreeNode(543, [t5, t43])

    t7 = TreeNode(7, [])
    t7543 = TreeNode(7543, [t7, t543])

    t6 = TreeNode(6, [])
    t2 = TreeNode(2, [])
    t26 = TreeNode(62, [t6, t2])

    t_all = TreeNode(0, [t26, t7543])
    t_all.isRootOrSubRoot = True
    t1_root = TreeNode(1, [t_all])
    t1_root.isRootOrSubRoot = True

    # T2
    t7 = TreeNode(6, [])
    t6 = TreeNode(7, [])
    t76 = TreeNode(76, [t7, t6])

    t5 = TreeNode(5, [])
    t576 = TreeNode(576, [t5, t76])

    t2 = TreeNode(2, [])
    t2576 = TreeNode(2576, [t2, t576])

    t3 = TreeNode(3, [])
    t4 = TreeNode(4, [])
    t34 = TreeNode(34, [t3, t4])

    t_all = TreeNode(0, [t2576, t34])
    t_all.isRootOrSubRoot = True
    t2_root = TreeNode(1, [t_all])
    t2_root.isRootOrSubRoot = True

    return t1_root, t2_root

def create_children(descendants):
    children = []
    for child in descendants:
        children.append(TreeNode(child.name, create_children(child.descendants)))
    return children

def parse_newick_to_tree(filename):
    with io.open(filename, encoding='utf8') as fp:
        trees = newick.load(fp)

    tree_root = TreeNode(trees[0].name, [])
    for child in create_children(trees[0].descendants):
        tree_root.addChild(child)
    return tree_root

def compute_rf_distance(t1_root, t2_root):
    # DFS number T1
    dfs_numbering(t1_root, 1)
    annotate_tree_intervals(t1_root)

    # Annotate T2 with numbers and subtreesizes
    annotate_tree_intervals(t2_root)
    annotate_tree_subsetsize(t2_root)

    # T1 Intervals
    t1ints = report_intervals(t1_root, False)
    n = t1_root.interval[1] + 1

    # T2 candidates and ints
    t2ints = report_intervals(t2_root, False)
    t2_cand_ints = report_intervals(t2_root, True)

    t1ints = radix_sort_intervals(t1ints, n)
    t2_cand_ints = radix_sort_intervals(t2_cand_ints, n)

    numShared = count_shared_intervals(t1ints, t2_cand_ints)

    return (len(t1ints) - numShared) + (len(t2ints) - numShared)

def generate_tree(numLeaves):
    root = TreeNode("id" + str(numLeaves), [])
    if(numLeaves == 1):
        return root
    elif(numLeaves >= 3):
        coin = random.randint(0, 1)
        if(coin):
            athird = (numLeaves / 3)
            if(numLeaves% 3 == 0):
                root.addChild(generate_tree(athird))
                root.addChild(generate_tree(athird))
                root.addChild(generate_tree(athird))
            elif (numLeaves % 3 == 1):
                root.addChild(generate_tree(math.floor(athird)))
                root.addChild(generate_tree(math.floor(athird)))
                root.addChild(generate_tree(math.ceil(athird)))
            elif (numLeaves % 3 == 2):
                root.addChild(generate_tree(math.floor(athird)))
                root.addChild(generate_tree(math.ceil(athird)))
                root.addChild(generate_tree(math.ceil(athird)))
        else:
            half = numLeaves / 2
            root.addChild(generate_tree(math.floor(half)))
            root.addChild(generate_tree(math.ceil(half)))

    elif (numLeaves >= 2):
        half = numLeaves / 2
        root.addChild(generate_tree(math.floor(half)))
        root.addChild(generate_tree(math.ceil(half)))
    #if (numLeaves >= 1):
        #root.addChild(generate_tree(numLeaves-1))
    return root


if __name__ == "__main__":

    results = []
    n = 1000
    for i in range(20):
        endSum = 0
        for j in range(5):
            tree1 = generate_tree(n)
            tree2 = generate_tree(n)
            start = time.clock()
            dist = compute_rf_distance(tree1, tree2)
            end = time.clock() - start

            endSum += end

        print("Completed iteration: {}".format(i))
        results.append((n, endSum / 5))
        n = int(n * 1.3)

    with open("output.txt", "w+") as fp:
        for (i, time) in results:
            fp.write(str(i) + "," + str(time) + "\n")
    '''
    folder = sys.argv[1]
    folder2 = sys.argv[2]
    allFiles = os.listdir(folder)
    allFilesPerm = os.listdir(folder2)


    combs = itertools.combinations(allFiles, 2)
    print(combs)
    i = 0
    for (file, filePerm) in zip(allFiles, allFilesPerm):
        #print((file,filePerm))
        t1 = parse_newick_to_tree(folder + "\\" + file)
        t2 = parse_newick_to_tree(folder2 + "\\" + filePerm)

        #print("Leaves in T1: {} Leaves in T2: {}".format(count_leaves(t1), count_leaves(t2)))
        dist = compute_rf_distance(t1, t2)

        print("Distance between: {} and {} is: {}.".format(file, filePerm, dist))
    
    for comb in combs:
        i += 1
        t1 = parse_newick_to_tree(folder + "\\" + comb[0])
        t2 = parse_newick_to_tree(folder + "\\" + comb[1])
        dist = compute_rf_distance(t1, t2)
        dist2 = compute_rf_distance(t2, t1)
        equal = dist == dist2
        print("Distance between: {} and {} is: {}. Symmetric? {}".format(comb[0], comb[1], dist, equal))
    print(i)
    #with open(folder + "\\" + allFiles[0]) as fp:

    #print(allFiles)

    t1_root, t2_root = generate_test_trees()
    t1_root.isRootOrSubRoot = True
    t2_root.isRootOrSubRoot = True
    leafcount = count_leaves(t2_root)

    #TODO: Root trees if they are not rooted
    t1 = parse_newick_to_tree('Testdata/tree1.new')

    t2 = parse_newick_to_tree('Testdata/tree2.new')
    '''




