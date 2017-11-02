import sys
class TreeNode:
    def __init__(self, id, children):
        self.isRootOrSubRoot = False
        self.children = children
        self.id = id
        self.interval = None
        self.subtree_size = 0

    def isLeaf(self):
        return len(self.children) == 0


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


if __name__ == "__main__":
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

    dfs_numbering(t1_root, 1)
    annotate_tree_intervals(t1_root)


    annotate_tree_intervals(t2_root)
    annotate_tree_subsetsize(t2_root)

    pp_tree(t2_root, 0)
    print('')
    print(report_intervals(t2_root, False))
    print(report_intervals(t2_root, True))


