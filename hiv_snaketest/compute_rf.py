from ete3 import Tree

# Load the trees from Newick files
tree1 = Tree("true.newick", format=1)
tree2 = Tree("out.newick", format=1)
print(tree1)
print(tree2)

# Compute the Robinson-Foulds distance between the two trees
rf = tree2.robinson_foulds(tree1, unrooted_trees=True)

# Output results
#print(f"Robinson-Foulds distance: {rf_distance}")
#print(f"Maximum possible RF distance: {max_rf}")
#print(f"Common leaves: {len(common_leaves)}")
#print(f"Partitions in tree1: {len(parts_t1)}")
#print(f"Partitions in tree2: {len(parts_t2)}")
print(rf)

