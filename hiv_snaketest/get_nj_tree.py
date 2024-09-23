from Bio import AlignIO
import io
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import sys

# Load the MSA
alignment = AlignIO.read(sys.argv[1], "fasta")

# Calculate the distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Construct the NJ tree
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(distance_matrix)

# Visualize the tree
Phylo.draw(nj_tree)

# Optionally, save the tree to a file
newick_io = io.StringIO()
Phylo.write(nj_tree, newick_io, "newick")
newick_format = newick_io.getvalue()
with open("out.newick" ,"w") as f:
    f.write(newick_format)

