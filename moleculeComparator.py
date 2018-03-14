# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 16:03:28 2018

@author: michal
"""

from fakeBiopython import myResidue
from cif_analyser import molecule2graph
from networkx.algorithms.isomorphism import GraphMatcher

class moleculeMatcher(GraphMatcher):
    def semantic_feasibility(self, G1_node, G2_node):
        return self.G1.node[G1_node]["element"] == self.G2.node[G2_node]["element"]
    

res1 = myResidue()
res1.read_xyz("xyz/isomorphism/szczawiowy.xyz")
res1G = molecule2graph(res1.get_atoms())

res2 = myResidue()
res2.read_xyz("xyz/isomorphism/szczawiowy_partial.xyz")
res2G = molecule2graph(res2.get_atoms())

gMatcher = moleculeMatcher( res1G, res2G)
print(gMatcher.is_isomorphic())
print(gMatcher.mapping)

print(gMatcher.subgraph_is_isomorphic())
print(gMatcher.mapping)