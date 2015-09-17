import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import depth_first_tree

from tools import *

"""
Makes a set of edges where the order of vertices is irrelevant
"""
def makeSet(edges):
   return set(map(lambda e: tuple(sorted(e)), edges))


"""
Returns edgeMap for graph algorithms
"""
def getEdgeMap(edges):
   edgeMap = dict([])
   for e in makeSet(edges):
      if e[0] in edgeMap:
         edgeMap[e[0]].add(e[1])
      else:
         edgeMap[e[0]] = set([e[1]])
      if e[1] in edgeMap:
         edgeMap[e[1]].add(e[0])
      else:
         edgeMap[e[1]] = set([e[0]])
   return edgeMap


"""
Return nodes of an undirected graph where f(degree) is true
"""
def getNodes(edgeMap, f):
   return [key for (key,val) in edgeMap.iteritems() if f(len(val))]

"""
Returns adjacency matrix of undirected graph with given weights
"""
def getUndirAdjMatrix(edges, N, weights=None):
   if weights is None: weights = np.ones((len(edges), 1))
   data = np.reshape(np.vstack([weights, weights]),(2*len(edges),))
   ij = np.transpose(np.vstack([edges, np.roll(edges, 1, 1)])) if len(edges)>0 else ([], [])
   return csr_matrix((data, ij), (N, N))


def traceBack(pred, start, connLimit = 20):
   node = start
   path = [] 
   i=0
   while i < connLimit+1 and i < 1000+1:
      i = i+1
      previous = pred[node]
      if previous > 0: # changed from from isfinite
         path.append(int(previous))
         if previous == start:
            break
         else:
            node = previous
      else:
         print "Warning: Tree path tracing has reached dead end"
         break
   else: print "Warning: Tree path tracing has timed out!"
   return path

"""
Recursion for depth-limited search
"""
def dlsRec(weightedGraph, pred, stopNode, node, pathlist, curLimit, connLimit, lengthLimit):
   # side-effect on pathlist
   if curLimit > 0:
      for v in weightedGraph.getrow(node).nonzero()[1]:
         remainingLength = lengthLimit - weightedGraph._get_single_element(node, v)
         if not (v == pred[node] or remainingLength < 0): # do not consider node far along the path
            if stopNode == v:
               pred[v] = node
               pathlist.append(traceBack(pred, v, connLimit))
            else:
               if pred[v] ==-1:
                  pred[v] = node
                  dlsRec(weightedGraph, pred, stopNode, v, pathlist, curLimit-1, connLimit, remainingLength)
                  pred[v] = -1
   return

"""
Return lists of paths describing cycles
"""
def getCycles(edges, vD, connLimit=np.inf, lengthLimit=np.inf):
   if len(edges)>0:
      N = len(vD.vertices)
      dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
      graph = getUndirAdjMatrix(edges, N, dists)
      dfsTree = depth_first_tree(graph, edges[0][0], False)
      cycleHints = makeSet(edges).difference(makeSet(np.transpose(dfsTree.nonzero())))
      paths = []
      for e in cycleHints:
         dlsRec(graph, np.ones((N,)) * -1, e[0], e[0], paths, connLimit, connLimit, lengthLimit)
      return repulseCycles(paths), cycleHints #map(lambda x: x[:-1], repulseCycles(paths)), cycleEdges
   else: return [], []

"""
Return cycle Edges from vertex list
"""
def getCycleEdges(vertices):
   return map(lambda x,y: tuple(sorted([x,y])), vertices, np.roll(vertices,1))

"""
Return list of unique cycles
"""
def repulseCycles(cycles):
   return map(lambda x: x[0], groupbyKeys(map(lambda c: tuple(sorted(getCycleEdges(c))), cycles), cycles))





