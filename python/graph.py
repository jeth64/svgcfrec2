

"""
maeks a set of edges where the order of vertices is irrelevant
"""
def makeSet(edges):
   return set(map(lambda e: tuple(sort(e)), edges))


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
      if e[0] in edgeMap:
         edgeMap[e[1]].add(e[0])
      else:
         edgeMap[e[1]] = set([e[0]])
  return edgeMap


"""
Return nodes of an undirected graph where f(degree) is true
"""
def getNodes(edgeMap, f):
   return [key for (key,val) in edgeMap.iteritems() if f(len(val))]
