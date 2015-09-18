import numpy as np
from scipy.spatial import Voronoi
from matplotlib.path import Path
from itertools import *
from operator import add
from priodict import priorityDictionary

# own modules
from tools import *
from graph import *
from geometry import *
from wedge import *


"""
Even-odd rule
"""
def pointsInShape(points, mplpaths):
   # compoundShape considers points in holes as in shape, thats why the insideness has to be checked differently
   
   return map(lambda point: np.sum(map(lambda p: p.contains_point(point), mplpaths))%2, points)

"""
Checks if one of the given mplpaths intersect a line
"""   
def lineIntersectsPath(line, mplpaths):
   mplLine = Path(line, [1, 2])
   return np.any(map(lambda p: p.intersects_path(mplLine, False), mplpaths))

"""
Calculate Voronoi skeleton from vertices 
"""
def skeleton(vertexLists, mplpaths): # allows intersections with path
   points = np.vstack(vertexLists)
   voronoiDiagram = Voronoi(points)
   vertexValidity = pointsInShape(voronoiDiagram.vertices, mplpaths)
   ridgeValidity = np.array(map(lambda x: -1 not in x and vertexValidity[x[0]] and vertexValidity[x[1]],  voronoiDiagram.ridge_vertices),bool)
   skeletonEdges = np.array(voronoiDiagram.ridge_vertices)[ridgeValidity]
   skeletonEdges = filter(lambda e: True != (lineIntersectsPath(voronoiDiagram.vertices[e,:], mplpaths)), skeletonEdges) # funktioniert nicht ueberall: eckpunkte
   return voronoiDiagram, skeletonEdges, ridgeValidity

"""
Calculate Voronoi skeleton from vertices 
"""
def skeleton2(vertexLists, mplpaths): # allows intersections with path
   accLengths = reductions(add, imap(len, vertexLists))
   points = np.vstack(vertexLists)
   voronoiDiagram = Voronoi(points)
   vertexValidity = pointsInShape(voronoiDiagram.vertices, mplpaths)
   if len(vertexLists) > 1:
      validPtEdges = set(map(lambda v: tuple(sorted([v-1, v%accLengths[-1]])), accLengths))
   else:
      validPtEdges = set([])
   invalidPtEdges = set(map(lambda v1, v2: tuple(sorted([v1%accLengths[-1], v2-1])), accLengths, np.roll(accLengths, -1)))
   ridgeValidity = np.array(map(lambda e1, e2: (abs(e1[0]-e1[1]) > 1 or tuple(sorted(e1)) in validPtEdges) and (not tuple(sorted(e1)) in invalidPtEdges) \
                                             and (-1 not in e2) and vertexValidity[e2[0]] and vertexValidity[e2[1]] and not lineIntersectsPath(voronoiDiagram.vertices[e2,:], mplpaths), \
                                voronoiDiagram.ridge_points, voronoiDiagram.ridge_vertices),bool) # without vertexValidity check: with exoskeleton
   skeletonEdges = np.array(voronoiDiagram.ridge_vertices)[ridgeValidity]
   return voronoiDiagram, skeletonEdges, ridgeValidity

"""
Splits line if it crosses the shape boundaries
"""
def condSplitLine(nodes, coordinates, mplpaths):
   if len(nodes) > 1:
      if lineIntersectsPath(coordinates[[nodes[0], nodes[-1]], :], mplpaths):
         if len(nodes) > 2:
            edges = condSplitLine(nodes[:(len(nodes)/2)+1], coordinates, mplpaths) + condSplitLine(nodes[len(nodes)/2:], coordinates, mplpaths)
         else: edges = [[]]
      else: edges = [[nodes[0], nodes[-1]]]
   else: edges = [[]]
   return edges

"""
Simplify Voronoi skeleton with splitting after connecting points
"""
def simplifySkeleton(skeletonEdges, coordinates, mplpaths):
   newEdges = []
   if len(skeletonEdges) > 0:
      edgeMap = getEdgeMap(skeletonEdges)
      leaves = getNodes(edgeMap, lambda deg: deg<2)
      joints = getNodes(edgeMap, lambda deg: deg>2)
      endpoints = set(joints).union(leaves)
      for endpoint in joints:
         if endpoint in edgeMap:
            neighbours = edgeMap.pop(endpoint)
            for n in neighbours:
               node = endpoint
               curNode = n
               if curNode in edgeMap and not lineIntersectsPath(coordinates[[node, curNode], :], mplpaths):
                  prevNode = endpoint
                  nodes = [prevNode, curNode] #new
                  while curNode not in endpoints:
                     cand = edgeMap.pop(curNode).difference([prevNode]).pop()
                     prevNode = curNode
                     curNode = cand
                     nodes.append(curNode) #new
                  edges = condSplitLine(nodes, coordinates, mplpaths)
                  newEdges.extend(edges)
   return newEdges
   

def Dijkstra(G,start,end=None):

   D = {}   # dictionary of final distances
   P = {}   # dictionary of predecessors
   Q = priorityDictionary()   # est.dist. of non-final vert.
   Q[start] = 0

   for v in Q: 
      D[v] = Q[v]  
      if v == end: 
         break
      for w in G.getrow(v).nonzero()[1]:
         vwLength = D[v] + G[v,w]
         if w in D:
            if vwLength < D[w]:
               raise ValueError, \
  "Dijkstra: found better path to already-final vertex"
         elif w not in Q or vwLength < Q[w]:
            Q[w] = vwLength
            P[w] = v
   
   return (D,P)
         
def shortestPath(G,start,end):
   """
   Find a single shortest path from the given start vertex
   to the given end vertex.
   The input has the same conventions as Dijkstra().
   The output is a list of the vertices in order along
   the shortest path.
   """

   D,P = Dijkstra(G,start,end)
   Path = []
   while 1:
      Path.append(end)
      if end == start: break
      end = P[end]
   Path.reverse()
   return Path

"""
Returns multiple possible arms
"""
def traceArm(edgeMap, traceStart, vD, dir1, dir2, paths, shortestPaths=None, G=None, lineStart=None):
   if lineStart is None:
      lineStart = traceStart
   dirvec = normalize(np.mean([normalize(dir1),normalize(dir2)],0))
   maxAngle = abs(angle(dir1,dir2))/2
   
   nodes = np.array(edgeMap.keys())
   angles = map(lambda x: abs(angle(vD.vertices[x,:]-vD.vertices[lineStart,:], dirvec, np.inf)), nodes)
   valid = filter(lambda x: not lineIntersectsPath([vD.vertices[x,:],vD.vertices[lineStart,:]], paths), nodes[angles < maxAngle])
   ps = []
   if len(valid) > 0:
      sValid = sorted(valid, None, lambda x: np.linalg.norm(vD.vertices[x,:]-vD.vertices[traceStart,:]), False)
      usedVerts = set([])
      while len(sValid) > 0:
         end = sValid.pop()
         if end not in usedVerts:
            usedVerts.add(end)
            if shortestPaths is None:
               p=shortestPath(G,end, traceStart)
            else:
               p = [traceStart]
               cur = traceStart
               while cur != end:
                  n = shortestPaths[end, cur]
                  p.append(n)
                  cur = n
            usedVerts.update(p)
            ps.append(p)
   else: ps.append([traceStart])

   return ps

"""
greedy version
"""
def traceArm2(edgeMap, traceStart, vD, dir1, dir2, paths, shortestPaths=None, G=None, lineStart=None):
   if lineStart is None:
      lineStart = traceStart
   dirvec = normalize(np.mean([normalize(dir1),normalize(dir2)],0))
   maxAngle = abs(angle(dir1,dir2))/2
   
   path = [start]
   oldlength = 0
   candidates = list(edgeMap[start].difference(path))
   while len(candidates)>0:
      angles = map(lambda x: abs(angle(vD.vertices[x,:]-vD.vertices[lineStart,:], dirvec, np.inf)), nodes)
      
      newCoordinate = vD.vertices[candidates[np.argmin(angles)]]
      newlength = np.linalg.norm(new-coordinate - vD.vertices[start])
      curEdge = Path([vD.vertices[start], new-coordinate], [1, 2])
      if np.any(map(lambda p: curEdge.intersects_path(p, False), paths)) or oldlength > newlength or min(angles) > maxAngle:
         break
      oldlength = newlength
      path.append(candidates[np.argmin(angles)])
      candidates = list(edgeMap[path[-1]].difference(path))
   return path


