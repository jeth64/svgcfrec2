import numpy as np
from scipy.spatial import Voronoi
from matplotlib.path import Path
from itertools import *

# own modules
from tools import *
from graph import *

"""
Even-odd rule
"""
def pointsInShape(points, mplpaths):
   # compoundShape considers points in holes as in shape, thats why the insideness has to be checked differently
   return map(lambda point: np.sum(map(lambda p: p.contains_point(point), mplpaths))%2, points)
   
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
   #accLengths = reductions(add, imap(len, vertexLists))
   #invalidConn = map(lambda x: tuple(sort(x-1,x%len(points))), accLengths) 
   #ridgeValidity2 = map(lambda x: x not in invalidConn, voronoiDiagram.ridge_points)
   skeletonEdges = np.array(voronoiDiagram.ridge_vertices)[ridgeValidity]
   return voronoiDiagram, skeletonEdges


"""
Simplify Voronoi skeleton
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
                  while curNode not in endpoints:
                     cand = edgeMap.pop(curNode).difference([prevNode]).pop()
                     if lineIntersectsPath(coordinates[[node, cand], :], mplpaths):
                        newEdges.append([node,curNode])
                        node = curNode
                     prevNode = curNode
                     curNode = cand
                  newEdges.append([node,curNode])
   return newEdges

