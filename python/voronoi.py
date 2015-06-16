import numpy as np
from scipy.spatial import Voronoi
from itertools import *


def skeleton(vertexLists, compoundShape): # allows intersections with path
   points = np.vstack(vertexLists)
   voronoiDiagram = Voronoi(points)
   vertexValidity = compoundShape.contains_points(voronoiDiagram.vertices)
   ridgeValidity1 = map(lambda x: vertexValidity[x[0]] and vertexValidity[x[1]], voronoiDiagram.ridge_vertices)
   #accLengths = reductions(add, imap(len, vertexLists))
   #invalidConn = map(lambda x: turple(sort(x-1,x%len(points))), accLengths) 
   #ridgeValidity2 = map(lambda x: x not in invalidConn, voronoiDiagram.ridge_points)
   skeletonEdges = voronoiDiagram[ridgeValidity]
   return voronoiDiagram, skeletonEdges


def simplifySkeleton(skeletonEdges, coordinates, compoundShape): #test
   newEdges = []
   if len(skeletonEdges) > 0:
      edgeMap = getEdgeMap(edges)
      while len(edgeMap) > 0: # works even if components not connected
         leaves = getNodes(edgemap, lambda deg: deg<2)
         isecs = getNodes(edgemap, lambda deg: deg>2)
         #visited = dict([]).fromkeys(edgeMap.iterkeys, False)
         if len(leaves) > 0: firstNode = leaves[0]
         elif len(isecs) > 0: firstNode = isecs[0]
         else: firstNode = edgeMap.keys[0]
         stack = [firstNode]
         while len(stack) > 0:
            vertex = stack.pop([0])
            neighbours = edgeMap.pop(vertex)
            for n in neighbours:
               node = vertex
               curNode = n
               print "coords", coordinates[[node, curNode], :]
               if curNode in edgeMap and not Path(coordinates[[node, curNode], :], [1, 2]).intersects_path(compoundShape):
                  prevNode = vertex
                  while curNode not in isecs.union(leaves):
                     cand = edgeMap.pop([curNode]).difference([prevNode])
                     line = Path(coordinates[[node, cand], :], [1, 2])
                     if line.intersects_path(compoundShape):
                        newEdges.append[(node,curNode)]
                        node = curNode
                     prevNode = curNode
                     curNode = cand
                  newEdges.append[(node,n)]
                     
                  
       
   return newEdges

