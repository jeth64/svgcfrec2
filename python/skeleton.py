import numpy as np
from scipy.spatial import Voronoi
from matplotlib.path import Path
from itertools import *
from operator import add

# own modules
from tools import *
from graph import *
from geometry import *

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
   #accLengths = reductions(add, imap(len, vertexLists))
   #invalidConn = map(lambda x: tuple(sort(x-1,x%len(points))), accLengths) 
   #ridgeValidity2 = map(lambda x: x not in invalidConn, voronoiDiagram.ridge_points)
   skeletonEdges = np.array(voronoiDiagram.ridge_vertices)[ridgeValidity]
   skeletonEdges = filter(lambda e: True != (lineIntersectsPath(voronoiDiagram.vertices[e,:], mplpaths)), skeletonEdges) # effekt?
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

def condSplitLine(nodes, coordinates, mplpaths):
   if len(nodes) > 1:
      if lineIntersectsPath(coordinates[[nodes[0], nodes[-1]], :], mplpaths):
         if len(nodes) > 2:
            edges = condSplitLine(nodes[:len(nodes)/2], coordinates, mplpaths) + condSplitLine(nodes[len(nodes)/2:], coordinates, mplpaths)
         else: edges = [[]]
      else: edges = [[nodes[0], nodes[-1]]]
   else: edges = [[]]
   return edges 

"""
Simplify Voronoi skeleton with splitting after connecting points
"""
def simplifySkeleton2(skeletonEdges, coordinates, mplpaths): # test!
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
                  print edges
                  newEdges.extend(edges)
   return newEdges

def detectContourWedges(edges, vD, paths, minTriangleSimilarity=0.7, maxHeadEdges=10, maxHeadLength=np.inf):
   cycles, cycleHints = getCycles(edges, vD, maxHeadEdges, maxHeadLength)
   if len(cycles) > 0:
      sim1 = map(lambda c: triangleSimilarity(vD.vertices[c,:], "area"), cycles)
      ws1 = np.array(map(lambda x: x[0], sim1))
      verts1 = np.array(map(lambda x,c: np.array(c)[x[1],], sim1, cycles))
      #ws2 = map(lambda c: triangleSimilarity(vD.vertices[c,:], "fourier", 5, True), cycles)
      #print "area:", ws1
      #print "fourier:", ws2

      idx = resolveConflicts(map(getCycleEdges, cycles), ws1, minTriangleSimilarity)
      wedgeHeads = verts1[idx,:]
      edgeMap = getEdgeMap(edges)
      #wedgeArms = map(lambda tri: map(lambda i: traceArmGreedy(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths, 3), range(len(tri))), wedgeHeads)
      wedgeArms = map(lambda tri: map(lambda i: traceArmGreedy(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths), range(len(tri))), wedgeHeads)
      
      print wedgeArms
   else:
      wedgeHeads = []
      wedgeArms = []
   return wedgeArms, wedgeHeads, cycles, cycleHints

def detectSolidWedges(edges, vD, maxWd = np.inf):
   edgeMap = getEdgeMap(edges)
   joints = getNodes(edgeMap, lambda deg: deg>2)
   #edgeTriples = map(lambda i: map(lambda v: (i, v), edgeMap[i]), isecs))]

   # chooses random point as distances should be equal:
   wd = map(lambda j: np.linalg.norm(vD.vertices[j]-vD.points[filter(lambda v: j in v[0], zip(vD.ridge_vertices, vD.ridge_points))[0][1][0]]), joints)
   print wd
   
   return wd

"""
Allows directions between dir1 and dir2
"""
def traceArmGreedy(edgeMap, start, vD, dir1, dir2, paths, nAllowedInvalid): #tested: good
   dirvec = np.mean([dir1,dir2],0)
   maxAngle = abs(angle(dir1,dir2))/2
   allowed = 3 
   path = [start]
   oldlength = 0
   candidates = list(edgeMap[start].difference(path))
   stepsSinceInvalid = 0
   while len(candidates)>0:
      angles = map(lambda x: abs(angle(dirvec,x-vD.vertices[start])), vD.vertices[candidates])
      newCoordinate = vD.vertices[candidates[np.argmin(angles)]]
      newlength = np.linalg.norm(newCoordinate-vD.vertices[start])
      if oldlength > newlength or min(angles) > maxAngle or lineIntersectsPath([vD.vertices[start], newCoordinate], paths):
         stepsSinceInvalid = stepsSinceInvalid + 1
      else: stepsSinceInvalid = 0
      oldlength = newlength
      path.append(candidates[np.argmin(angles)])
      if stepsSinceIsec > nAllowedInvalid:
         break
      candidates = list(edgeMap[path[-1]].difference(path))
   return path[:-stepsSinceIsec] if stepsSinceIsec != 0 else path

def traceArmCompleteRec(edgeMap, pred, start, node, pathlist, vD):
   # side-effect on pathlist
   for v in edgeMap[node]:
      if v != pred[node]:
         if lineIntersectsPath([vD.vertices[start],vD.vertices[v]],paths):
            pathlist.append(traceBack(pred, node, connLimit))
         else:
            if pred[v] ==-1:
               pred[v] = node
               traceArmCompleteRec(edgeMap, pred, start, v, pathlist, vD)
               pred[v] = -1


def traceArmComplete(edgeMap, start, vD, center, paths):
   traceArmCompleteRec(edgeMap, pred, start, node, pathlist, vD)
   return pathlist
   

"""
Returns 3 limbs of wedge skeleton, center as coordinate
"""
def wedgeSkeleton(vertices, edgeMap, vD, paths, center, greedy = True):
   return map(lambda ind: traceArmGreedy(edgeMap, vD, ind, center, paths), vertices)


"""
Limbs consisting of vertex lists, modes :double or :quad :vert :vert-quad
limbs may not start with same vertex if modes :vert or :vert-quad are used
"""
def wedgeFromSkeleton(limbs, vD, center = None, mode ="double"):
   corners = map(lambda x:x[-1],limbs)
   startpts = map(lambda x:x[0], limbs)
   if center is None:
      center = np.mean(vD.vertices[startpts], 0)
   if mode == "double":
      matrix = np.array(map(lambda c1, c2: [vD.vertices[c1], center, center, vD.vertices[c2]], \
                            corners, np.roll(corners, -1)))
   elif mode == "quad":
      matrix = np.array(map(lambda c1, c2: elevateBezier([vD.vertices[c1], center, vD.vertices[c2]]), \
                            corners, np.roll(corners, -1)))
   elif mode == "vert":
      matrix = np.array(map(lambda c1, c2, s: [vD.vertices[c1], vD.vertices[s], vD.vertices[s], vD.vertices[c2]], \
                            corners, np.roll(corners, -1), np.roll(startpts, 1)))
   elif mode == "vert-quad":
      matrix = np.array(map(lambda c1, c2, s: elevateBezier([vD.vertices[c1], vD.vertices[s], vD.vertices[c2]]), \
                            corners, np.roll(corners, -1), np.roll(startpts, 1)))
   elif mode == "mid":
      matrix = np.array(map(lambda c1, c2, s: [vD.vertices[c1], np.mean([vD.vertices[c1], vD.vertices[s]], 0), np.mean([vD.vertices[c2], vD.vertices[s]], 0), vD.vertices[c2]], \
                            corners, np.roll(corners -1), np.roll(startpts, 1)))
   elif mode == "mid-cross":
      matrix = np.array(map(lambda c1, c2, s: [vD.vertices[c1], np.mean([vD.vertices[c2], vD.vertices[s]], 0), np.mean([vD.vertices[c1], vD.vertices[s]], 0), vD.vertices[c2]], \
                            corners, np.roll(corners, -1), np.roll(startpts, 1)))
   else:
      print "unknown reconstruction mode"
      matrix = []
   return matrix
