import numpy as np
from scipy.spatial import Voronoi
from matplotlib.path import Path
from itertools import *
from operator import add
from scipy.sparse.csgraph import shortest_path


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
   #if len(voronoiDiagram.vertices) > 23279:
      #print map(lambda p: p.contains_point(voronoiDiagram.vertices[23279]), mplpaths)      
      #print "ps", pointsInShape(voronoiDiagram.vertices[range(23273, 23280),:], mplpaths)
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
   
"""
def getContourWedges(cycle, vD, edgeMap, paths, shortestPaths, minSimilarity=0, minProb=0.5):
   wedges = []
   tri = fitTriangle(vD.vertices[cycle,:])
   similarity = 1-tri["rel-err"]
   if similarity >= minSimilarity:
      verts = np.array(cycle)[tri["triangle"],]
      arms1 = traceArm(edgeMap, verts[0], vD, vD.vertices[verts[0]]-vD.vertices[verts[1]], vD.vertices[verts[0]]-vD.vertices[verts[2]], paths, shortestPaths)
      arms2 = traceArm(edgeMap, verts[1], vD, vD.vertices[verts[1]]-vD.vertices[verts[2]], vD.vertices[verts[1]]-vD.vertices[verts[0]], paths, shortestPaths)
      arms3 = traceArm(edgeMap, verts[2], vD, vD.vertices[verts[2]]-vD.vertices[verts[0]], vD.vertices[verts[2]]-vD.vertices[verts[1]], paths, shortestPaths)
      wedges = map(lambda arms: ContourWedge(cycle, list(arms), vD.vertices, similarity),product(arms1, arms2, arms3))
   return filter(lambda w: w.angleProb() >= minProb, wedges)

def detectContourWedges(edges, vD, paths, minTriangleSimilarity=0.7, maxHeadEdges=20, maxHeadLength=400):
   cycles, cycleHints = getCycles(edges, vD, maxHeadEdges, maxHeadLength)
   N = len(vD.vertices)
   dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
   graph = getUndirAdjMatrix(edges, N, dists)
   a, shortestPaths = shortest_path(graph,'auto', False, True)
   
   if len(cycles) > 0:
      edgeMap = getEdgeMap(edges)

      allWedges = imap(lambda c: sorted(getContourWedges(c, vD, edgeMap, paths, shortestPaths, minTriangleSimilarity), None, lambda x: x.angleProb(), True), cycles)
      wedges = map(first, ifilter(notEmpty, allWedges))  # choose best of all wedges for a circle

      usedHeadEdges = set([])
      chosenWedges = []
      for w in sorted(wedges, None, lambda x: x.angleProb()*x.probability, True):
         if len(usedHeadEdges.intersection(w.headEdges)) == 0:
            chosenWedges.append(w)
            usedHeadEdges.update(w.headEdges)
      
   else:
      chosenWedges = []
   return chosenWedges, cycles, cycleHints

def getSolidWedges(joint, vD, edgeMap, paths, shortestPaths, minWd = 0, minProb = 0.5):
   wedges = []
   wd = np.linalg.norm(vD.vertices[joint]-vD.points[filter(lambda v: joint in v[0], zip(vD.ridge_vertices, vD.ridge_points))[0][1][0]])
   #print wd
   if wd >= minWd:
      #print "                       ",joint
      verts = list(edgeMap[joint])
      #print verts
      n = len(verts)
      # TODO: whaat of the following is best
      armsListList = map(lambda i: traceArm(edgeMap, verts[i], vD, vD.vertices[verts[i]]-vD.vertices[verts[(i+1)%n]], vD.vertices[verts[i]]-vD.vertices[verts[(i+2)%n]], paths, shortestPaths,joint), range(n))
      #armsListList = map(lambda i: traceArm(edgeMap, verts[i], vD, vD.vertices[joint]-vD.vertices[verts[(i+1)%n]], vD.vertices[joint]-vD.vertices[verts[(i+2)%n]], paths, shortestPaths,joint), range(n))
      wedges = list(chain(*map(lambda armsList: map(lambda arms: SolidWedge(joint, list(arms), vD.vertices, wd), \
                                                    product(*armsList)), \
                               combinations(armsListList,3))))
      
   return filter(lambda w: w.angleProb() >= minProb, wedges)

def detectSolidWedges(edges, exclude, usedEdges, vD, paths, minWd = 0, minFreeSides = 1):
   edgeMap = getEdgeMap(edges)
   N = len(vD.vertices)
   dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
   graph = getUndirAdjMatrix(edges, N, dists)
   a, shortestPaths = shortest_path(graph,'auto', False, True)
   
   joints = set(getNodes(edgeMap, lambda deg: deg>2)).difference(exclude)
   
   allWedges = imap(lambda j: sorted(getSolidWedges(j, vD, edgeMap, paths, shortestPaths, minWd), None, lambda x: x.angleProb(), True), joints)
   wedges = map(first,ifilter(notEmpty, allWedges))  # choose best of all wedges for a center

   chosenWedges = []
   for w in sorted(wedges, None, lambda x: x.angleProb(), True):
      if len(usedEdges.intersection(w.headEdges)) <= (3-minFreeSides):
         #print len(usedEdges.intersection(w.headEdges))
         #print w.headEdges
         #print w.probability
         #print w.angleProb()
         chosenWedges.append(w)
         usedEdges.update(w.headEdges)
         usedEdges.update(w.armEdges)

   return chosenWedges
"""

"""
Returns multiple possible arms
"""
def traceArm(edgeMap, traceStart, vD, dir1, dir2, paths, shortestPaths, lineStart=None):
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
