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
   #accLengths = reductions(add, imap(len, vertexLists))
   #invalidConn = map(lambda x: tuple(sort(x-1,x%len(points))), accLengths) 
   #ridgeValidity2 = map(lambda x: x not in invalidConn, voronoiDiagram.ridge_points)
   skeletonEdges = np.array(voronoiDiagram.ridge_vertices)[ridgeValidity]
   skeletonEdges = filter(lambda e: True != (lineIntersectsPath(voronoiDiagram.vertices[e,:], mplpaths)), skeletonEdges) # effekt?
   return voronoiDiagram, skeletonEdges, ridgeValidity


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
                  newEdges.extend(edges)
   return newEdges

"""
Returns list of best indices for items in setList, such that reduce(add, [setList[i] for i in idx]) returns list of unique items
"""
def resolveWedgeConflicts(head, arms, quality, minQuality=1): # quality between 0 and 1, 1 is best, TODO
   sIndQ = sorted(zip(range(len(head)), quality), None, lambda x: x[1], True)
   idx = []
   usedEdges = set([])
   
   for i, q in sIndQ:
      if q < minQuality: break
      
      if len(usedEdges.intersection(setList[i])) == 0:
         idx.append(i)
         usedEdges = usedEdges.union(setList[i])
      else:
         print ""
   return idx

def getContourWedges(cycle, vD, edgeMap, paths, shortestPaths, minSimilarity):
   wedges = []
   tri = fitTriangle(vD.vertices[cycle,:])
   similarity = 1-tri["rel-err"]
   if similarity >= minSimilarity:
      verts = cycle[tri["triangle"]]
      arms1 = traceArmComplete2(edgeMap, verts[0], vD, vD.vertices[verts[0]]-vD.vertices[verts[1]], vD.vertices[verts[0]]-vD.vertices[verts[2]], paths, shortestPaths)
      arms2 = traceArmComplete2(edgeMap, verts[1], vD, vD.vertices[verts[1]]-vD.vertices[verts[2]], vD.vertices[verts[1]]-vD.vertices[verts[0]], paths, shortestPaths)
      arms3 = traceArmComplete2(edgeMap, verts[2], vD, vD.vertices[verts[2]]-vD.vertices[verts[0]], vD.vertices[verts[2]]-vD.vertices[verts[1]], paths, shortestPaths)
      wedges = map(lambda arms: ContourWedge(cycle, list(arms), vD.vertices),product(arms1, arms2, arms3))
   return wedges

def detectContourWedges(edges, edges2, vD, paths, minTriangleSimilarity=0.7, maxHeadEdges=20, maxHeadLength=400):
   cycles, cycleHints = getCycles(edges, vD, maxHeadEdges, maxHeadLength)
   N = len(vD.vertices)
   dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
   graph = getUndirAdjMatrix(edges, N, dists)
   a, shortestPaths = shortest_path(graph,'auto', False, True)
   
   if len(cycles) > 0:
      edgeMap = getEdgeMap(edges)
      
      sim1 = map(lambda c: triangleSimilarity(vD.vertices[c,:], "area"), cycles)
      ws1 = np.array(map(lambda x: x[0], sim1))
      verts1 = np.array(map(lambda x,c: np.array(c)[x[1],], sim1, cycles))
      #ws2 = map(lambda c: triangleSimilarity(vD.vertices[c,:], "fourier", 5, True), cycles)
      #print "area:", ws1
      #print "fourier:", ws2
      #arms = map(lambda tri: map(lambda i: traceArmComplete(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths,shortestPaths), range(len(tri))), verts1)
      #idx = resolveWedgeConflicts(map(getCycleEdges, cycles), arms, ws1, minTriangleSimilarity)

      idx = resolveConflicts(map(getCycleEdges, cycles), ws1, minTriangleSimilarity)
      
      wedgeCycles = replace(idx, cycles)
      wedgeHeads = verts1[idx,:]
      
      
      #wedgeArms = map(lambda tri: map(lambda i: traceArmGreedy(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths, 3), range(len(tri))), wedgeHeads)
      wedgeArms = map(lambda tri: map(lambda i: traceArmComplete(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths,shortestPaths), range(len(tri))), wedgeHeads)
      #wedgeArms = map(lambda tri: map(lambda i: traceArmComplete2(edgeMap, tri[i],vD, vD.vertices[tri[i]]-vD.vertices[tri[(i+1)%3]], vD.vertices[tri[i]]-vD.vertices[tri[(i+2)%3]], paths,shortestPaths), range(len(tri))), wedgeHeads)
      
   else:
      wedgeHeads = []
      wedgeArms = []
      wedgeCycles = []
   return wedgeArms, wedgeHeads, cycles, cycleHints, wedgeCycles

def getSolidWedges(joint, vD, edgeMap, paths, shortestPaths, minWd):
   wedges = []
   wd = np.linalg.norm(vD.vertices[joint]-vD.points[filter(lambda v: joint in v[0], zip(vD.ridge_vertices, vD.ridge_points))[0][1][0]])

   if wd >= minWd:
      verts = list(edgeMap[joint])
      n = len(verts)
      armsListList = map(lambda i: traceArmComplete2(edgeMap, verts[i], vD, vD.vertices[verts[i]]-vD.vertices[verts[(i+1)%n]], vD.vertices[verts[i]]-vD.vertices[verts[(i+2)%n]], paths, shortestPaths), range(n))
      wedges = list(chain(*map(lambda armsList: map(lambda arms: SolidWedge(joint, list(arms), vD.vertices), \
                                                    product(*armsList)), \
                               combinations(armsListList,3))))
   return wedges


def detectSolidWedges(edges, exclude, usedEdges, vD, paths, minWd = 0, minFreeSides = 1):
   edgeMap = getEdgeMap(edges)
   N = len(vD.vertices)
   dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
   graph = getUndirAdjMatrix(edges, N, dists)
   a, shortestPaths = shortest_path(graph,'auto', False, True)
   
   joints = set(getNodes(edgeMap, lambda deg: deg>2)).difference(exclude)
   
   wedges = map(lambda j: sorted(getSolidWedges(j, vD, edgeMap, paths, shortestPaths, minWd), None, lambda x: x.angleProb(), True)[0], joints)

   chosenWedges = []
   for w in sorted(wedges, None, lambda x: x.angleProb(), True):
      print w.angleProb()
      print w.headEdges
      print w.armEdges
      if len(usedEdges.intersection(w.headEdges)) <= (3-minFreeSides):
         chosenWedges.append(w)
         usedEdges.update(w.headEdges)
         usedEdges.update(w.armEdges)

   return chosenWedges


"""
Allows directions between dir1 and dir2
Returns complete path
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
      if stepsSinceInvalid > nAllowedInvalid:
         break
      candidates = list(edgeMap[path[-1]].difference(path))
   return path[:-stepsSinceInvalid] if stepsSinceInvalid != 0 else path


"""
Does not return complete path (apply shortes path afterwards if needed)
"""
def traceArmComplete(edgeMap, start, vD, dir1, dir2, paths, shortestPaths):
   ends = [start]
   dirvec = np.mean([dir1,dir2],0)
   maxAngle = abs(angle(dir1,dir2))/2

   nodes = np.array(edgeMap.keys())
   angles = map(lambda x: abs(angle(vD.vertices[x,:]-vD.vertices[start,:], dirvec, np.inf)), nodes)
   
   valid = filter(lambda x: not lineIntersectsPath([vD.vertices[x,:],vD.vertices[start,:]], paths), nodes[angles < maxAngle])
   if len(valid) > 0:
      lengths = map(lambda x: np.linalg.norm(vD.vertices[x,:]-vD.vertices[start,:]), valid)
      ends.append(sorted(valid, None, lambda x: np.linalg.norm(vD.vertices[x,:]-vD.vertices[start,:]), True)[0])
   
   p = [start]
   cur = start
   while cur != ends[-1]:
      n = shortestPaths[ends[-1], cur]
      p.append(n)
      cur = n

   return p
   
"""
Does not return complete path (apply shortes path afterwards if needed)
"""
def traceArmComplete2(edgeMap, start, vD, dir1, dir2, paths, shortestPaths):
   dirvec = np.mean([dir1,dir2],0)
   maxAngle = abs(angle(dir1,dir2))/2
   
   nodes = np.array(edgeMap.keys())
   angles = map(lambda x: abs(angle(vD.vertices[x,:]-vD.vertices[start,:], dirvec, np.inf)), nodes)
   valid = filter(lambda x: not lineIntersectsPath([vD.vertices[x,:],vD.vertices[start,:]], paths), nodes[angles < maxAngle])

   ps = []
   if len(valid) > 0:
      sValid = sorted(valid, None, lambda x: np.linalg.norm(vD.vertices[x,:]-vD.vertices[start,:]), False)
      
      usedVerts = set([])
      while len(sValid) > 0:
         end = sValid.pop()
         if end not in usedVerts:
            usedVerts.add(end)
            p = [start]
            cur = start
            while cur != end:
               n = shortestPaths[end, cur]
               p.append(n)
               cur = n
            usedVerts.update(p)
            ps.append(p)
   else: ps.append([start])

   return ps

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
