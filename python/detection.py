import numpy as np
from itertools import *
from scipy.sparse.csgraph import shortest_path
from sys import exit

from geometry import *
from skeleton import *
from graph import *


def getContourWedgesFor(cycle, vD, edgeMap, paths, shortestPaths, minSimilarity, minProb):
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

def getContourWedges(vD, edges, edgeMap, paths, shortestPaths, minSimilarity, minProb, maxHeadEdges, maxHeadLength):
   cycles, cycleHints = getCycles(edges, vD, maxHeadEdges, maxHeadLength)
   if len(cycles) > 0:
      cWedges = filter(notEmpty, imap(lambda c: sorted(getContourWedgesFor(c, vD, edgeMap, paths, shortestPaths, minSimilarity, minProb), \
                                                       None, lambda x: x.angleProb(), True), cycles))
   else: cWedges = []
   return cWedges, cycles, cycleHints

   
def getSolidWedgesFor(joint, vD, edgeMap, paths, shortestPaths, minWS, minProb, distRangeS):
   wedges = []
   dist = np.linalg.norm(vD.vertices[joint]-vD.points[filter(lambda v: joint in v[0], zip(vD.ridge_vertices, vD.ridge_points))[0][1][0]])
   if dist < distRangeS[0]: WS = 0
   elif dist > distRangeS[1]: WS = 1
   else: WS = (dist-distRangeS[0])/(distRangeS[1]-distRangeS[0])

   print "range", distRangeS
   print dist
   print WS
   
   if WS >= minWS:
      verts = list(edgeMap[joint])
      n = len(verts)
      # TODO: whaat of the following is best
      armsListList = map(lambda i: traceArm(edgeMap, verts[i], vD, vD.vertices[verts[i]]-vD.vertices[verts[(i+1)%n]], vD.vertices[verts[i]]-vD.vertices[verts[(i+2)%n]], paths, shortestPaths,joint), range(n))
      #armsListList = map(lambda i: traceArm(edgeMap, verts[i], vD, vD.vertices[joint]-vD.vertices[verts[(i+1)%n]], vD.vertices[joint]-vD.vertices[verts[(i+2)%n]], paths, shortestPaths,joint), range(n))
      wedges = list(chain(*map(lambda armsList: map(lambda arms: SolidWedge(joint, list(arms), vD.vertices, WS), \
                                                    product(*armsList)), \
                               combinations(armsListList,3))))
      
   return filter(lambda w: w.angleProb() >= minProb, wedges)

def getSolidWedges(vD, edgeMap, paths, shortestPaths, minWd, minProb, distRangeS, exclude = []):
   joints = set(getNodes(edgeMap, lambda deg: deg>2)).difference(exclude)
   sWedges = filter(notEmpty, imap(lambda j: sorted(getSolidWedgesFor(j, vD, edgeMap, paths, shortestPaths, minWd, minProb, distRangeS), None, lambda x: x.angleProb(), True), joints))
   return sWedges
   

"""
Find wedges and choose a subset

parameters:
- minAngleProb
- noContours
- noSolids

 for wedge contours:
  - minWC
  - maxHeadLength
  - maxHeadEdges
(last two only to avoid outliers)

 for solid wedges:
  - minWS
  - minFreeSides
  
"""
def getWedges(edges, vD, paths,  strategy, parameters):
   contourWedges = []
   solidWedges = []
   cycles = []
   cycleHints = []

   cWedges = []
   sWedges = []

   edgeMap = getEdgeMap(edges)
   N = len(vD.vertices)
   dists = map(lambda e: np.linalg.norm(vD.vertices[e[0],:]-vD.vertices[e[1],:]), edges)
   graph = getUndirAdjMatrix(edges, N, dists)
   a, shortestPaths = shortest_path(graph,'auto', False, True)

   # (Get all wedges with a minimum triangle similarity or center distance to contour for contour and solid wedges respectively
   # in each case result is a list of possible wedges for each center or triangle; list can be empty) old text

   # apply given strategy 
   if strategy == "threshold":
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         contourWedges = list(chain(*cWedges))
      if not parameters["noSolids"]:
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"])
         solidWedges = list(chain(*sWedges))
         
   elif strategy == "threshold-repulse": # one wedge per center or triangle
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         contourWedges = map(first, cWedges)
      if parameters["noSolids"]:
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"])
         solidWedges = map(first, sWedges)
   
   elif strategy == "contour-fill": # wedges sorted by product of angle probs and other prob, contour wedges preferred, greedy
      usedHeadEdges = set([])
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         for w in sorted(list(chain(*cWedges)), None, lambda x: x.angleProb()*x.weight, True):
            if len(usedHeadEdges.intersection(w.headEdges)) == 0: # two heads may not use the same edge
               contourWedges.append(w)
               usedHeadEdges.update(w.headEdges)

      if not parameters["noSolids"]:
         usedArmEdges = reduce(lambda x, y: x.union(y), imap(lambda w: w.armEdges, contourWedges), set([]))
         usedEdges = usedArmEdges.union(usedHeadEdges)
         usedVerts = list(chain(*usedHeadEdges))      
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"], usedVerts)

         if len(sWedges) > 0:
            for w in sorted(list(chain(*sWedges)), None, lambda x: x.angleProb()*x.weight, True):
               if len(usedEdges.intersection(w.headEdges)) <= (3-parameters["minFreeSides"]): # at least a certain amount of the edges of the head must be free
                  solidWedges.append(w)
                  usedEdges.update(w.headEdges)
                  usedEdges.update(w.armEdges)
   
      
   elif strategy == "contour-fill-angles": # wedges sorted by angle probs
      usedHeadEdges = set([])
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         for w in sorted(list(chain(*cWedges)), None, lambda x: x.angleProb(), True):
            if len(usedHeadEdges.intersection(w.headEdges)) == 0: # two heads may not use the same edge
               contourWedges.append(w)
               usedHeadEdges.update(w.headEdges)

      if not parameters["noSolids"]:
         usedArmEdges = reduce(lambda x, y: x.union(y), imap(lambda w: w.armEdges, contourWedges), set([]))
         usedEdges = usedArmEdges.union(usedHeadEdges)
         usedVerts = list(chain(*usedHeadEdges))      
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"], usedVerts)

         if len(sWedges) > 0:
            for w in sorted(list(chain(*sWedges)), None, lambda x: x.angleProb(), True):
               if len(usedEdges.intersection(w.headEdges)) <= (3-parameters["minFreeSides"]): # at least a certain amount of the edges of the head must be free
                  solidWedges.append(w)
                  usedEdges.update(w.headEdges)
                  usedEdges.update(w.armEdges)

   elif strategy == "contour-fill-weight": # wedges sorted by weights
      usedHeadEdges = set([])
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         for w in sorted(list(chain(*cWedges)), None, lambda x: x.weight, True):
            if len(usedHeadEdges.intersection(w.headEdges)) == 0: # two heads may not use the same edge
               contourWedges.append(w)
               usedHeadEdges.update(w.headEdges)

      if not parameters["noSolids"]:
         usedArmEdges = reduce(lambda x, y: x.union(y), imap(lambda w: w.armEdges, contourWedges), set([]))
         usedEdges = usedArmEdges.union(usedHeadEdges)
         usedVerts = list(chain(*usedHeadEdges))      
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"], usedVerts)

         if len(sWedges) > 0:
            for w in sorted(list(chain(*sWedges)), None, lambda x: x.weight, True):
               if len(usedEdges.intersection(w.headEdges)) <= (3-parameters["minFreeSides"]): # at least a certain amount of the edges of the head must be free
                  solidWedges.append(w)
                  usedEdges.update(w.headEdges)
                  usedEdges.update(w.armEdges)

   elif strategy == "contour-fill-threshold": # like threshold-repulse but exclude vertices for center that are part of a wedge head
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])
         contourWedges = map(first, cWedges)
      if parameters["noSolids"]:
         usedHeadEdges = map(lambda w: w.headEdges, contourWedges)
         usedVerts = list(chain(*usedHeadEdges))  
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"], usedVerts)
         solidWedges = map(first, sWedges)

   elif strategy == "balanced-angles": # wedges sorted by angle probs, types of wedges equally weighted (other than angles senseless since not comparable)
      usedHeadEdges = set([])
      if not parameters["noContours"]:
         cWedges, cycles, cycleHints = getContourWedges(vD, edges, edgeMap, paths, shortestPaths, parameters["minWC"], parameters["minAngleProb"], parameters["maxHeadEdges"], parameters["maxHeadLength"])

      if not parameters["noSolids"]:
         sWedges = getSolidWedges(vD, edgeMap, paths, shortestPaths, parameters["minWS"], parameters["minAngleProb"], parameters["distRangeS"], usedVerts)
   
      for w in sorted(list(chain(chain(*cWedges),chain(*sWedges))), None, lambda x: x.angleProb(), True):
         if len(usedHeadEdges.intersection(w.headEdges)) == 0: # two heads may not use the same edge
            if isinstance(w, ContourWedge):
               contourWedges.append(w)
            else:
               solidWedges.append(w)

            usedHeadEdges.update(w.headEdges)

   else:
      print "Reduction strategy unknown!"
      print "Ending program..."
      exit
      
   return contourWedges, solidWedges, cycles, cycleHints


