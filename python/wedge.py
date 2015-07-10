import numpy as np
from itertools import *
from operator import add
from geometry import angle

class Wedge:
   
   def __init__(self, center, headEdges, arms, vertexMapping):
      self.headEdges = set(headEdges)
      self.armEdges = set(chain(*map(lambda arm: map(lambda v1, v2: tuple(sorted([v1, v2])),arm[:-1],arm[1:]), arms)))      
      self.deepestPoint = Wedge.__deepestPoint(self, arms, vertexMapping, center)
      # order vertices clockwise
      verts = map(lambda arm: vertexMapping[arm[-1],:], arms)
      if angle(verts[0]-self.deepestPoint,verts[1]-self.deepestPoint,0,True) < angle(verts[0]-self.deepestPoint,verts[2]-self.deepestPoint, 0, True):
         self.vertices = np.array(verts)
      else:
         self.vertices = np.array(list(reversed(verts)))
      
      print arms
      print "deepest point", self.deepestPoint
      print "v",verts,angle(verts[0]-self.deepestPoint,verts[2]-self.deepestPoint, 0, True)

      print "vertices", self.vertices
      print "dirs", self.directions()
      print "lengths", self.lengths()
      print "angles", self.angles()
      print "angleProb", self.angleProb()
      print reduce(add, self.angles())
   
   def __intersection(self, l1, l2):
      A = np.vstack((l1[1]-l1[0], l2[1]-l2[0])).T
      x = np.linalg.solve(A,(l2[1]-l1[0]).T)
      return l1[0] + x[0] * (l1[1]-l1[0])

   def __deepestPoint(self, arms, vertexMapping, center):
      wedgeArms = np.array(map(lambda a: [vertexMapping[a[0]],vertexMapping[a[-1]]], arms))
      nonDegen = map(lambda arm: len(arm) > 1, arms)
      nNonDegen = reduce(add, nonDegen) 
      if nNonDegen == 3: # mean of all 3 intersections
         deepest = np.mean(map(lambda l: self.__intersection(*l), [wedgeArms[[0,1]],wedgeArms[[1,2]], wedgeArms[[2,0]]]), 0)
      elif nNonDegen == 2: # intersection of non-degenerate arms
         deepest = self.__intersection(*np.compress(nonDegen,wedgeArms,0))
      elif nNonDegen == 1:
         lines = map(lambda nD, a: a if nD else [center, a[-1]], nonDegen, wedgeArms)
         deepest = np.mean(map(lambda l1,l2: self.__intersection(l1,l2), lines, np.roll(lines,1,0)),0)
      else:
         deepest = center # also possible: barycenter
      return deepest

   def directions(self):
      return map(lambda v: v-self.deepestPoint, self.vertices)
      
   def lengths(self):
      return map(lambda v: np.linalg.norm(v-self.deepestPoint), self.vertices)

   def angles(self):
      return np.degrees(map(lambda v1,v2: angle(v1-self.deepestPoint, v2-self.deepestPoint,0,True), self.vertices, np.roll(self.vertices,-1,0))) # -1?

   def angleProb(self):
      err = reduce(add, map(lambda a: abs(120-a), self.angles()))
      return 1- (err / (360+err)) # 0<= x<=1, aber fehlerwinkel sollten reichen

   def reconstruct(self, mode = "double"):
      if mode == "double":
         matrix = np.array(map(lambda v1, v2: [v1, self.deepestPoint, self.deepestPoint, v2], \
                               self.vertices, np.roll(self.vertices, -1,0)))
      elif mode == "quad":
         matrix = np.array(map(lambda v1, v2: elevateBezier([v1, self.deepestPoint, v2]), \
                               self.vertices, np.roll(self.vertices, -1,0)))
      else:
         print "unknown reconstruction mode"
         matrix = []
      return matrix


class ContourWedge(Wedge):
   def __init__(self, cycle, arms, vertexMapping):
      center = np.mean(map(lambda a: vertexMapping[a[0]], arms),0)
      headEdges = map(lambda v1, v2: tuple(sorted([v1, v2])),cycle[:-1],np.roll(cycle,1))      
      Wedge.__init__(self, center, headEdges, arms, vertexMapping)
      

      
class SolidWedge(Wedge):
   def __init__(self, joint, arms, vertexMapping):
      center = vertexMapping[joint]
      headEdges = map(lambda arm: tuple(sorted([joint,arm[0]])),arms)
      Wedge.__init__(self, center, headEdges, arms, vertexMapping)

      
      
   
   
