import numpy as np
from itertools import *
from operator import add
from tools import *

"""
Calculates area for simple polygons with shoelace formula
"""
def shoelace(vertices):
   if len(vertices) > 2:
      XY = np.transpose(vertices)
      down = np.dot(XY[0,:], np.roll(XY[1,:], -1, 0))
      up = np.dot(XY[0,:], np.roll(XY[1,:], 1, 0))
      # return abs(reduce(add, map(lambda xyi, xyinci: np.det(np.vstack([xyi xyinci])), vertices, )) / 2.0
      return abs(down-up)/2.0
   else: return 0

def fitTriangle(polygon):
   res = {} 
   if len(polygon) > 3:
      entireArea = shoelace(polygon)
      triVerts = list(combinations(polygon, 3))
      triEnum = list(combinations(range(len(polygon)), 3))
      triAreas = map(shoelace, triVerts)
      partsList = map(lambda idx: map(lambda s, e: list(chain(takewhile(lambda y: e != y, dropwhile(lambda x: s != x, cycle(range(len(polygon))))), [e])), \
                                                   idx, np.roll(idx, -1)), triEnum)
      absErrs = map(lambda parts: reduce(add, map(lambda part: shoelace(polygon[part,:]), parts)), partsList)
      minInd = np.argmin(absErrs)
      relErr = absErrs[minInd] / (entireArea + absErrs[minInd])
      bestInds = triEnum[minInd]
      res = {"triangle": bestInds, "original": polygon, "rel-err": relErr, "abs-err": absErrs[minInd]}
   elif len(polygon) < 3:
      res = {"triangle": None, "original": polygon, "rel-err": 1, "abs-err": np.inf}
   else:
      res = {"triangle": range(len(polygon)), "original": polygon, "rel-err": 0, "abs-err": 0}
   return res

"""
Calculate angle between two vectors
"""
def angle(u, v, forZeros = 0, oneDir = False):
   normalization = np.linalg.norm(u)*np.linalg.norm(v)
   if normalization == 0:
      return forZeros
   sig = np.sign(np.cross(u,v))
   if oneDir and sig < 0:
      alpha = 360 - np.arccos(np.dot(u,v) / normalization)
   else:
      alpha = sig * np.arccos(np.dot(u,v) / normalization)
   return alpha
   
   
"""
Returns Fourier descriptor of a polygon after Zahn Roskies 72
FD first discussed by cosgriff
"""
def polygonFourier(polygon, nDescriptors, polar=True):
   dirs = np.roll(polygon,-1, 0)-polygon
   deltaL = np.linalg.norm(dirs, 2,1)
   deltaPhi = map(angle, dirs, np.roll(dirs, -1,0))
   l = np.array(list(reductions(add, deltaL)))
   L = l[-1]
   a = np.empty((nDescriptors,))
   b = np.empty((nDescriptors,))
   A = np.empty((nDescriptors,))
   alpha = np.empty((nDescriptors,))
   for n in range(1,nDescriptors+1):
      a[n-1] = np.sum(np.multiply(deltaPhi,np.sin((2*np.pi*n/L)*l)))/(-n*np.pi)
      b[n-1] = np.sum(np.multiply(deltaPhi,np.cos((2*np.pi*n/L)*l)))/(n*np.pi)
      p = [a[n-1], b[n-1]]
      A[n-1] = np.linalg.norm(p) 
      alpha[n-1] = angle([1,0],p)
   # mu = - np.pi - np.sum(np.multiply(deltaPhi,l))/L
   return A,alpha

def triangleSimilarity(polygon, method, nDescriptors=5, polar=True):
   if method == "area":
      fT = fitTriangle(polygon)
      e = 1-fT["rel-err"]
      t = fT["triangle"]
   elif method == "fourier":
      tri = np.array([[0, 0],[2,0],[1,-1]])
      atri,btri = polygonFourier(tri, nDescriptors, polar)
      a,b = polygonFourier(polygon, nDescriptors, polar)
      ea = np.linalg.norm(atri-a)
      eb = np.linalg.norm(btri-b)
      print "ea", ea
      print "eb", eb
      e = 1-ea
      t = [0, np.floor(len(polygon)/3), floor(2*len(polygon)/3)] #stub
   else: print "Unknown method for triangle similarity"
   return e, t
