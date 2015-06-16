import numpy as np
from time                   import time

# own modules
from bezier import *


"""
Calculates sum of distances of consecutive control points
"""
def polypath(segment):
   return reduce(add, imap(dist, segment[:-1], segment[1:]))


"""
Splits a cubic Bezier curve until cut-off attribute evaluates smaller than given threshold
If split is dynamic, one recursion divides a curve into a dynamic number of curves,
if split is static, one recursion divides a curve into exactly 2 curves
"""
def splitSegmentTo(cutOffAttribute, threshold, segment,split):
   if cutOffAttribute == "distance":
      s = dist(segment[0], segment[-1])
   elif cutOffAttribute == "polypath":
      s = polypath(segment)
   else:
      print "Unknown discretization attributes"
      exit(1)
   if s < threshold:
      newSegments = [segment]
   else:
      if split == "dynamic": # TODO: whats the perfect split?
         segments = list(deCasteljauSplitInto(ceil(s/threshold), segment))
      elif split == "static":
         segments = deCasteljauSplit(segment, 0.5)
      newSegments = reduce(add, imap(lambda s: splitSegmentTo(cutOffAttribute, threshold, s,split), segments))
   return newSegments


"""
Splits the curves of a cubic Bezier spline until cut-off attribute evaluates smaller than given threshold
If split is dynamic, one recursion divides a curve into a dynamic number of curves,
if split is static, one recursion divides a curve into exactly 2 curves
"""
def splitSegmentsTo(cutOffAttribute, threshold, path, split = "static"):
  return np.array(reduce(add, imap(lambda x: splitSegmentTo(cutOffAttribute, threshold, x, split), path)))


def discretize(matrix, method, threshold):
   #test_discretization_methods()
   #exit(0)
   return splitSegmentsTo(method, threshold, matrix)[:,0,:]


def test_discretization_methods():
   method = "polypath"
   threshold = 0.1
   for i in range(3):
      mat =  np.random.rand(4,2)
      
      t1 = time()
      n1 = len(map(polypath,splitSegmentTo(method, threshold, mat, "dynamic")))
      end1 = time()-t1
      print "dynamic: "
      print " number of segments:", n1
      print " time:", end1, "\n"
      
      t2 = time()
      n2 = len(map(polypath,splitSegmentTo(method, threshold, mat, "static")))
      end2 = time()-t2
      print "static: "
      print " number of segments:", n2
      print " time:", end2, "\n"
      # produces less segments in less time
      
   return True

