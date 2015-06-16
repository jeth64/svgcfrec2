import numpy as np

from sys import exit
from operator import add
from itertools import *

#own modules
from tools import *
from math import ceil

"""
Calculates cubic bezier curve for a line with evenly distributed control points
"""
def line2bezier(p1, p2):
   return [p1, (2*p1+p2)/3, (p1+2*p2)/3, p2]


"""
Elevates a quatratic bezier curve to cubic
"""
def elevateBezier(qbezier):
   return np.array([qbezier[0],(qbezier[0]+2*qbezier[1])/3, \
                   (qbezier[2]+2*qbezier[1])/3, qbezier[2]])

"""
Takes a list of point lists of the recursion steps and adds the list for the next recursion
"""
def deCasteljauRec(controls, t):
   coefs = np.array(controls[-1])
   if len (coefs) < 2:
      return controls
   else:
      return controls + deCasteljauRec([map(lambda Pi, Pinci: (1-t) * Pi + t * Pinci, coefs[:-1], coefs[1:])], t)

"""
Splits Bezier curve at position t into two curves
"""
def deCasteljauSplit(controls, t):
   cPointList = deCasteljauRec([list(controls)], t)
   return [np.array(list(map(first, cPointList))),
           np.array(list(reversed(map(lambda x: x[-1], cPointList))))]

# testvals: seg = np.array([[0,0],[1,1],[2,1],[3,0]])

"""
Calculates the point on the curve for position parameter t
"""
def deCasteljau(controls, t):
   cPointList = deCasteljauRec([list(controls)], t)
   return np.array(cPointList[-1][0])

"""
Generator that splits a cubic Bezier curve into n curves
"""
def deCasteljauSplitInto(n, segment):
   i = n      
   while i > 0:
      segments = deCasteljauSplit(segment, 1.0/i)
      yield segments[0]
      segment = segments[1]
      i -= 1



