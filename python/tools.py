import numpy as np
from itertools import *

def first(coll):
   return coll[0]

def partition(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def reductions(f, coll, init=None):
   l = list(coll if init is None else chain([init],coll))
   res = [l[0]]
   for item in l[1:]:
      res.append(f(res[-1], item))
   return res

def groupbyKeys(keys, coll):
   return imap(lambda x: map(lambda y: y[1], x[1]), \
               groupby(sorted(zip(keys, coll), None, first), first))

def dist(x,y):
   return np.linalg.norm(x-y)


def copy(listOfLists):
   return map(list, listOfLists)

"""
Returns list of best indices for items in setList, such that reduce(add, [setList[i] for i in idx]) returns list of unique items
"""
def resolveConflicts(setList, quality, minQuality=1): # quality between 0 and 1, 1 is best
   sIndQ = sorted(zip(range(len(setList)), quality), None, lambda x: x[1], True)
   idx = []
   usedEdges = set([])
   for i, q in sIndQ:
      if q < minQuality: break
      print q
      if len(usedEdges.intersection(setList[i])) == 0:
         idx.append(i)
         usedEdges = usedEdges.union(setList[i])
   return idx
