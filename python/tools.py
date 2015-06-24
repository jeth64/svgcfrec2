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

   
