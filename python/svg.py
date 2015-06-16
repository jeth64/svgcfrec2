import numpy as np

from xml.etree.ElementTree  import parse, SubElement
from re                     import finditer, split, match
from itertools import *
from matplotlib.path import Path
from matplotlib.transforms import Transform, Affine2D, composite_transform_factory

# own modules
from tools import *
from bezier import *


"""
Return element tree, root and namespace of svg document
"""
def getSVG(filename):
   tree = parse(filename)
   root = tree.getroot()
   if root.tag[-3:] != "svg":
      print "Tag of SVG file root must be 'svg'. Give valid SVG file."
      exit(1)
   namespace = root.tag[:-3]
   return tree, root, namespace


"""
Delete an element from tree
"""
def deleteElement(root, element):
   root.remove(element)

"""
Return list of n lists of mi kijx4x2 matrices representing n shapes (= path objects),
where mi the number of contours (= splines) a shape i consists of
and kij the number of curves a spline ij consists of
"""
def getPaths(tree, node, namespace, isUnclean = False): # TODO: implement unclean
   pathGroups = []
   parentMap = {c:p for p in tree.iter() for c in p}
   tag = (namespace + "path") if (namespace is not None) else "path"
   for path in node.iter(tag):
      transform = getTransform(parentMap, path)
      paths = pathstr2mplpaths(path.get("d"))
      newPaths = map(lambda p: transform.transform_path(p), paths)
      mats = map(lambda p: mplpath2matrix(p), newPaths)
      pathGroups.append(zip(mats, newPaths))
   return pathGroups

def getTransform(parentMap, node): # test, speed up by caching transform for node
   n = node
   transform = Affine2D()
   while True:
      tString = n.get("transform")
      if tString is not None:
         for m in reversed(list(finditer(r"(?P<method>\w+)\(\s*(?P<args>[^)]*)\)", tString))):
            args = m.group('args').replace(",", " ").split()
            method = m.group('method')
            if method == "matrix":
               transform = composite_transform_factory(transform, Affine2D.from_values(*args))
            elif method == "translate":
               transform.translate(args[0], 0 if len(args)<2 else args[1])
            elif method == "scale":
               transform.scale(*args)
            elif method == "rotate":
               if len(args) == 1:
                  transform.rotate_deg(args[0])
               else:
                  transform.rotate_deg_around(*args)
            elif method == "skewX":
               transform.skew_deg(args[0], 0)
            elif method == "skewY":
               transform.skew_deg(0, args[0])
      if n in parentMap:
         n = parentMap[n]
      else: break
   return transform
         
      
   

"""
Devides path string at 'moveto' commands (M/m) and returns a matrix for each substring
"""
def pathstr2mplpaths(pathstr):
   mats = []
   paths = []
   offset = np.array([0,0])
   for match in finditer("([Mm])[^Mm]+", pathstr):
      shift = np.array([0,0]) if match.group(1).isupper() else offset
      mplpath = splinestr2mplpath(match.group(0), shift)
      if len(mplpath) > 0:
         paths.append(mplpath)
         offset = paths[-1].vertices[-1,:]
   return paths
   


"""
For each command part calculates 4 points representing it
"""
def mplpath2matrix(path):
   matrix = []
   last = None
   for coords, code in path.iter_segments():
      vertices = coords.reshape(len(coords)/2,2)
      if code == 1:
         if last is not None:
            matrix.append(line2bezier(last, vertices[0]))
      elif code == 2:
         matrix.append(line2bezier(last, vertices[0]))
      elif code == 3:
         parts = partition(vertices, 2)
         matrix.append(elevateBezier(np.vstack([[last], vertices])))
      elif code == 4:
         matrix.append(np.vstack([[last], vertices]))
      elif code == 79:
         if np.all(last != matrix[0][0]):
            matrix.append(line2bezier(last, matrix[0][0]))
      else:
         print "Unknown Matplotlib command"
      last = vertices[-1]
   return np.array(matrix)


"""
For a string representing a contour returns matplotlib Path
"""
def splinestr2mplpath(pathstr, initialOffset):
   offset = initialOffset
   vertices = []
   codes = []
   for it in finditer('([A-DF-Za-df-z])([^A-DF-Za-df-z]+)',pathstr):
      cmd = it.group(1)
      nrstr = it.group(2).replace("-"," -").replace("e -", "e-").strip()
      numbers = imap(float, list(ifilter(lambda x: len(x) > 0, split(" |,", nrstr))))
      code, newVertices = translateToMplCmd(cmd, numbers, vertices)
      if len(newVertices) > 0:
         vertices.extend(newVertices)
         codes.extend([code]*len(newVertices))
         offset = vertices[-1][-1]
   return Path(vertices, codes)


"""
Translations from path commands to matplotlib codes
"""
mplCode = { "M": 1, "L": 2, "Q": 3, "C": 4, "Z": 79}


"""
Translate command-points pair to one suitable for matplotlib path
'predecing' is expected to be list of preceding spline point list
"""
def translateToMplCmd(cmd, numbers, preceding):
   offset = preceding[-1] if len(preceding) > 0 else np.array([0, 0])
   if cmd.upper() == "M":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["M"]
      vertices = list(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
   elif cmd.upper() == "L":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["L"]
      vertices = list(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
   elif cmd.upper() == "H":
      xs = numbers if cmd.isupper() else reductions(add, numbers, offset[0])[1:]
      code = mplCode["L"]
      vertices = zipwith(np.array, xs, repeat(offset[1]))
   elif cmd.upper() == "V":
      ys = numbers if cmd.isupper() else reductions(add, numbers, offset[1])[1:]
      code = mplCode["L"]
      vertices = zipwith(np.array, repeat(offset[0]), ys)
   elif cmd.upper() == "Q":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["Q"]
      vertices = list(points) if cmd.isupper() else makeAbsolute(points, offset, 2)
   elif cmd.upper() == "T":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["Q"]
      vertices = list(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
      vertices = addFirstControlforSmooth(vertices, preceding, 1)
   elif cmd.upper() == "C":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["C"]
      vertices =  list(points) if cmd.isupper() else makeAbsolute(points, offset, 3)
   elif cmd.upper() == "S":
      points = imap(np.array, partition(numbers,2))
      code = mplCode["C"]
      vertices = list(points) if cmd.isupper() else makeAbsolute(points, offset, 2)
      vertices = addFirstControlforSmooth(vertices, preceding, 2)
   elif cmd.upper() == "Z":
      code = mplCode["Z"]
      vertices = [preceding[0]]
   else:
      print "Path command unknown:" + cmd
      code = mplCode["L"]
      vertices = [] # offset
      print "Added point instead"
   return code, vertices

"""
Determine the first control point for smooth-curve-to commands
"""
def addFirstControlforSmooth(points, preceding, step):
   parts = list(chain(preceding,partition(points, step)))
   firstPts = imap(lambda i: 2*parts[i][-1]-parts[i][-2], xrange(len(parts)-1))
   return imap(lambda p, part: [p] + part[1:], firstPts, parts)

"""
Translate a list of stepwise relative points to a list of only absolute points
"""
def makeAbsolute(points, offset, step):
   parts = list(partition(points, step))
   firstPts = reductions(lambda x,y: x+y[-1], parts[:-1], np.array(offset))
   return np.vstack(map(lambda p, part: p + part, firstPts, parts))


"""
Determine path string for cubic Bezier spline
"""
def matrix2path(path):
   p = np.around(path,decimals=3)
   anker = p[0,0,:]
   pathString = "M{0[0]},{0[1]}".format(anker)
   for i in range(len(p)):
      for j in range(1,4):
         pathString = pathString + (" c" if (j == 1) else " ")  \
                  + "{0[0]},{0[1]}".format(p[i,j,:]-p[i,0,:])
   return pathString

"""
Add group to element tree 
"""
def addGroup(namespace, parent, attributes):
   tag = (namespace+"g") if (namespace is not None) else "g"
   return SubElement(parent, tag, attributes)

"""
Add path to element tree 
"""
def addPath(namespace, parent, path, attributes):
   tag = (namespace+"path") if (namespace is not None) else "path"
   attributes["d"] = matrix2path(path)
   return SubElement(parent, tag, attributes)

"""
Add line to element tree 
"""
def addLine(namespace, parent, line, attributes):
   tag = (namespace +"line") if (namespace is not None) else "line"
   attributes.update(izip(["x1", "y1", "x2", "y2"], imap(str, line.flatten)))
   return SubElement(parent, tag, attributes)

"""
Add circle to element tree 
"""
def addCircle(namespace, parent, center, attributes):
   tag = (namespace+"circle") if namespace is not None else "circle"
   attributes.update(izip(["cx", "cy"], imap(str, center)))
   return SubElement(parent, tag, attributes)

