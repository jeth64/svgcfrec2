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
Delete paths from tree
"""
def deletePaths(tree, namespace, layer):
   parentMap = getParentMap(tree)
   tag = (namespace + "path") if (namespace is not None) else "path"
   for element in layer.iter(tag):
      parent = parentMap[element]
      parent.remove(element)

"""
Get parent map for specified tree
"""
def getParentMap(tree):
   return {c:p for p in tree.iter() for c in p}

"""
Return list of n lists of mi kijx4x2 matrices representing n shapes (= path objects),
where mi the number of contours (= splines) a shape i consists of
and kij the number of curves a spline ij consists of
"""
def getPaths(tree, node, namespace, isUnclean = False): # TODO: implement unclean
   pathGroups = []
   parentMap = getParentMap(tree)
   tag = (namespace + "path") if (namespace is not None) else "path"
   for path in node.iter(tag):
      transform = getTransform(parentMap, path)
      paths = pathstr2mplpaths(path.get("d"))
      newPaths = map(lambda p: transform.transform_path(p), paths)
      mats = map(lambda p: mplpath2matrix(p), newPaths)
      pathGroups.append(zip(mats, newPaths))
   return pathGroups

"""
Get transform to be applied on node
"""
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
def pathstr2mplpaths(pathstr): # TODO: check
   paths = []
   offset = np.array([0,0])
   for match in finditer("([Mm])[^Mm]+", pathstr):
      shift = np.array([0,0]) if match.group(1).isupper() else offset
      mplpath = splinestr2mplpath(match.group(0), shift)
      if len(mplpath) > 0:
         # print mplpath.vertices
         paths.append(mplpath)
         offset = paths[-1].vertices[0,:]
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
      code, newVertices = translateToMplCmd(cmd, numbers, vertices, initialOffset)
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
def translateToMplCmd(cmd, numbers, preceding, initialOffset):
   # print cmd
   offset = preceding[-1] if len(preceding) > 0 else initialOffset
   if cmd.upper() == "M":
      points = map(np.array, partition(numbers,2))
      code = mplCode["M"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
   elif cmd.upper() == "L":
      points = map(np.array, partition(numbers,2))
      code = mplCode["L"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
   elif cmd.upper() == "H":
      xs = numbers if cmd.isupper() else reductions(add, numbers, offset[0])[1:]
      code = mplCode["L"]
      vertices = zipwith(np.array, xs, repeat(offset[1]))
   elif cmd.upper() == "V":
      ys = numbers if cmd.isupper() else reductions(add, numbers, offset[1])[1:]
      code = mplCode["L"]
      vertices = zipwith(np.array, repeat(offset[0]), ys)
   elif cmd.upper() == "Q":
      points = map(np.array, partition(numbers,2))
      code = mplCode["Q"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 2)
   elif cmd.upper() == "T":
      points = map(np.array, partition(numbers,2))
      code = mplCode["Q"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 1)
      vertices = addFirstControlforSmooth(vertices, preceding, 1)
   elif cmd.upper() == "C":
      points = map(np.array, partition(numbers,2))
      code = mplCode["C"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 3)
   elif cmd.upper() == "S":
      points = map(np.array, partition(numbers,2))
      code = mplCode["C"]
      vertices = np.array(points) if cmd.isupper() else makeAbsolute(points, offset, 2)
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
def cubic2path(matrix):
   ptstrMatrix = map(lambda part: map(lambda p: ",".join([str(x) for x in p]), part) ,np.around(matrix, 3))
   M = "M" + ptstrMatrix[0][0]
   Cs = " ".join(map(lambda part:  "C" + " ".join(part[1:]), ptstrMatrix))
   return M+" "+Cs

"""
Determine path string for quadratic Bezier spline
"""
def quadratic2path(matrix):
   ptstrMatrix = map(lambda part: map(lambda p: ",".join([str(x) for x in p]), part) ,np.around(matrix, 3))
   M = "M" + ptstrMatrix[0][0]
   Qs = " ".join(map(lambda part:  "Q" + " ".join(part[1:]), ptstrMatrix))
   return M+" "+Qs

"""
Determine path string for line group spline
"""
def lines2path(matrix):
   ptstrMatrix = map(lambda part: map(lambda p: ",".join([str(x) for x in p]), part), np.around(matrix, 3))
   return " ".join(map(lambda part:  "M " + part[0] + " L "+ part[1], ptstrMatrix))

"""
Add group to element tree 
"""
def addGroup(namespace, parent, attributes):
   tag = (namespace+"g") if (namespace is not None) else "g"
   return SubElement(parent, tag, attributes)

"""
Add path to element tree 
"""
def addCubicPath(namespace, parent, matrix, attributes):
   tag = (namespace+"path") if (namespace is not None) else "path"
   attributes["d"] = cubic2path(matrix)
   return SubElement(parent, tag, attributes)

"""
Add line to element tree 
"""
def addLine(namespace, parent, line, attributes):
   tag = (namespace +"line") if (namespace is not None) else "line"
   attributes.update(izip(["x1", "y1", "x2", "y2"], imap(str, np.array(line).flatten())))
   return SubElement(parent, tag, attributes)

"""
Add circle to element tree 
"""
def addCircle(namespace, parent, center, attributes):
   tag = (namespace+"circle") if namespace is not None else "circle"
   attributes.update(izip(["cx", "cy"], imap(str, center)))
   return SubElement(parent, tag, attributes)

"""
Add path consisting of non-connected lines to element tree 
"""
def addLinesPath(namespace, parent, matrix, attributes):
   tag = (namespace+"path") if (namespace is not None) else "path"
   attributes["d"] = lines2path(matrix)
   return SubElement(parent, tag, attributes)

"""
Add path with quadratic beziers to element tree 
"""
def addQuadraticPath(namespace, parent, matrix, attributes):
   tag = (namespace+"path") if (namespace is not None) else "path"
   attributes["d"] = quadratic2path(matrix)
   return SubElement(parent, tag, attributes)


"""
Add text at location  
"""
def addText(namespace, parent, location, text, attributes):
   tag = (namespace+"text") if namespace is not None else "text"
   attributes.update(izip(["x", "y"], imap(str, location)))
   el = SubElement(parent, tag, attributes)
   el.text = text
   return el

"""
Add text and circle at location  
"""
def addLabel(namespace, parent, location, text, attributes, textoffset = np.array([0.2, 0])):
   addCircle(namespace, parent, location, attributes)
   addText(namespace, parent, textoffset + location, text, attributes)
   return None

"""
Add polygon to element tree 
"""
def addPolygon(namespace, parent, vertices, attributes): 
   tag = (namespace +"polygon") if (namespace is not None) else "polygon"
   attributes["points"] = " ".join(map(lambda x: ",".join(map(str, x)), vertices))
   return SubElement(parent, tag, attributes)
  
"""
Write SVG data to png file
"""
def writeAsPNG(root, filename):
   from cairo   import ImageSurface, Context, FORMAT_ARGB32
   from rsvg    import Handle
   img = ImageSurface(FORMAT_ARGB32, float(root.attrib["width"]),float(root.attrib["height"])) #test
   ctx = Context(img)
   handle = Handle(None, tostring(root))
   handle.render_cairo(ctx)
   img.write_to_png(filename)
   return filename
