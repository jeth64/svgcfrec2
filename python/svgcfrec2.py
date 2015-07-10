import numpy as np
from sys                    import exit
from optparse               import OptionParser, TitledHelpFormatter
from time                   import time, asctime
from os.path                import splitext
from scipy.sparse.csgraph   import connected_components
from xml.etree.ElementTree  import tostring
from operator import add
import warnings
from math import ceil, floor
from matplotlib.path import Path
from scipy.spatial.distance import pdist, squareform

# own modules
from bitmaptrace import traceBitmap
from svg import *
from discrete import *
from skeleton import *
from geometry import *


def findCuneiforms(infile, options):
   
   if options.verbose: print "\nInput: ", infile,"\n"
   fileBaseName, fileExtension = splitext(infile)

   if fileExtension != ".svg":
      filename = "."+fileBaseName + ".svg"
      filename = traceBitmap(infile, filename, options)
   else:
      filename = infile

   tree, root, ns = getSVG(filename)

   outlayer = root.find(".//*[@id='"+options.outLayerID+"']")
   if outlayer is None:
      outlayer = addGroup(ns, root, {"id":options.outLayerID})
         
   if options.layerID is None:
      findCuneiformsInSVG(tree, root, outlayer, ns, options)
   else:
      inlayer = root.find(".//*[@id='"+options.layerID+"']")
      if inlayer is None:
         print "No layer with id='"+options.layerID+"' found"
         print "Program aborting..."
         exit(1)
      
      cf = findCuneiformsInSVG(tree, inlayer, outlayer, ns, options)

   if options.remove and options.layerID is not None:
      if options.verbose: print "Delete input layer for recognition..."
      deleteElement(root,inlayer)
      if options.verbose: print " Done.\n"

   if options.removepaths:
      if options.verbose: print "Delete paths used for recognition..."
      if options.layerID is None:
         deletePaths(tree, ns, root)
      else: deletePaths(tree, ns, inlayer)
      if options.verbose: print " Done.\n"

   if  options.clean and fileExtension is not ".svg":
      if options.verbose: print "Delete created input svg..."
      call(["rm", filename])
      if options.verbose: print " Done.\n"

   return tree

colors = ["red", "magenta", "cyan", "green" ]

def findCuneiformsInSVG(tree, inlayer, outlayer, namespace, options):
   pathGroups = getPaths(tree, inlayer, namespace, options.unclean)
   if options.group:
      pathGroups = groupPaths(list(chain(*pathGroups)))

   i=0
   for paths in pathGroups:
      matrices, mplpaths = zip(*paths) # unzip
      
      obj = Path.make_compound_path(*mplpaths).get_extents().get_points()
      height = obj[1,1]-obj[0,1] # height or width?

      maxHeadLength = 4*height # avoid using, oder:longest edge of simplified mal 3 oder 4
      
      ds = sorted(pdist(np.vstack(matrices)[:,0,:], 'euclidean'))
      maxDist = ds[0]
      i=1
      while maxDist < 10**(-3):
         maxDist = ds[i]
         i = i+1
      
      if options.verbose:
         print "Starting discretization..."
         print "  Maximum node distance along path is", maxDist
      t1 = time()
      vertexLists = map(lambda m: discretize(m, options.discrete, maxDist), matrices)
      t2 = time()
      if options.verbose:
         print "  Finished in", t2-t1, "seconds\n" 

      if options.verbose:
         print "Starting skeletonization..."
      t1 = time()
      vD, skel, validIdx = skeleton(vertexLists, mplpaths)
      t2 = time()
      if options.verbose:
         print "  Finished in", t2-t1, "seconds\n"
         
      if options.nosimplify:   
         newSkel = skel
      else:
         if options.verbose:
            print "Starting skeleton simplification..."
         t1 = time()
         newSkel = simplifySkeleton(skel, vD.vertices, mplpaths)
         t2 = time()
         if options.verbose:
            print "  Finished in", t2-t1, "seconds\n"

      if True: #paradigm == "contour-fill":
         
         if options.verbose:
            print "Detecting contour wedges..."
         t1 = time()
         wedgeArms1, wedgeHeads1, cycles, cycleHints, wedgeCycles = detectContourWedges(newSkel, skel, vD, mplpaths, maxHeadEdges=10,maxHeadLength=maxHeadLength)
         t2 = time()
         if options.verbose:
            print "  Finished in", t2-t1, "seconds\n"

         dSite2VertDists = list(chain(*map(lambda eV, eP: [np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[0]]), np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[1]])], \
                                           np.array(vD.ridge_vertices)[validIdx], np.array(vD.ridge_points)[validIdx])))
         minWd = np.percentile(dSite2VertDists, 75)
         minFreeSides = 1
         usedVerts = set(chain(chain(*wedgeCycles)))
         usedEdges = set([])
         if options.verbose:
            print "Detecting solid wedges..."
            print "  Minimum joint distance to contour is", minWd
            print "  Minimum number of initial free edges is", minFreeSides
         t1 = time()
         
         wedges = detectSolidWedges(newSkel, usedVerts, usedEdges, vD, mplpaths, minWd, minFreeSides)
         t2 = time()
         if options.verbose:
            print "  Finished in", t2-t1, "seconds\n"
        

      
      if True:
         for p in matrices:
            addPath(namespace, outlayer, p, {"fill":"none", "stroke": colors[i%len(colors)], "stroke-width": "0.3"} )
      i = i+1

      if True:
         for p in matrices:      
            for mat in p:
               addCircle(namespace, outlayer, mat[0,:], {"fill": "cyan", "r": "0.4"} )

      if True:
         for v in np.vstack(vertexLists):
            addCircle(namespace, outlayer, v, {"fill": "blue", "r": "0.2"} )
      
      if False: 
         for e in vD.ridge_vertices:
            if not e[0] < 0 or e[1] < 0: 
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "yellow", "stroke-width": "0.3"})
      if False:
         for e in skel:
            if not e[0] < 0 or e[1] < 0: 
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "cyan", "stroke-width": "0.3"})

      
      if False: # original skeleton
         verts = list(chain(*map(lambda mat: mat[:,0,:], matrices)))
         vD2, skel2, validIdx2 = skeleton(verts, mplpaths)
         
         for e in vD2.ridge_vertices:
            if -1 not in e:
               addLine(namespace, outlayer, vD2.vertices[e,:], {"fill":"none", "stroke": "green", "stroke-width": "0.3"})
         for e in skel2:
            addLine(namespace, outlayer, vD2.vertices[e,:], {"fill":"none", "stroke": "cyan", "stroke-width": "0.3"})
            
            
      if True:
         for e in newSkel:
            addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "blue", "stroke-width": "0.3"})
      
      if False:
         for c in cycles:
            poly = vD.vertices[c,:]
            addPolygon(namespace, outlayer, poly, {"fill":"none", "stroke": "magenta", "stroke-width": "0.3"})
      if False: 
         for e in cycleHints:
            addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "cyan", "stroke-width": "0.3"})

      if True:
         for c in wedgeHeads1:
            poly = vD.vertices[c,:]
            addPolygon(namespace, outlayer, poly, {"fill":"none", "stroke": "yellow", "stroke-width": "0.3"})
      
      if False:
         for a in chain(*wedgeArms1):
            for e in zip(a[:-1],a[1:]):
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "grey", "stroke-width": "0.3"})

      if True:
         for a in chain(*wedgeArms1):
            addLine(namespace, outlayer, vD.vertices[[a[0],a[-1]],:], {"fill":"none", "stroke": "orange", "stroke-width": "0.3"})
               
      if False:
         for c,w in zip(cycles,ws1): # can be tried with thresholding or conflict sets (conflict worse with longer common edges)
            if w < 0.2:
               poly = vD.vertices[c,:]
               addPolygon(namespace, outlayer, poly, {"fill":"none", "stroke": "red", "stroke-width": "0.3"})


      if True:
         for w in wedges:
            print "arms"
            for e in w.armEdges:
               print e
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "orange", "stroke-width": "0.3"})
            print "head"
            for e in w.headEdges:
               print e
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "yellow", "stroke-width": "0.3"})
            for v in w.vertices:
               addLine(namespace, outlayer, [w.deepestPoint, v], {"fill":"none", "stroke": "magenta", "stroke-width": "0.3"})
           
      if True:
         for p in set(chain(*newSkel)):
            addLabel(namespace, outlayer, vD.vertices[p,:], str(p)+": " +str(np.around(vD.vertices[p,:],2)), {"fill":"magenta", "font-size":"0.5px", "r":"0.2"}, np.array([0.3, 0]))
      
   return outlayer # stub


def groupPaths(paths): # divide problem in var.hy
   p1InP2 = np.array(map(lambda p1: map(lambda p2: p1[1].contains_path(p2[1]), paths), paths))
   connComp = connected_components(p1InP2, False)
   return groupbyKeys(connComp[1], paths)
   

if __name__ == '__main__':
   floatType = np.float32

   try:
      start_time = time()
      parser = OptionParser(formatter=TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
      parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
      parser.add_option ('-o', '--ofile', metavar="FILE", action='store', default="out.svg", help='output file',
                                          dest="outputfile")
      parser.add_option ('-l', '--layerid', metavar="ID", action='store', help='ID of SVG layer to find cuneiforms in, if not given, all layers are considered', dest="layerID")
      parser.add_option ('-L', '--outputlayerid', metavar="ID", action='store', default="cuneiforms", help='ID of SVG layer to save cuneiforms to', dest="outLayerID")

      # Bitmap tracing
      parser.add_option ('-a', '--autotrace', action='store_true', default=False, help='Use AutoTrace instead of potrace for bitmap tracing')
      parser.add_option ('-t', '--turdsize', action='store', default=0, help='For bitmap tracing: suppress speckles of sizes smaller than this')
      # same with potrace and autotrace?
        
      # Interpretation
      parser.add_option ('-g', '--group', action='store_true', default=False, help='Automate grouping of paths')
      parser.add_option ('-u', '--unclean', action='store_true', default=False, help='Simplify paths before continuing')

      # Discretization
      parser.add_option ('-d', '--discrete', action='store', default="polypath", help='Discretization method')

      # Skeletonization
      parser.add_option ('-S', '--nosimplify', action='store', help='Turn off skeleton simplification')

      # Detection
      parser.add_option ('-p', '--paradigm', metavar="PARADIGM", action='store', default="contour-fill", help='paradigme for wedge detectio', dest="strokecolor")
      
      # Representation
      parser.add_option ('-w', '--strokewidth', metavar="WIDTH", action='store', default="0.5", help='stroke-width of paths', dest="strokewidth")
      parser.add_option ('-s', '--strokecolor', metavar="COLOR", action='store', default="blue", help='stroke-color of paths', dest="strokecolor")

      parser.add_option ('-r', '--remove', action='store_true', default=False, help='Remove input layer from document')
      parser.add_option ('-R', '--removepaths', action='store_true', default=False, help='Remove paths used for recognition from document')
      parser.add_option ('-c', '--clean', action='store_true', default=False, help='Remove all created temporary files from disk')

      (options, args) = parser.parse_args()

      # gapSize = float(options.gapSize)
   
      pathAttributes = {"fill":"none", "stroke": options.strokecolor, "stroke-width": options.strokewidth}

      if options.verbose: print asctime()
      if len(args) < 1:
         print "No input file given"
         exit(0)

      with warnings.catch_warnings(): # to suppress visible deprecation warning of np.rank
         warnings.simplefilter("ignore") 
         tree = findCuneiforms(args[0], options)

      if options.verbose: print "\nWriting result to: ", options.outputfile, "\n"
      tree.write(options.outputfile)

      if options.verbose: print asctime()
      if options.verbose: print 'TOTAL TIME IN SECONDS:',
      if options.verbose: print (time() - start_time)
      exit(0)
   except KeyboardInterrupt, e:
      raise e
   except SystemExit, e:
      raise e


