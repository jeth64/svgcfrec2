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
from detection import *


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

   # TODO: file 6: 10885 connects wrong: very short edge with similar angle: compare within a 5deg radius? but error actually occurs because of weird skeleton edge

   minAngleProb = float(options.minAngleProb)
   maxHeadEdges = int(options.maxHeadEdges) 
   maxHeadLength = float(options.maxHeadLength)
   minFreeSides = int(options.minFreeSides)
      
   groupnr=0
   nSolids = 0
   nContours = 0
   for paths in pathGroups:
      if options.verbose:
         print "\nPart", groupnr,"\n"
      
      matrices, mplpaths = zip(*paths) # unzip

      if options.maxDist is None:
         ds = sorted(pdist(np.vstack(matrices)[:,0,:], 'euclidean'))
         maxDist = ds[0]
         j=1
         while maxDist < 1: 
            maxDist = ds[j]
            j = j+1

         dn = map(lambda d, idx: np.inf if abs(idx[0]-idx[1])==1 else d, pdist(np.vstack(matrices)[:,0,:], 'euclidean'), combinations(range(len(np.vstack(matrices))),2)) # almost no difference
         print sorted(dn)[0], ds[0]
      else: maxDist = float(options.maxDist)

      
      if options.verbose:
         print "Starting discretization..."
         print "  Maximum node distance along path is", maxDist
      t1 = time()
      vertexLists = map(lambda m: discretize(m, options.discrete, maxDist), matrices)
      t2 = time()
      if options.verbose:
         print "  Calculated", len(np.vstack(vertexLists)), "discrete points"
         print "  Finished in", t2-t1, "seconds\n" 

      if options.verbose:
         print "Starting skeletonization..."
      t1 = time()
      #vD, skel, validIdx = skeleton(vertexLists, mplpaths)
      vD, skel, validIdx = skeleton2(vertexLists, mplpaths)
      t2 = time()
      if options.verbose:
         print "  Calculated", len(skel), "edges"
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
            print "  Reduced edges to", len(newSkel)
            print "  Finished in", t2-t1, "seconds\n"
     
      """
      dSite2VertDists = list(chain(*map(lambda eV, eP: [np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[0]]), np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[1]])], np.array(vD.ridge_vertices)[validIdx], np.array(vD.ridge_points)[validIdx])))
      site2siteDists = map(lambda eP: np.linalg.norm(vD.points[eP[0]]-vD.points[eP[1]]), np.array(vD.ridge_points)[validIdx])
      minWd = np.percentile(dSite2VertDists, 75) # TODO: normalize wd? take 99 percentile as maximum and set above to 1? Problem: wahrschkt 1 muss nicht sicherwedge sein

      print minWd
      print np.percentile(dSite2VertDists, 99)
      print np.percentile(dSite2VertDists, 1)
      print 0.9 * (np.percentile(dSite2VertDists, 99)-np.percentile(dSite2VertDists, 1)) 
      print map(lambda x:np.percentile(dSite2VertDists, x), range(0, 100, 10))
      print "wd", np.percentile(dSite2VertDists, 25), np.percentile(dSite2VertDists, 50), np.percentile(dSite2VertDists, 75)
      print map(lambda x:np.percentile(site2siteDists, x), range(0, 100, 10))
      print "wdP", np.percentile(site2siteDists, 25), np.percentile(site2siteDists, 50), np.percentile(site2siteDists, 75)
      """
      # Determine detection parameters

      distRangeS = []
      if not options.noSolids:
         site2vertDists = list(chain(*map(lambda eV, eP: [np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[0]]), np.linalg.norm(vD.points[eP[0]]-vD.vertices[eV[1]])], np.array(vD.ridge_vertices)[validIdx], np.array(vD.ridge_points)[validIdx])))
         distRangeS = [np.percentile(site2vertDists, 1),np.percentile(site2vertDists, 99)]
      print distRangeS
      
      if options.minWS is None:
         #site2siteDists = map(lambda eP: np.linalg.norm(vD.points[eP[0]]-vD.points[eP[1]]), np.array(vD.ridge_points)[validIdx])
         #minWS = np.percentile(site2siteDists, 15) # TODO: add percentile to detection parameters
         minWS = 0.7
      else:
         minWS = float(options.minWS)

      if options.minWC is None:
         minWC = 0.7
      else:
         minWC = float(options.minWC)
        

      if options.verbose:
         print "Starting detection..."
         print "  Minimum wedge probability concerning angles (minAngleProb) is", minAngleProb
         print "  Minimum weight for contour wedges (minWC) is", minWC
         print "  Minimum weight for solid wedges (minWS) is", minWS 
         print "  Maximum head length (maxHeadLength) is", maxHeadLength
         print "  Maximum number of head edges (maxHeadEdges) is", maxHeadEdges
         print "  Minimum number of free initial edges (minFreeSides) is", minFreeSides
         
      parameters = {"minAngleProb": minAngleProb, "noContours": options.noContours, "noSolids": options.noSolids, "maxHeadLength": maxHeadLength, "maxHeadEdges": maxHeadEdges, "minWC": minWC, "minWS": minWS, "minFreeSides": minFreeSides, "distRangeS": distRangeS}
      t1 = time()
      contourWedges, solidWedges, cycles, cycleHints=getWedges(newSkel, vD, mplpaths, options.strategy, parameters)
      t2 = time()
      if options.verbose:
         print "\n  Selected", (len(contourWedges)+len(solidWedges)), "wedges:"
         print "  -", len(contourWedges), "contour wedges"
         print "  -", len(solidWedges), "solid wedges"
         print "\n  Finished in", t2-t1, "seconds\n\n"

      nContours = nContours + len(contourWedges)
      nSolids = nSolids + len(solidWedges)
      
      if False:
         for p in matrices:
            addPath(namespace, outlayer, p, {"fill":"none", "stroke": colors[groupnr%len(colors)], "stroke-width": "0.3"} )

      if False:
         for p in matrices:      
            for mat in p:
               addCircle(namespace, outlayer, mat[0,:], {"fill": "cyan", "r": "0.4"} )

      if True:
         for v in np.vstack(vertexLists):
            addCircle(namespace, outlayer, v, {"fill": "black", "r": "0.4"} )
      
      if False: 
         for e in vD.ridge_vertices:
            if not e[0] < 0 or e[1] < 0: 
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "yellow", "stroke-width": "0.2"})
      if False:
         for e in skel:
            if not e[0] < 0 or e[1] < 0: 
               addLine(namespace, outlayer, vD.vertices[e,:], {"fill":"none", "stroke": "blue", "stroke-width": "0.3"})

      
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
           
      if False:
         for p in set(chain(*newSkel)):
            addLabel(namespace, outlayer, vD.vertices[p,:], str(p)+": " +str(np.around(vD.vertices[p,:],2)), {"fill":"magenta", "font-size":"0.5px", "r":"0.2"}, np.array([0.3, 0]))
      if False:
         for p in set(chain(*skel)):
            addLabel(namespace, outlayer, vD.vertices[p,:], str(p)+": " +str(np.around(vD.vertices[p,:],2)), {"fill":"magenta", "font-size":"0.5px", "r":"0.2"}, np.array([0.3, 0]))

      reconstruct(contourWedges, solidWedges, options, namespace, outlayer)

      groupnr = groupnr+1

   if options.verbose:
      print "Total:", nSolids+nContours, "wedges"
      print "      -", nContours, "contour wedges"
      print "      -", nSolids, "solid wedges"
   return outlayer # stub


def groupPaths(paths): # divide problem in var.hy
   p1InP2 = np.array(map(lambda p1: map(lambda p2: p1[1].contains_path(p2[1]), paths), paths))
   connComp = connected_components(p1InP2, False)
   return groupbyKeys(connComp[1], paths)

def reconstruct(contourWedges, solidWedges, options, namespace, outlayer):
   if options.svgRep == "skeleton-all": # shows all edges of skeleton
      for w in contourWedges:
         matrix = map(lambda e: vD.vertices[e,:], chain(w.headEdges, w.armEdges))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.contourStroke, "stroke-width": options.strokeWidth})
      for w in solidWedges:
         matrix = map(lambda e: vD.vertices[e,:], chain(w.headEdges, w.armEdges))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.solidStroke, "stroke-width": options.strokeWidth})
            
   elif options.svgRep == "skeleton": # shows only simplified edges
      for w in contourWedges:
         firsts = map(lambda arm: vD.vertices[arm[0],:], w.arms)
         lasts = map(lambda arm: vD.vertices[arm[-1],:], w.arms)
         head = map(lambda x,y: [x, y],firsts, np.roll(firsts, -1, 0))
         arms = map(lambda x,y: [x, y],firsts, lasts)
         matrix= np.array(list(chain(head, arms)))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.contourStroke, "stroke-width": options.strokeWidth})
      for w in solidWedges: 
         firsts = map(lambda arm: vD.vertices[arm[0],:], w.arms)
         lasts = map(lambda arm: vD.vertices[arm[-1],:], w.arms)
         head = map(lambda x: [w.deepestPoint, x], firsts)
         arms = map(lambda x,y: [x, y], firsts, lasts)
         matrix= np.array(list(chain(head, arms)))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.solidStroke, "stroke-width": options.strokeWidth})

   elif options.svgRep == "minimal": # shows only simplified edges
      for w in contourWedges: 
         matrix = np.array(map(lambda x: [w.deepestPoint,x], w.vertices))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.contourStroke, "stroke-width": options.strokeWidth})
      for w in solidWedges: 
         matrix = np.array(map(lambda x: [w.deepestPoint,x], w.vertices))
         addLinesPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.solidStroke, "stroke-width": options.strokeWidth})
         
   elif options.svgRep == "cubic":
      for w in contourWedges: 
         matrix = np.array(map(lambda v1, v2: [v1, w.deepestPoint, w.deepestPoint, v2], w.vertices, np.roll(w.vertices, -1,0)))
         addCubicPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.contourStroke, "stroke-width": options.strokeWidth})
      for w in solidWedges: 
         matrix = np.array(map(lambda v1, v2: [v1, w.deepestPoint, w.deepestPoint, v2], w.vertices, np.roll(w.vertices, -1,0)))
         addCubicPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.solidStroke, "stroke-width": options.strokeWidth})
      
   elif options.svgRep == "quadratic":
      for w in contourWedges: 
         matrix = np.array(map(lambda v1, v2: [v1, w.deepestPoint, v2], w.vertices, np.roll(w.vertices, -1,0)))
         addQuadraticPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.contourStroke, "stroke-width": options.strokeWidth})
      for w in solidWedges: 
         matrix = np.array(map(lambda v1, v2: [v1, w.deepestPoint, v2], w.vertices, np.roll(w.vertices, -1,0)))
         addQuadraticPath(namespace, outlayer, matrix, {"fill":"none", "stroke": options.solidStroke, "stroke-width": options.strokeWidth})
   return outlayer
         
if __name__ == '__main__':

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
      parser.add_option ('-D', '--maxDist', action='store', default=None, help='Maximum node distance for discretization', dest="maxDist")

      # Skeletonization
      parser.add_option ('-S', '--nosimplify', action='store', help='Turn off skeleton simplification')

      # Detection and Repulsion
      parser.add_option ('-s', '--strategy', action='store', default="contour-fill", help='strategy for wedge repulsion')
      parser.add_option ('-F', '--noSolids', action='store_true', default=False, help='determines if solid wedges are detected')
      parser.add_option ('-C', '--noContours', action='store_true', default=False, help='determines if contour wedges are detected')
      parser.add_option ('-x', '--minWS', action='store', default=None, help='minimum distance of a wedge center to the shape contour for solid wedges')
      parser.add_option ('-y', '--minWC', action='store', default=None, help='minimum triangle similarity for contour wedges')
      parser.add_option ('-z', '--minAngleProb', action='store', default="0.5", help='minimum angle probability for wedges')
      parser.add_option ('-f', '--minFreeSides', action='store', default="1", help='minimum of free head edges for selection of solid wedges in contour-fill strategies')

      # Optimization
      parser.add_option ('-e', '--maxHeadEdges', action='store', default="10", help='maximum head edges of contour wedges')
      parser.add_option ('-i', '--maxHeadLength', action='store', default="inf", help='maximum head length of contour wedges')
      
      # Representation
      parser.add_option ('-q', '--svgRep', action='store', default="cubic", help='representation of wedge in svg file')
      parser.add_option ('-w', '--strokeWidth', action='store', default="0.3", help='stroke-width of paths')
      parser.add_option ('-p', '--contourStroke', action='store', default="orange", help='stroke-color of paths')
      parser.add_option ('-P', '--solidStroke', action='store', default="magenta", help='stroke-color of paths')

      parser.add_option ('-r', '--remove', action='store_true', default=False, help='Remove input layer from document')
      parser.add_option ('-R', '--removepaths', action='store_true', default=False, help='Remove paths used for recognition from document')
      parser.add_option ('-c', '--clean', action='store_true', default=False, help='Remove all created temporary files from disk')

      (options, args) = parser.parse_args()
   
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


