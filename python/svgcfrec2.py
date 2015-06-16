import numpy as np
from sys                    import exit
from optparse               import OptionParser, TitledHelpFormatter
from time                   import time, asctime
from os.path                import splitext
from scipy.sparse.csgraph   import connected_components
from operator import add

# own modules
from bitmaptrace import traceBitmap
from svg import *
from discrete import *
from voronoi import *


def convertFromSVG(svgData, filename): #TODO:Test
   from cairo   import ImageSurface, Context, FORMAT_ARGB32
   from rsvg    import Handle
   img = ImageSurface(FORMAT_ARGB32, 640,480)
   ctx = Context(img)
   handle = Handle(None, str(svgData))
   handle.render_cairo(ctx)
   img.write_to_png(filename)
   return filename

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

   if  options.clean and fileExtension is not ".svg":
      if options.verbose: print "Delete created input svg..."
      call(["rm", filename])
      if options.verbose: print " Done.\n"


   return tree

def findCuneiformsInSVG(tree, inlayer, outlayer, namespace, options):
   pathGroups = getPaths(tree, inlayer, namespace, options.unclean)

   if options.group:
      pathGroups = groupPaths(list(chain(pathGroups)))
   
   for paths in pathGroups:
      matrices, mplpaths = zip(*paths) # should unzip
      compoundShape = Path.make_compound_path(*mplpaths) # for point testing, test if exact path is used or vertices as polygons
      
      vertices = np.vstack(map(lambda m: discretize(m, options.discrete, 3.0), matrices))

      for p in matrices:
         addPath(namespace, outlayer, p, {"fill":"none", "stroke": "red", "stroke-width": "0.3"} )

      #for v in vertices:
         #addCircle(namespace, outlayer, v, {"fill": "blue", "r": "0.2"} )
      
   return outlayer # stub


def groupPaths(paths): # divide problem in var.hy
   p1InP2 = np.array(map(lambda p1: map(lambda p2: p1[1].contains_path(p2[1]), paths), paths))
   connComp = connected_components(p1InP2, False)
   return groupbykeys(connComp, paths)
   

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
      
      # Representation
      parser.add_option ('-w', '--strokewidth', metavar="WIDTH", action='store', default="0.5", help='stroke-width of paths', dest="strokewidth")
      parser.add_option ('-s', '--strokecolor', metavar="COLOR", action='store', default="blue", help='stroke-color of paths', dest="strokecolor")

      parser.add_option ('-r', '--remove', action='store_true', default=False, help='Remove input layer from document')
      parser.add_option ('-c', '--clean', action='store_true', default=False, help='Remove all created temporary files from disk')

      (options, args) = parser.parse_args()

      # gapSize = float(options.gapSize)
   
      pathAttributes = {"fill":"none", "stroke": options.strokecolor, "stroke-width": options.strokewidth}

      if options.verbose: print asctime()
      if len(args) < 1:
         print "No input file given"
         exit(0)

      tree = findCuneiforms(args[0], options)

      if options.verbose: print "Writing result to: ", options.outputfile, "\n"
      tree.write(options.outputfile)

      if options.verbose: print asctime()
      if options.verbose: print 'TOTAL TIME IN SECONDS:',
      if options.verbose: print (time() - start_time)
      exit(0)
   except KeyboardInterrupt, e:
      raise e
   except SystemExit, e:
      raise e

