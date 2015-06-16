from sys        import exit
from subprocess import call
from os.path    import splitext

def traceBitmap(infile, outfile, options):
   fileBaseName, fileExtension = splitext(infile)
   tempfile = "."+fileBaseName + ".pgm"
   try:
      if options.verbose: print "Create bitmap file for tracing..."
      call(["convert", "-flatten",infile, tempfile])
      if options.verbose: print " Done."
   except:
      print " You need 'imagemagick' to be installed to use bitmap tracing"
      exit(1)
   if options.verbose: print "Trace bitmap..."
   if options.autotrace:
      try:
         call(["autotrace",tempfile,"--output-format=svg","--output-file="+outfile,
                                    "--despeckle-level="+str(options.turdsize)])
         # on linux only if autotrace installed
      except:
         print " Using 'autotrace' for bitmap tracing failed!"
         print " Check if commandline tool 'autotrace' is installed or skip autotrace option."
         exit(1)
   else: # use potrace, generally faster say sources 
      try:
         call(["potrace",tempfile,"--svg","-o"+outfile,"-t "+str(options.turdsize)])
         # on linux only if potrace installed
      except:
         print " Using 'potrace' for bitmap tracing failed!"
         print " Check if commandline tool 'potrace' is installed"
         exit(1)
   if options.verbose: print " Done."
   if options.clean:
      if options.verbose: print "Delete bitmap file..."
      call(["rm", tempfile])
      if options.verbose: print " Done."
   return outfile

