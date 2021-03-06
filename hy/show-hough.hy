(import [xml.etree.ElementTree [parse SubElement]]
        [svg [*]]
        [bezier [*]]
        [hough [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]]
        [scipy.spatial.distance [pdist squareform]]
        [matplotlib.pyplot [*]]
        [pylab [savefig]]
        )

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn config [key]
  (case key
        :td 4.0
        :tf 1.3
        :hough-bins [180 200]
        :axis [[0 360] [-70 70]]))

;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv outfile (+ filebasename "-hough.png"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]

              [td (config :td)]
              [tf (config :tf)]
              [dists (squareform (pdist (array nodes)))]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [edges (replace2d nodes ind-edges)]

              [edge-means (mean (array (list edges)) 1)]

              [hough-pts (hough-transform edge-means)]
              [pts (transpose (list (full-hough-range (list hough-pts))))]]
          (figure)
          (title (+ "Hough transform of " filebasename))
          (ylabel "Distance from origin (r)")
          (xlabel "Angle (theta in °)")
          (hist2d (first pts) (second pts) (config :hough-bins) (config :axis))
          (savefig outfile)))))
