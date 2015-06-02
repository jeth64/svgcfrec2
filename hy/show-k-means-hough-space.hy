(import [xml.etree.ElementTree [parse SubElement]]
        [var [*]]
        [svg [*]]
        [bezier [*]]
        [hough [*]]
        [kmeans [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]]
        [scipy.sparse.csgraph [connected_components]]
        [matplotlib.pyplot [figure xlabel ylabel title hist2d]]
        [pylab [savefig]])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn config [key]
  (case key
        :td 4.0
        :tf 1.3
        :theta-range [0]
        :hough-bins [180 200]
        :axis [[0 360] [-70 70]]))



;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv outfile (+ filebasename "-k-means-components-hough.svg"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]
              [node-dirs-v (squeeze (get-node-dirs-v new-path :gradient :forward))]
              [node-dirs (squeeze (get-node-dirs new-path :gradient :bi))]
              [n (len nodes)]

              [a (print (shape node-dirs))]

              [td (config :td)]
              [tf (config :tf)]
              [dists (get-pdists nodes)]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [find-edges (filter (fn [e] (inner-edge? nodes node-dirs-v e))
                                  ind-edges)]
              [edges (list (replace2d nodes find-edges))]
              [edge-means (mean (array (list edges)) 1)]

              [thetas (list (map (fn [pair] (hstack [(get node-dirs (first pair))
                                                     (get node-dirs (second pair))]))
                                 ind-edges))]

              [a (print "thetas" (shape thetas))]

              [cc (second (connected_components graph False))]
              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc (range n)))))) ]
              [grouped-ind-edges (list (map (fn [x] (.intersection ind-edges
                                                                   (make-set (combinations x 2))))
                                        conn-comp))]
              [grouped-edges (map (fn [g] (list (replace2d nodes g)))
                                  grouped-ind-edges)]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [point edge-means] (add-circle root {"fill" "orange" "r" "0.2"} point))
          (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                      line))
          (for [ge grouped-edges]
            (let [[edge-means (mean (list ge) 1)]
                  [hough-pts (hough-transform edge-means thetas)]
                  [a (print "h-pts" (shape hough-pts))]]

              (for [point (get-thetas hough-pts 1 30)]
                (add-hough-line root {"fill" "none" "stroke" "yellow" "stroke-width" "0.1"} point))

              ))))
      (.write tree outfile)))
