(import [xml.etree.ElementTree [parse SubElement]]
        [svg [*]]
        [bezier [*]]
        [hough [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]]
        [itertools [combinations]]
        [scipy.spatial.distance [pdist squareform]]
        [scipy.sparse.csgraph [connected_components]]
        [matplotlib.pyplot [*]]
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
        :hough-bins [180 200]
        :axis [[0 360] [-70 70]]))


(defn angle [v]
  "Angle of v in degrees as in polar coordinate system"
  (% (- (degrees (apply arctan2 (list (reversed v)))) 90) 180))

(defn get-node-dirs [path &optional [angle-range [0]]]
  "Get polar angle of curve gradients of nodes in path with a given range of angles around result"
  (let [[node-dirs (list (map (fn [segment]
                                (list (map (fn [x] (% (+ (angle (apply sub (slice segment 0 None 3)))
                                                         x)
                                                      180))
                                           angle-range)))
                              path))]]
    (hstack [node-dirs (roll node-dirs 1)])))

;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]
              [node-dirs (get-node-dirs new-path)]

              [n (len nodes)]

              [td (config :td)]
              [tf (config :tf)]
              [dists (squareform (pdist (array nodes)))]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [edges (list (replace2d nodes ind-edges))]

              [cc (second (connected_components graph False))]

              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc (range n)))))) ]


              [grouped-ind-edges (list (map (fn [x] (.intersection ind-edges
                                                                   (make-set (combinations x 2))))
                                        conn-comp))]
              [grouped-edges (map (fn [g] (list (replace2d nodes g)))
                                  grouped-ind-edges)]]
          (setv i 0)
          (for [ge grouped-edges]
            (let [[edge-means (mean (list ge) 1)]
                  [hough-pts (hough-transform edge-means)]
                  [pts (transpose (list (full-hough-range (list hough-pts))))]]
              (figure)
              (title (+ "Hough transform of " filebasename ", Component " (str i)))
              (ylabel "Distance from origin (r)")
              (xlabel "Angle (theta in °)")
              (hist2d (first pts) (second pts) (config :hough-bins) (config :axis))
              (savefig (+ filebasename "-hough(comp" (str i) ").png")))
            (setv i (inc i))))))))
