(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]]
        [itertools [combinations]]
        [matplotlib.pyplot [*]]
        [pylab [savefig]]
        )

(defn show-hough [points name file]
  (let [[pts (transpose (list points))]]
    (do (figure)
        (title name)
        (ylabel "Distance from origin (r)")
        (xlabel "Angle (theta in °)")
        (hist2d (first pts) (second pts) [180 200])
        (savefig file))
    None))

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn hough-transform [nodes &optional theta-list r-range]
  (let [[points (if theta-list
                  (map (fn [p thetas]
                         (list (apply concat
                                      (map (fn [theta]
                                             (list (map (fn [offset]
                                                          [theta
                                                           (+ (* (first p) (cos (deg2rad theta)))
                                                              (* (second p) (sin (deg2rad theta)))
                                                              offset)])
                                                        (if r-range r-range [0]))))
                                           thetas))))
                       nodes
                       theta-list)
                  (map (fn [theta] (map (fn [p] [theta
                                                 (+ (* (first p) (cos (deg2rad theta)))
                                                    (* (second p) (sin (deg2rad theta))))])
                                        nodes))
                       (range 0 180)))]]
    (apply concat points)))

(defn angle [v]
  "Angle of v in degrees as in polar coordinate system"
  (% (- (degrees (apply arctan2 (list (reversed v)))) 90) 180))

(defn get-node-dirs [path &optional angle-range]
  "Get polar angle of curve gradients of nodes in path with a given range of angles around result"
  (let [[node-dirs (list (map (fn [segment]
                                (list (map (fn [x] (% (+ (angle (apply sub (slice segment 0 None 3)))
                                                         x)
                                                      180))
                                           (if angle-range angle-range [0]))))
                              path))]]
    (hstack [node-dirs (roll node-dirs 1)])))

(defn replace [coll keys]
  "Replaces given keys by items in collection"
  (map (fn [x] (get coll x)) keys))

(defn replace2d [coll key-matrix]
  "Replaces keys in matrix by items in collection"
  (map (fn [key] (list (replace coll key))) key-matrix))

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
              [node-dirs (get-node-dirs new-path
                                        ;(range -15 16)
                                        )]

              [n (len nodes)]


              [td 4.0]
              [tf 1.3]
              [dists (squareform (pdist (array nodes)))]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]


              [ind-edges (make-set (transpose (nonzero graph)))]
              [edges (list (replace2d nodes ind-edges))]

              [a (print (len edges) (first edges))]

              [thetas (map (fn [pair] (hstack [(get node-dirs (first pair))
                                               (get node-dirs (second pair))]))
                           ind-edges)]

              [edge-means (mean (array (list edges)) 1)]

              [hough-pts (list (hough-transform edge-means
                                                ;thetas (arange -1 1 0.05)
                                                ))]
              [cc (second (connected_components graph False))]

            ;  [a (print cc)]
              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc
                                                                (range n)))))) ]

            ;  [a (print conn-comp)]
              [separate-2edge-comp (group-by (map (fn [y] (= (len y) 2))
                                                    conn-comp)
                                               conn-comp)]

              [two-edge-comp (if (in True separate-2edge-comp)
                               (list (get separate-2edge-comp True))
                               [])]

;              [a (print two-edge-comp)]

              [multi-edge-comp (if (in False separate-2edge-comp)
                                 (list (get separate-2edge-comp False))
                                 [])]

            ;  [a (for [x multi-edge-comp] (print x))]

              [grouped-ind-edges (list (map (fn [x] (.intersection ind-edges
                                                               (make-set (combinations x 2))))
                                        multi-edge-comp))]
              [grouped-edges (map (fn [g] (list (replace2d nodes g)))
                                  grouped-ind-edges)]
            ;  [a (print "lge" (len grouped-edges) "fge" (first grouped-edges) )]

              ]

          (setv i 0)
          (for [ge grouped-edges]
            (show-hough (hough-transform (mean (array (list ge)) 1))
                        (+ "Hough transform of " filebasename ", Component " (str i))
                        (+ filebasename "-hough(comp" (str i) ").png"))
            (setv i (inc i)))

          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                      line))
          (for [point edge-means] (add-circle root {"fill" "orange" "r" "0.2"} point))

          ))
    ;;  (.write tree outfile)

      ))
