(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn get-simple-curvature [nodes]
  "Returns [0-2-curve 1-3-curve ... (n-3)-(n-3)-curve (n-2)-0-curve (n-1)-1-curve"
  (list (map (fn [x y z] (cross (- x y) (- x z)))
             nodes
             (roll nodes -1)
             (roll nodes -2))))

(defn get-curvature [curves bounds] ;; only one way -> problem?
  (let [[gap (dec (abs (apply sub bounds)))]
        [start (min bounds)]]
    (assert (> gap 0))
    (reduce add (take gap (drop start curves)))))



;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-edges(dist+ratio+curve).svg"))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]

              [td 4.0]
              [tf 1.3]
              [dists (squareform (pdist (array nodes)))]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]
              [ind-edges (make-set (transpose (nonzero graph)))]

              [curves (get-simple-curvature nodes)]


              [find-edges (filter (fn [pair] (neg? (get-curvature curves pair))) ind-edges)]
              [edges (map (fn [pair] [(get nodes (first pair))
                                      (get nodes (second pair))])
                            find-edges)]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                      line))))
      (.write tree outfile)))

;; -> 'curvature' constraint doesnt achieve goal of eliminating concave edges and eliminates valid edges(see test4)
