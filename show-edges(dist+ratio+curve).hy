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


(defn get-curvature [nodes e]
  (let [[edge-dir (apply sub (replace nodes e))]
        [path-dir (apply sub (replace nodes (, (first e) (% (inc (first e))
                                                            (len nodes)))))]]
    (cross path-dir edge-dir)))

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

              [find-edges (filter (fn [e] (pos? (get-curvature nodes e)))
                                  ind-edges)]

              [edges (replace2d nodes find-edges)]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                      line))))
      (.write tree outfile)))
