(import [xml.etree.ElementTree [parse SubElement]]
        [var [*]]
        [svg [*]]
        [bezier [*]]
        [hough [*]]
        [kmeans [*]]
        [wedges [*]]
        [geometry [*]]
        [distance [*]]
        [tools [*]]
        [numpy [*]]
        [numpy :as np]
        [operator [add sub]]
        [itertools [combinations]]
        [scipy.sparse.csgraph [connected_components]]
        [matplotlib.pyplot [figure xlabel ylabel title hist2d]]
        [pylab [savefig]])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (range 1 5)))
;;  (list (map (fn [x] (.join "" ["testfiles_size/test" (.__str__ x) ".svg"])) (reversed (list (range 1 11)))))
  ;;  ["more-strings.svg"]
  ;["testfiles/test2.svg"]
)

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn config [key]
  (case key
        :td 4.0
        :tf 1.3
        :hough-bins [180 200]
        :axis [[0 360] [-70 70]]))

(if (= __name__ "__main__"))
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv outfile (+ filebasename "-entire-test.svg"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (setv all-lines [])
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]
              [node-dirs (squeeze (get-node-dirs-v new-path :gradient :forward))]
              [n (len nodes)]

              [td (config :td)]
              [tf (config :tf)]
              [dists (get-pdists nodes)]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [find-edges (set (filter (fn [e] (inner-edge? nodes node-dirs e))
                                       ind-edges))]
              [edges (list (replace2d nodes find-edges))]

              [cc (second (connected_components graph False))]
              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc (range n)))))) ]
              [grouped-ind-edges (list (remove empty?
                                               (map (fn [x]
                                                      (.intersection find-edges
                                                                     (make-set (combinations x 2))))
                                                    conn-comp)))]
              [grouped-edges (list (map (fn [g] (list (replace2d nodes g)))
                                        grouped-ind-edges))]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          ;;(for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"} line))

          (print filebasename)
          (for [ge grouped-edges]
            (let [[ge-means (mean (array ge) 1)]
                  [(, alphas groups e) (get-lines ge 1.0 10)]

                  [boxes (list (map (fn [g alpha] (bounding-box (list (replace ge-means g)) alpha))
                                    groups alphas))]
                  ;; different strategy when n < 3/4/5? (manchmal box-ecken gleicher punkt)

                  [lines (list (map (fn [box] [(mean [(get box 0) (get box 1)] 0)
                                               (mean [(get box 2) (get box 3)] 0)])
                                boxes))]
                  ]
              (.extend all-lines lines)
              (for [point ge-means] (add-circle root {"fill" "orange" "r" "0.2"} point))
              (for [line lines]
                (add-line root {"fill" "none" "stroke" "yellow" "stroke-width" "0.1"} line)))))

      ;  (for [point (get-isecs all-lines)] (add-circle root {"fill" "red" "r" "0.2"} point))
        (setv wedges (get-wedges all-lines))
        (for [wedge wedges]
          ;(add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.5"} (first wedge))
          (for [isect (second wedge)] (add-circle root {"fill" "red" "r" "0.2"} isect))
          (add-circle root {"fill" "firebrick" "r" "0.2"} (mean (second wedge) 0))
          (for [vertex (last wedge)]
            (add-circle root {"fill" "magenta" "r" "0.2"} vertex)
            None)
          None)
        (setv chosen (set (reduce-wedges (list (map last wedges)))))
        (for [wedge chosen]
          (for [vertex wedge]
            (add-circle root {"fill" "cyan" "r" "0.2"} vertex)
            None)
          None)
        (setv ws (list (filter (fn [w] (in  (tuple (map tuple (last w))) chosen)) wedges)))
        (for [wedge ws]
          (add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.5"} (first wedge))
          None))
      None
      (.write tree outfile)
      None))
