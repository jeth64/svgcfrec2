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
        :hough-bins [180 200]
        :axis [[0 360] [-70 70]]))

;;-> dann rotated box

;; beta entspricht: alpha/gamma, beta/gamma oder cos(theta)/r, sin(theta)/r


;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv outfile (+ filebasename "-k-means.svg"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]
              [node-dirs (squeeze (get-node-dirs-v new-path :gradient :forward))]

              [td (config :td)]
              [tf (config :tf)]
              [dists (get-pdists nodes)]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [find-edges (filter (fn [e] (inner-edge? nodes node-dirs e))
                                  ind-edges)]
              [edges (list (replace2d nodes find-edges))]

              [edge-means (mean (array (list edges)) 1)]

              [(, betas groups e) (ind-get-thetas edges 0.2 30)]
              [a (print filebasename " l: " (len groups))]
           ;   [a (print groups)]
              [hough-pts (map (fn [b] (let [[theta (arctan (/ (second b) (first b)))]]
                                        [(/ (cos theta) (first b)) theta]))
                              betas)]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [point edge-means] (add-circle root {"fill" "orange" "r" "0.2"} point))
          (for [point hough-pts] (add-hough-line root {"fill" "none" "stroke" "yellow" "stroke-width" "0.1"} point))
          (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                      line))))
      (.write tree outfile)))
