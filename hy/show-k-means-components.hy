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
;  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (range 1 5)))
  (list (map (fn [x] (.join "" ["testfiles_size/test" (.__str__ x) ".svg"]))
             (reversed (list (range 1 10)))))
;  ["more-strings.svg"]
  )

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
      (setv outfile (+ filebasename "-k-means-components.svg"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [a (print filebasename)]

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
              [a (print "ie" (shape (list ind-edges)) "fie" (shape (list find-edges))
                        "e" (shape edges))]
              [a (print 1 )]
              [edge-means (if (< 0 (len edges))
                            (mean (array (list edges)) 1)
                            [])]
              [a (print 1)]

              [cc (second (connected_components graph False))]
              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc (range n)))))) ]
              [grouped-ind-edges (list (map (fn [x] (.intersection find-edges
                                                                   (make-set (combinations x 2))))
                                        conn-comp))]
              [grouped-edges (map (fn [g] (list (replace2d nodes g)))
                                  grouped-ind-edges)]]
          (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
          (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
          (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
          (for [point edge-means] (add-circle root {"fill" "orange" "r" "0.2"} point))
        ;  (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"} line))

          (for [ge grouped-edges] ; Fehler?
            (let [[ge-means (mean (array (list ge)) 1)]
                  [(, betas groups e) (ind-get-thetas ge 0.3 100)]
                  [hough-pts (list (map (fn [b] (let [[theta (arctan (/ (second b) (first b)))]]
                                                  [(/ (cos theta) (first b)) theta]))
                                        betas))]
                  [alphas (list (map (fn [b] (- 90 (rad2deg (arctan (/ (second b) (first b))))))
                                     betas))]
                  [boxes (list (map (fn [g alpha] (bounding-box (list (replace ge-means g)) alpha))
                                    groups alphas))]
                  ]
           ;   (print " l: " (len groups))
           ;   (print " g: " groups e)
           ;   (print "alphas " alphas)
              (for [box boxes]
               ; (add-polygon root {"fill" "none" "stroke" "red" "stroke-width" "0.1"} box)
                (add-line root {"fill" "none" "stroke" "yellow" "stroke-width" "0.1"}
                          [(mean [(get box 0) (get box 1)] 0)
                           (mean [(get box 2) (get box 3)] 0)]))


             ; (print (shape betas))
             ; (print "hpt" hough-pts)
            ;  (for [point hough-pts] (add-hough-line root {"fill" "none" "stroke" "yellow" "stroke-width" "0.1"} point))
              ))))
      (.write tree outfile)))
