(import [xml.etree.ElementTree [parse SubElement]]
;        [misc [*]]
        [var [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [graph [*]]
        [tools [*]]
        [numpy [*]]
        [time [time]]
        [operator [add sub]]
        [scipy.spatial [Voronoi]]
        [scipy.sparse.csgraph [connected_components]]
        [matplotlib.path [Path]]
        [itertools :as it])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles_whole/test" (.__str__ x) ".svg"]))
             (reversed (list (range 4 5))))))

;; (defn prepare2 [path] (split-segments-to2 3.0 1.0 path)) ;; unsinnig

;;(defn prepare [path] (split-segments-to "distance" 3.0 (split-segments-to "mid-dist" 1.0 path)))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(def config (dict [[:min-path-node-number 0] ;; ignore paths with few nodes
                   [:simplify True]
                   [:voronoi-cells False]
                   [:labels False]]))

;;main
(for [infile testfiles] ;; test
  (do (print "\nProcessing" infile "...\n")
      (setv start-time (time))
      (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-voronoi.svg"))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))

      (print "Calculating path matrices..." )
      (setv paths-list (get-paths2 root namespace))
      (print "Took..." (- (time) start-time) "seconds\n")

      (for [paths paths-list];; for each svg path node ...
        (setv time8 (time))
     ;;   (print "" time8 (do (map prepare2 paths) None) (- (time) time8)) ;;2.4e5 +2.3e5
        ;; unterschied bei test2 wholefiles: 3.5e-5 zu 9e-6 nur aus ordnung zu verdanken: letztes immer schneller
        (print "" time8 (do (map prepare paths) None) (- (time) time8)) ;; 2.3e5 2.3e5

        (setv mpl-paths (list (map (fn [p] (Path (list (map first p))))
                                   (filter (fn [p] (> (len p) (:min-path-node-number config)))
                                           (map prepare paths)))))

        (setv (, n labels) (divide-problem mpl-paths))

        (for [current-group (range n)]
          (let [[cur-paths (list (replace mpl-paths (first (where (= current-group (array labels))))))]

                [nodes (vstack (list (map (fn [x] x.vertices) cur-paths)))]

              ;;  [a (print "Calculating Voronoi Tesselation..." )]
                [vor (Voronoi nodes)]
           ;;     [a (print "Took..." (- (time) start-time) "seconds")]
            ;;    [a (print (len vor.ridge_vertices) "ridges in total\n")]

                ;; filter segments and merge when there is no intersections
          ;;      [a (print "Filter and merge segments..." )]
                [(, filtered-edges ridge-points) (filter-edges3 vor cur-paths)]
                [ind-edges (if (:simplify config)
                             (simplify-graph filtered-edges vor cur-paths)
                             filtered-edges)]

           ;;     [a (print "Took..." (- (time) start-time) "seconds")]
            ;;    [a (print (len filtered-edges) "edges found")]
          ;;      [a (print (len ind-edges) "edges left after simplifying\n")]

                ;; indices to coordinates
                [edges (list (replace2d vor.vertices ind-edges))]
                [all-edges (list (replace vor.vertices
                                          (list (remove (fn [e] (any (= (array e) -1)))
                                                        vor.ridge_vertices))))]]
            ;;  (for [point nodes] (add-circle root {"fill" "grey" "r" "0.2"} point))
            (print (len edges))
            (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} line))
            (when (:voronoi-cells config)
              (for [line all-edges]
                (add-line root {"fill" "none" "stroke" "black" "stroke-width" "0.1"} line)))
            (when (:labels config)
              (for [(, point text) (map (fn [i c] [c (str i)])
                                        (range (len vor.vertices))
                                        vor.vertices)]
                (add-text root {"fill" "magenta" "font-size" "0.3px" } point text))))))
      (.write tree outfile)
      (print "\n TOTAL TIME IN SECONDS:\n" (- (time) start-time))))
