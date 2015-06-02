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
        [matplotlib.path [Path]]
        [itertools [chain]]
        [itertools :as it]
        )

(def testfiles
;;  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (range 2 3)))
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (reversed (list (range 6 7)))))
  )

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(def config (dict [[:min-path-node-number 0] ;; ignore paths with few nodes
                   [:simplify True]
                   [:max-hole-length 50000.0] ;500 or 1000, critical
                   [:max-hole-connectors 8]
                   [:hole-indications False]
                   [:repulse-edges True]]))
;;main
(for [infile testfiles]
  (do (print "\nProcessing" infile "...\n")
      (setv start-time (time))
      (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-holes.svg"))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))

      (print "Calculating path matrices..." )
      (setv paths-list (get-paths2 root namespace))
      (print "Took..." (- (time) start-time) "seconds\n")

      (for [paths paths-list];; for each svg path node ...
        (setv mpl-paths (list (map (fn [p] (Path (list (map first p))))
                                   (filter (fn [p] (> (len p) (:min-path-node-number config)))
                                           (map prepare paths)))))
        (setv (, n labels) (divide-problem mpl-paths))

        (for [current-group (range n)]
          (let [[cur-paths (list (replace mpl-paths (first (where (= current-group
                                                                     (array labels))))))]
                [nodes (vstack (list (map (fn [x] x.vertices) cur-paths)))]

                [a (print "Calculating Voronoi Tesselation..." )]
                [vor (Voronoi nodes)]
                [a (print "Took..." (- (time) start-time) "seconds")]
                [a (print (len vor.ridge_vertices) "ridges in total\n")]

                ;; filter segments and merge when there is no intersections
                [a (print "Filter and merge segments..." )]
                [filtered-edges (filter-edges2 vor cur-paths)]
                [ind-edges (if (:simplify config)
                             (simplify-graph filtered-edges vor cur-paths)
                             filtered-edges)]

                [a (print "Took..." (- (time) start-time) "seconds")]
                [a (print (len filtered-edges) "edges found")]
                [a (print (len ind-edges) "edges left after simplifying\n")]

                ;; indices to coordinates
                [edges (list (replace2d vor.vertices ind-edges))]

                [(, holes cedges) (get-holes ind-edges vor
                                             (:max-hole-connectors config)
                                             (:max-hole-length config))]
                [repulsed (repulse-holes holes)]]
            (for [point nodes] (add-circle root {"fill" "yellow" "r" "0.2"} point))
            (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} line))
            (for [pol (replace vor.vertices (if (:repulse-edges config) repulsed holes))]
              (add-polygon root {"fill" "none" "stroke" "cyan" "stroke-width" "0.2"}
                           (list (butlast pol))))
            (when (:hole-indications config)
              (for [line (replace2d vor.vertices cedges)]
                (add-line root {"fill" "none" "stroke" "blue" "stroke-width" "0.2"} line))))))
      (.write tree outfile)
      (print "\n TOTAL TIME IN SECONDS:\n" (- (time) start-time))))
