(import [xml.etree.ElementTree [parse SubElement]]
;        [misc [*]]
        [var [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [graph [*]]
        [tools [*]]
        [numpy [*]]
        [numpy :as np]
        [time [time]]
        [operator [add sub]]
        [scipy.spatial [Voronoi]]
        [matplotlib.path [Path]]
        [itertools [chain]]
        [itertools :as it]
        )

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (reversed (list (range 8 9))))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(def config (dict [[:min-path-node-number 0] ;; ignore paths with few nodes
                   [:simplify True]
                   [:max-hole-length 50000.0] ;500 or 1000, critical
                   [:max-hole-connectors 8]
                   [:hole-indications False]
                   [:repulse-edges True]
                   [:hole-choser :dist-angle]
                   [:hole-choser-max-d 3]
                   [:hole-choser-min-angle 130]
                   [:hole-choser-area-error 0.6]]))


;;main
(for [infile testfiles]
  (do (print "\nProcessing" infile "...\n")
      (setv start-time (time))
      (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-chosen-holes.svg"))))
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

                [a (print "Part" current-group "of" infile "\n")]
                [a (print "Calculating Voronoi Tesselation..." )]
                [time1 (time)]
                [vor (Voronoi nodes)]
                [a (print "Took..." (- (time) time1) "seconds")]
                [a (print (len vor.ridge_vertices) "ridges in total\n")]

                ;; filter segments and merge when there is no intersections
                [a (print "Filter and merge segments..." )]
                [time2 (time)]
                [(, filtered-edges ridge-points) (filter-edges3 vor cur-paths)]
                [ind-edges (if (:simplify config)
                             (simplify-graph filtered-edges vor cur-paths)
                             filtered-edges)]

                [a (print "Took..." (- (time) time2) "seconds")]
                [a (print (len filtered-edges) "inner edges found")]
                [a (print (len ind-edges) "edges after simplifying\n")]

                ;; indices to coordinates
                [edges (list (replace2d vor.vertices ind-edges))]
                [(, repulsed cedges) (get-holes ind-edges vor
                                                (:max-hole-connectors config)
                                                (:max-hole-length config)
                                                True)]
                [(, chosen whole) (get-wedge-holes2 repulsed vor config)]
                [a (print chosen)]
                ]
            ;;(for [point nodes] (add-circle root {"fill" "yellow" "r" "0.2"} point))
           ;; (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} line))
            (for [pol (replace2d vor.vertices (if (:repulse-edges config)
                                              repulsed
                                              (list (map butlast holes))))]
              (add-polygon root {"fill" "none" "stroke" "black" "stroke-width" "0.1"} pol))
            (for [pol (replace2d vor.vertices chosen)]
              (add-polygon root {"fill" "none" "stroke" "black" "stroke-width" "0.2"} pol)))))
      (.write tree outfile)
      (print "\n TOTAL TIME IN SECONDS:\n" (- (time) start-time))))
