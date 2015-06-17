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
        [itertools :as it])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (reversed (list (range 2 3 ))))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(def config (dict [[:min-path-node-number 0] ;; ignore paths with few nodes
                   [:min-group-length 4] ;; minimum for voronoi to work
                   [:simplify True]
                   [:skeletons-only False]
                   [:max-hole-length 50000.0] ;500 or 1000, critical
                   [:max-hole-connectors 8]
                   [:hole-indications False]
                   [:repulse-edges True]
                   [:hole-choser :dist-angle]
                   [:hole-choser-max-d 3]
                   [:hole-choser-min-angle 130]
                   [:hole-choser-area-error 0.25]
                   [:pt-weight-threshold 2]]))


;;main
(for [infile testfiles]
  (do (print "\nProcessing" infile "...\n")
      (setv start-time (time))
      (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-all-wedges.svg"))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))

      (print "Calculating path matrices..." )
      (setv paths-list (get-paths2 root namespace))
      (print "Took..." (- (time) start-time) "seconds\n")

      (for [paths paths-list];; for each svg path node ...
        ;; like inkscape: 1 path= 1 object
        ;; TODO: allow different paradigm: automatic grouping across path objects
        (setv mpl-paths (list (map (fn [p] (Path (list (map first p)))) ;; ggf take mpl path before preparation -> speedup
                                   (filter (fn [p] (> (len p) (:min-path-node-number config)))
                                           (map prepare paths)))))
        (setv (, n labels) (divide-problem mpl-paths))

        (for [current-group (range n)]
          (setv cur-paths (list (replace mpl-paths (first (where (= current-group
                                                                    (array labels)))))))
          (setv nodes (vstack (list (map (fn [x] x.vertices) cur-paths))))
          (if (< (len nodes) (:min-group-length config))
            (do (print "Part" current-group "of" infile "skipped due to insufficient number of nodes")
                None
             )
            ;; causes weird behaviour of for-loop?
            (let [[a (print "Part" current-group "of" infile)]
                  [a (print "Calculating Voronoi Tesselation..." )]
                  [vor (Voronoi nodes)]
                  [a (print "Took..." (- (time) start-time) "seconds")]
                  [a (print (len vor.ridge_vertices) "ridges in total\n")]

                  ;; filter segments and merge when there is no intersections
                  [a (print "Filter and merge segments..." )]
                  [(, filtered-edges ridge-points) (filter-edges3 vor cur-paths)]
                  [ind-edges (if (:simplify config)
                               (simplify-graph filtered-edges vor cur-paths)
                               filtered-edges)]

                 [a (print "Took..." (- (time) start-time) "seconds")]
                 [a (print (len filtered-edges) "edges found")]
                 [a (print (len ind-edges) "edges left after simplifying\n")]

                  ;; indices to coordinates
                  [edges (list (replace2d vor.vertices ind-edges))]
                  [edge-map (get-edge-map ind-edges)]

                  [(, repulsed cedges) (get-holes ind-edges vor
                                                  (:max-hole-connectors config)
                                                  (:max-hole-length config)
                                                  True)]
                  [(, chosen whole) (get-wedge-holes2 repulsed vor config)]

                  [centers1 (list (map (fn [hole] (np.mean (list (replace vor.vertices hole)) 0))
                                       chosen))]
                  [skeletons1 (list (map (fn [hole center whole-hole]
                                           (list (wedge-skeleton hole edge-map vor cur-paths
                                                                 center whole-hole)))
                                         chosen centers1 whole))]

                  ;;    [(, holes cedges) (get-holes ind-edges vor (:max-hole-connectors config) (:max-hole-length config))]
                  ;; [repulsed (list (map (fn [x] (list (butlast x))) (repulse-holes holes)))]
                  ;; [chosen (get-wedge-holes repulsed vor config)]
                  ;;  [centers1 (list (map (fn [hole] (np.mean (list (replace vor.vertices hole)) 0)) chosen))]
                  ;; [wedge-skeletons1 (list (map (fn [hole center] (list (wedge-skeleton hole edge-map vor cur-paths center))) chosen centers1))]

                  ;;
                  ;; new, todo: use whole
                  ;;

                  ;;[(, skeletons2 centers2) (find-wedges ind-edges vor filtered-edges ridge-points wedge-skeletons1 cur-paths)]
                  [(, skeletons2 centers2) (find-wedges2 ind-edges vor filtered-edges ridge-points skeletons1 whole cur-paths (:pt-weight-threshold config))]

                  ]
          ;;    (for [point nodes] (add-circle root {"fill" "yellow" "r" "0.2"} point))
          ;;    (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} line))
              (if (:skeletons-only config)
                (do (for [pol (replace2d vor.vertices chosen)]
                      (add-polygon root {"fill" "none" "stroke" "blue" "stroke-width" "0.2"} pol))
                    (for [sk skeletons1]
                      (for [limb sk]
                        (add-line root {"fill" "none" "stroke" "blue" "stroke-width" "0.2"}
                                  [(get vor.vertices (first limb)) (get vor.vertices (last limb))])))
                    (for [(, sk c) (zip skeletons2 centers2)]
                      (for [limb sk]
                        (add-line root {"fill" "none" "stroke" "cyan" "stroke-width" "0.2"}
                                  [c (get vor.vertices (last limb))]))))
                (do (for [wedge (list (map (fn [sk center] (wedge-from-skeleton sk vor center :double))
                                           skeletons1 centers1))]
                      (add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.3"} wedge))
                    (for [wedge (list (map (fn [sk center] (wedge-from-skeleton sk vor center :double))
                                           skeletons2 centers2))]
                      (add-path root {"fill" "none" "stroke" "cyan" "stroke-width" "0.3"} wedge))))
              None))))
      (.write tree outfile)
      (print "\n TOTAL TIME IN SECONDS:\n" (- (time) start-time))))
