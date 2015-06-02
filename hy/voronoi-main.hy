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
        [scipy.spatial [Voronoi Delaunay]]
        [scipy.sparse.csgraph [connected_components depth_first_tree shortest_path]]
        [scipy.sparse [csr_matrix]]
        [matplotlib.path [Path]]
        [itertools [chain]]
        [itertools :as it]
        )

(def testfiles
;;  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (range 2 3)))
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"])) (reversed (list (range 5 7)))))
  )

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))


(def config (dict [[:min-path-node-number 0] ;; ignore paths with few nodes
                   [:min-angle 80]
                   [:max-angle 150]
                   [:max-cycle-length 50]
                   [:max-cycle-nbh-length 1]
                   [:min-edge-ratio 0.05] ;; percentage min edge has to have
                    ]))



(defn get-cycles [edges N];; should start at leave!
  "Returns set of nodes that are contained in cycles"
  (let [[graph (csr_matrix (, (ones (len edges))
                              (transpose edges))
                           (, N N))]
        [dfs-tree (depth_first_tree graph (first (first edges)) false)]
        [cycle-edges (.difference (make-set edges)
                                  (make-set (transpose (.nonzero dfs-tree))))]
        ]
    (print "ce" cycle-edges)

    (list (apply concat cycle-edges))))

;; (get-cycles [[0 1] [1 2] [2 3] [3 1] [4 3]] 5)

;; (= [3 4] (get-cycles [[0 1] [1 2] [1 3] [2 3] [3 0] [4 2]] 5))
;;(= [3] (get-cycles [[0 1] [1 2] [1 3] [2 3] [3 0]] 5))


(defn find-cycles-rec [node dest dist-matrix paths]
  (if (= dest node)
    paths
    (for [ind (second (nonzero (get dist-matrix node)))]

      )))

(defn find-cycles [edge dist-matrix]
  (setv node (first edge))
  (setv paths [])

  )

(defn get-path [pred-matrix nodes]
  (setv node (first nodes))
  (setv dest (second nodes))
  (setv pred (get pred-matrix node))
  (setv path [node pred])
  (if (pos? pred)
    (do
     (while (not (= pred dest))
       (setv pred (get pred-matrix pred))
       (.append path pred))
     path)
    None))

(defn get-holes [cycle-edges dist-matrix] ;as csr_matrix
  (setv paths [])
  (for [edge cycle-edges]
    (setv d-matrix dist-matrix)
    (assoc d-matrix edge 0)
    (setv (, sp predecessors) (shortest_path d-matrix "D" False True))
    (when (< (get sp edge) (:max-cycle-length config))
      (.append paths (get-path predecessors edge)))))


(defn get-edge-length [vor inds] ;test
  (linalg.norm (apply sub (replace vor.vertices inds)) 2))

(defn reconstruct-wedges [wedges] ;; wedges being an object
  (array (list (map (fn [wedge] (list (map (fn [x y z] (elevate-bezier [x y z]))
                                           (:vertices wedge)
                                           (list (repeat (:center wedge) 3))
                                           (roll (:vertices wedge) -1 0))))

                    wedges))))

(defn dls [edges N start-node limit] ;; test : newest extension !!!!
  (let [[graph (csr_matrix (, (ones (len edges))
                              (transpose edges))
                           (, N N))]
        [pred (* (ones (, N)) inf)]] ; empty matrix
    (dls-rec graph pred start-node start-node limit)))

(defn trace-back [pred start]
  (setv node start)
  (setv path [start])
  (while (isfinite (get pred node))
    (.append path (get pred node)))
  path)

(defn dls-rec [graph pred start-node node limit]
  (when (> limit 0)
    (for [v (.nonzero (second (.getrow graph node)))]
      (if (= start-node v)
        (do
         (assoc pred v node)
         (print (trace-back pred v)))
        (when (isinf (get pred v))
          (do (assoc pred v node)
              (dls-rec graph pred start-node v (dec limit))))))))

;;main
(for [infile testfiles] ;; test
  (do (print "\nProcessing" infile "...\n")
      (setv start-time (time))
      (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-voronoi-all.svg"))))
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
          (let [[cur-paths (list (replace mpl-paths (first (where (= current-group (array labels))))))]
                [nodes (vstack (list (map (fn [x] x.vertices) cur-paths)))]

                [a (print "Calculating Voronoi Tesselation..." )]
                [vor (Voronoi nodes)]
                [a (print "Took..." (- (time) start-time) "seconds\n")]

                ;; filter segments and merge when there is no intersections
                [a (print "Filter and merge segments..." )]
                [filtered-edges (filter-edges2 vor cur-paths)]
                [ind-edges (simplify-graph filtered-edges)]

                [a (print "Took..." (- (time) start-time) "seconds\n")]
                [a (print (len filtered-edges) "edges found")]
                [a (print (len ind-edges) "edges left after simplifying")]

                ;; indices to coordinates
                [edges (replace2d vor.vertices ind-edges)]


                [emap (get-edge-map ind-edges)]

                [leaves (get-graph-leaves emap)]


                [mwedges (list (filter (fn [entry] (> (len (second entry)) 2))
                                       (.iteritems emap)))] ;; highly probable wedges

                [centers (list (map first mwedges))]
                [vertices (list (map (fn [x] (list (second x))) mwedges))]

                [wedge-ridges (list (map (fn [w] (list (map (fn [v] [(first w) v])
                                                            (second w))))
                                         mwedges))]

                [nleaves (list (map (fn [w] (len (.intersection leaves (second w))))
                                    mwedges))]

                [wedge-angle (list (map (fn [w]
                                          (list (map (fn [pair] (int (rad2deg (arccos (apply dot pair)))))
                                                     (it.combinations
                                                      (list
                                                       (map (fn [v] (normalize
                                                                     (- (get vor.vertices (first w))
                                                                        (get vor.vertices v))))
                                                            (second w)))
                                                      2))))
                                        mwedges))]

                [vertex-pairs (list (map (fn [vs] (list (it.combinations vs 2))) vertices))]

                [wedge-angle2
                 (list (map (fn [vpairs c]
                              (list (map (fn [pair]
                                           (int (rad2deg (arccos (apply dot (list
                                                                             (map (fn [v] (normalize (- v c)) )
                                                                                  (replace vor.vertices pair))))))))
                                         vpairs)))
                            vertex-pairs
                            (replace vor.vertices centers)))]
                [outer-edge? (list (map (fn [vs] (list (map (fn [v] (in v leaves)) vs))) vertices))]
                [mergeable?
                 (list (map (fn [vpairs]
                              (list (map (fn [pair] (> 2 (len (.intersection leaves pair)))) vpairs)))
                            vertex-pairs))]

                                ;   [a (print " amax?" wedge-angle wedge-angle2)]

                [lengths (list (map (fn [c vs] (list (map (fn [v] (get-edge-length vor [c v])) vs)))
                                    centers vertices))]

                ;; [dist-matrix (csr_matrix (, lengths (transpose edges)) (, N N))] ;; not symmetric
                [XY (transpose edges)]
              [a (print (shape XY) )]
              [a (print (shape edges) )]
              ;;  [dist-matrix (csr_matrix (, (list (repeat lengths 2)) (hstack (, XY (roll XY 1 0)))) (, N N))] ;; solve for empty "edges"

                [a (print "dls:")]
                [c (dls edges (len nodes) (first leaves) 8)]

                [min-angles (if (empty? wedge-angle) [] (amin wedge-angle 1))]
                [max-angles (if (empty? wedge-angle) [] (amax wedge-angle 1))]


              [a (print "centers" centers )]
              [a (print "vertices" vertices)]
              [a (print "nleaves" nleaves)]
              [a (print "vertex-pairs" vertex-pairs)]
              [a (print "angles " wedge-angle2)]
              [a (print "mergeable?" mergeable?)]
              [a (print "outer-edge?" outer-edge?)]
              [a (print "leaves " leaves)]
              [a (print "min-angles " min-angles)]
              [a (print "max-angles " max-angles)]

                [cycles (get-cycles ind-edges (len vor.vertices))]
            ;    [a (print (.intersection leaves (set  cycles)))]
                [wedge-nr (range (len centers))]

                [wedges (map (fn [nr c vs nl mina maxa l]
                               (dict [[:id nr] [:ind-center c] [:ind-vertices vs]
                                      [:center (get vor.vertices c)] [:vertices (get vor.vertices vs)]
                                      [:nleaves nl] [:min-angle mina] [:max-angle maxa]
                                      [:lengths l]
                                      [:leaves (list (map (fn [v] (in v leaves)) vs))]]))
                          wedge-nr centers vertices nleaves min-angles max-angles lengths)]

                [wedges-sure (list (filter (fn [x] (and (>= (:nleaves x) 2)
                                                        (>= (:min-angle x) (:min-angle config))
                                                        (<= (:max-angle x) (:max-angle config))
                                                         ;; small outer edges not with wedges
                                                        (all (logical_or (= False (:leaves x))
                                                                         (< (:min-edge-ratio config)
                                                                            (normalize (:lengths x)))))))
                                           wedges)) ]
                [a (print "sure-wedges" wedges-sure)]
                [explained-vertices ()]
                [wedges-rec (reconstruct-wedges wedges-sure)]
                [a (print "rec" wedges-rec)]


                ]

            (for [point nodes] (add-circle root {"fill" "yellow" "r" "0.2"} point))
            (for [line edges] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.1"}
                                        line))
            (for [point (replace vor.vertices cycles)]
              (add-circle root {"fill" "cyan" "r" "0.2"} point))
            (for [wedge wedges-rec]
              (add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.5"} wedge))
            ))


        )
      (.write tree outfile)
      (print "\n TOTAL TIME IN SECONDS:\n" (- (time) start-time))))
