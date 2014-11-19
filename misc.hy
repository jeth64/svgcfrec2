(import [xml.etree.ElementTree [parse SubElement]]
        [svg [*]]
        [bezier [*]]
        [modulo [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [itertools [groupby combinations]]
        [itertools :as it]
        [operator [sub mul add]]
        [scipy.spatial [ConvexHull]]
        [scipy.spatial.distance [pdist squareform]]
        [scipy.sparse.csgraph [connected_components]]
     ;   [matplotlib.pyplot [*]]
        )


(defn get-edge-maps [edges]
  (let [[vertices (transpose (array (list edges)))]]
    [(group-by (first vertices) (second vertices))
     (group-by (second vertices) (first vertices))]))

(defn straight-v? [p33 p66]
  "Uses the 33th and 66th percentile of a list of direction vectors of nodes to determine if segment is approximately straight"
  (> 0.5 (first (- p33 p66))))

(defn straight? [p33 p66]
  "Uses the 33th and 66th percentile of a list of angles to determine if segment is approximately straight"
  (> 20 (mod-dist 180 p33 p66)))

(defn get-directions-v [connected-comp node-dirs nodes]
  (remove none? (map
                 (fn [l] (let [[dirs (list (map (fn [i] (get node-dirs i)) l))]
                               [pts (list (map (fn [i] (get nodes i)) l))]
                               [p33 (percentile dirs 33 0)]
                               [p66 (percentile dirs 66 0)]
                            ;   [b (print "\n\npts: " pts)]
                             ;  [a (print "\ndirs: " (array dirs))]
                              ; [c (print "\np33: " p33 "  p66: " p66)]
                             ;  [x (print "\nstraight?  " (straight? p33 p66))]
                             ;  [d (print "\nbox: " (bounding-box pts (mean [p33 p66] 0)))]
                               ]
                           (if (straight? p33 p66)
                            (mean (partition 2 2 (bounding-box pts (mean [p33 p66] 0))) 1) ;changed
                             None)))
                 connected-comp)))

(defn get-directions [connected-comp node-dirs nodes] ;; returns line
  (remove none? (map
                 (fn [l] (let [[dirs (list (map (fn [i] (get node-dirs i)) l))]
                               [pts (list (map (fn [i] (get nodes i)) l))]
                               [p33 (percentile dirs 33)]
                               [p66 (percentile dirs 66)]
                             ;  [b (print "\n\npts: " (array pts))]
                             ;  [a (print "\ndirs: " dirs)]
                             ;  [c (print "\np33: " p33 "  p66: " p66)]
                             ;  [x (print "\nstraight?  " (straight? p33 p66))]
                             ;  [d (print "\nbox: " (bounding-box pts (mean [p33 p66])))]
                               ]
                           (if (straight? p33 p66)
                            (mean (partition 2 2 (bounding-box pts (mean [p33 p66]))) 1)
                             None)))
                 connected-comp)))

(defn get-directions2 [connected-comp node-dirs nodes] ;; zum testen
  (list (map
         (fn [l] (let [[dirs (list (map (fn [i] (get node-dirs i)) l))]
                       [pts (list (map (fn [i] (get nodes i)) l))]
                       [p33 (percentile dirs 33)]
                       [p66 (percentile dirs 66)]
                                ;   [b (print "\n\npts: " (array pts))]
                                ;  [a (print "\ndirs: " dirs)]
                                ;  [c (print "\np33: " p33 "  p66: " p66)]
                                ;   [x (print "\nstraight?  " (straight? p33 p66))]
                                ;   [d (print "\nbox: " (bounding-box pts (mean [p33 p66])))]
                       ]
                   [(array pts) dirs (bounding-box pts (mean [p33 p66]))]))
         connected-comp)))


(defn ext-normalize [v]
  "Extends vector normalization s.t. first element in returned vector is always positive"
  (if (pos? (first v))
    (/ v (linalg.norm v))
    (/ (negative v) (linalg.norm v))))

(defn normalize [v]
  "Simple vector normalization"
  (/ v (linalg.norm v)))

(defn get-node-dirs-old [path]
  "Get curve gradients of nodes in path"
  (let [[node-dirs (list (map (fn [segment] (ext-normalize (apply sub (slice segment 0 None 3))))
                              path))]
        [b (print (array node-dirs))]
        [a (print (mean [node-dirs (roll node-dirs 1)] 0))]
        ]
    (mean [node-dirs (roll node-dirs 1)] 0)))

(defn angle [v]
  "Angle of v in degrees as in polar coordinate system"
  (% (degrees (apply arctan2 (list (reversed v)))) 180))

(defn get-node-dirs [path]
  "Get curve gradients of nodes in path, as angle to x-axis"
  (let [[node-dirs (list (map (fn [segment] (angle (apply sub (slice segment 0 None 3))))
                              path))]]
    (list (map (fn [x y] (mod-mean 180 x y)) node-dirs (roll node-dirs 1)))))

(defn find-thinnings [path max-dist min-neigh-dist]
  (let [
      ;  [a (print "path2" path2)]
        [nodes (list (map first path))]
      ;  [node-dirs-old (get-node-dirs-v path)]
        [node-dirs (get-node-dirs path)]

        [n (len nodes)]
        [d (squareform (pdist (array nodes)))]

        [inds (- (tile (range n) [n 1])
                 (transpose (tile (range n) [n 1])))]
       ; [edges (make-set (filter (fn [x] (< 1 (apply sub x) (dec n)))
        ;                         (transpose (where (< d max-dist)))))]

                                ;   [a (print edges "\n")]
     ;   [ind-dists (list (map (fn [x] (mod-dist n (first x) (second x))) edges))]
                                ;  [d-edge-map (group-by ind-dists edges)]

        [graph (& (< d max-dist) (< min-neigh-dist inds) (< inds (- n min-neigh-dist)))]
        [cc (connected_components graph False)]
        [cc2 (list (filter (fn [y] (> (len y) 1))
                           (.values (group-by (second cc) (range 50))))) ]
     ;   [b (print "nodes" nodes)]
     ;   [b (print "node-dirs" node-dirs)]
     ;   [a (print "cc2" cc2)]
        [dirs (get-directions cc2 node-dirs nodes)]]
    dirs
    ))

(defn find-thinnings2 [path max-dist min-neigh-dist];; zum testen
  (let [[nodes (list (map first path))]
        [node-dirs (get-node-dirs path)]

        [n (len nodes)]
        [d (squareform (pdist (array nodes)))]

        [inds (- (tile (range n) [n 1])
                 (transpose (tile (range n) [n 1])))] ;schwachsinn

        [graph (& (< d max-dist) (< min-neigh-dist inds) (< inds (- n min-neigh-dist)))]
        [cc (connected_components graph False)]
        [cc2 (list (filter (fn [y] (> (len y) 1))
                           (.values (group-by (second cc) (range 50))))) ]
        [dirs (get-directions2 cc2 node-dirs nodes)]]
    dirs
    ))




(defn get-wedge [index-pairs inds isecs vertices]
  (let [[gisecs (list (map (fn [x] (get isecs x)) index-pairs))]
        [vs (list (map (fn [x] (get vertices x)) inds))]]
    (if (points-in-triangle? vs gisecs)
      [(array (list (map (fn [start end] (elevate-degree start (mean gisecs 0) end))
                         vs (roll vs 1))))
       (mean gisecs 0)]
      None)))


(defn make-set [pairings]
  (set (list (map (fn [x] (tuple (sort x))) pairings))))


(defn check-thinning [edge distance-edge-map n-nodes]
  (let [[left [(first edge)]]
        [right [(second edge)]]]
    (for [d (range 3 (inc (/ n-nodes 2)))]
      (get distance-edge-map d))))


;; min-node dist is minimum distance between neighbouring nodes
(defn find-thinnings-new [path max-dist];test
  (let [
        [nodes (list (map first path))]
       ; [a (print "n" (array nodes))]
        [node-dirs (get-node-dirs path)]

    ;    [a (print "n-d" node-dirs)]
        [n (len nodes)]

        [node-dists (pdist (array nodes))]
        [d (squareform node-dists)]

        [min-node-dist (min node-dists)]
    ;    [a (print "min-node-dist" min-node-dist)]
     ;   [a (print "max-node-dist" (max node-dists))]

        [mid-el (int (ceil (/ n 2)))]
        [half-range (range 1 mid-el)]
        [first-row (list (if (even? n)
                           (concat [0] half-range [mid-el] (reversed half-range))
                           (concat [0] half-range (reversed half-range))))]
        [ind-dist (array (list (map (fn [x] (roll (list first-row) x))
                                    (range n))))]
     ;   [a (print ind-dist)]
   ;     [a (print "ishape" (shape ind-dist) "d" (shape d))]

        [trace-dist (* ind-dist min-node-dist)]
        [limit (choose (> trace-dist max-dist) [trace-dist max-dist])]
        [graph (< d limit) ]

;        [b (print "graph" graph)]
        [cc (connected_components graph False)]
        [cc2 (list (filter (fn [y] (> (len y) 1))
                           (.values (group-by (second cc) (range 50))))) ]

        [dirs (get-directions cc2 node-dirs nodes)]]
    dirs
    ))
