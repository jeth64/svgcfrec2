 (import [xml.etree.ElementTree [parse SubElement]]
        [svg [*]]
        [bezier [*]]
        [modulo [*]]
        [geometry [*]]
        [graph [*]]
        [tools [*]]
        [numpy [array argmin transpose sum size negative arctan2 degrees deg2rad dot diag roll hstack shape cross]]
        [numpy.linalg [norm]]
        [numpy.random :as ran]
        [scipy.sparse.csgraph [connected_components]]
        [operator [sub mul add]]
        [itertools [groupby combinations]]

        [itertools :as it]
        [scipy.spatial.distance [pdist cdist squareform]])
;; new misc.hy...

;;
;; Helper functions
;;

(defn angle [v] ;; clash with numpy.angle
  "Angle of v in degrees as in polar coordinate system"
  (assert (= (size v) 2))
  (let [[radian (apply arctan2 (reversed (.tolist v)))]]
    (assert (not (coll? radian)) )
    (% (+ (degrees radian) 90) 180)))


(defn dirvec [polar-theta]
  "Inverse function of 'angle'"
  (dot (rotation-matrix (deg2rad (% (- polar-theta 90) 180)))
       [1 0]))

(defn ext-normalize [v]
  "Extends vector normalization s.t. second element in returned vector is always positive"
  (if (neg? (second (flatten v)))
    (/ (negative v) (norm v))
    (/ v (norm v))))

(defn normalize [v]
  "Simple vector normalization"
  (/ v (norm v)))


;;
;; Other functions
;;
(defn get-node-dirs [path &optional [mode :gradient] [direction :bi] [angle-range [0]]] ;test
  "Get polar angle of curve gradients of nodes in path with a given range of angles around result"
  ;; modes: :node2node - direction vectors between consecutive nodes is used
  ;;        :gradient - curve gradients
  ;; directions: :forward :backward :bi
  (let [[node-dirs1 (when (not (= direction :backward))
                      (list (map (fn [segment]
                                   (let [[theta (angle (apply sub (if (= mode :node2node)
                                                                    (slice segment 0 None 3)
                                                                    (slice segment 0 2))))]]
                                     (list (map (fn [x] (% (+ theta x) 180)) angle-range))))
                                 path)))]
        [node-dirs2 (when (not (= direction :forward))
                      (roll (if (= mode :node2node)
                              node-dirs1
                              (list (map (fn [segment]
                                           (let [[theta (angle (apply sub (if (= mode :node2node)
                                                                            (slice segment 0 None 3)
                                                                            (slice segment 2 None))))]]
                                             (list (map (fn [x] (% (+ theta x) 180)) angle-range))))
                                         path)))
                            1))]]
    (if (= :bi direction)
      (hstack [node-dirs1 node-dirs2])
      (case direction
            :forward node-dirs1
            :backward node-dirs2))))

(defn get-node-dirs-v [path &optional [mode :gradient] [direction :bi]] ;test
  "Get curve gradients of nodes  in path as vectors"
  ;; modes: :node2node - direction vectors between consecutive nodes is used
  ;;        :gradient - curve gradients
  ;; directions: :forward :backward :bi
  (let [[node-dirs1 (when (not (= direction :backward))
                      (list (map (fn [segment]
                                   (apply sub (if (= mode :node2node)
                                                  (slice segment 0 None 3)
                                                  (slice segment 0 2))))
                                 path)))]
        [node-dirs2 (when (not (= direction :forward))
                      (roll (if (= mode :node2node)
                              node-dirs1
                              (list (map (fn [segment]
                                           (apply sub (if (= mode :node2node)
                                                        (slice segment 0 None 3)
                                                        (slice segment 2 None))))
                                         path)))
                            1))]]
    (if (= :bi direction)
      (hstack [node-dirs1 node-dirs2])
      (case direction
            :forward node-dirs1
            :backward node-dirs2))))


(defn get-curvature [nodes node-dirs e] ;test - reaction erratic
  (cross (get node-dirs (first e))
         (apply sub (replace nodes e)))) ;;ext-normalize?

(defn inner-edge? [nodes node-dirs e] ;;dirty hack: curvature returns array instead of scalar value
  (let [[c (get-curvature nodes node-dirs e)]]
    (pos? (first (array [c])))))

(defn get-tiling-bins [ps nTiles overlap]
  (let [[(, minX minY) (amin ps 0)]
        [(, maxX maxY) (amax ps 0)]
        [(, nTilesX nTilesY) nTiles]
        [(, overlapX overlapY) overlap]
        [sizeX (/ (+ (- maxX minX) (* (dec nTilesX) overlapX))
                  nTilesX)]
        [sizeY (/ (+ (- maxY minY) (* (dec nTilesY) overlapY))
                  nTilesY)]
        [binsX (list (rest (reductions add
                                       (cons (- sizeX overlapX)
                                             (map (fn [i] (if (= (% i 2) 0) overlapX
                                                              (- sizeX (* 2 overlapX))))
                                                  (range (- (* nTilesX 2) 3))))
                                       minX)))]
        [binsY (list (rest (reductions add
                                       (cons (- sizeY overlapY)
                                             (map (fn [i] (if (= (% i 2) 0) overlapY
                                                              (- sizeY (* 2 overlapY))))
                                                  (range (- (* nTilesY 2) 3))))
                                       minY)))]]
    (, binsX binsY)))

(defn tiling [ps bins]
  "Every second bin is overlapping part of consecutive bins"
  (let [[(, binsX binsY) bins]
        [xy (transpose ps)]
        [indsX (digitize (first xy) binsX)]
        [indsY (digitize (second xy) binsY)]]
    (, indX indY)))

;; voronoi:

(defn filter-edges [vor paths];; vor is voronoi object ;; check numpy stuff
  ;; nur invalid edges verwenden?
  (let [[acc-lengths (list (reductions add (map (fn [p] (len p.vertices)) mpl-paths)))]

        ;; vertices that are not within the compound path (in hole or outside any path)
        ;; (like raycasting: counts number of paths it is in; entire paths are needed)
        [invalid-vertices (set (map first
                                    (filter
                                     (fn [y] (even? (sum (second y))))
                                     (zip (range (len vor.vertices))
                                          (transpose
                                           (array
                                            (list (map (fn [p] (.contains_points p vor.vertices))
                                                       paths))))))))]

         ;; edge not crossing path with orthogonal edge having sequentiel indices
        [inter-path-edges (if (= 1 (len paths))
                            (set [])
                            (set (zipwith (fn [x y] (tuple (sort [x (dec y)])))
                                          (cons 0 (butlast acc-lengths))
                                          (roll acc-lengths 1))))]

        ;; edge crossing path with orthogonal edge not having sequentiel indices
        [invalid-edges (set (map (fn [x y] (tuple (sort [x (dec y)])))
                                 (cons 0 (butlast acc-lengths))
                                 acc-lengths))]]
    (list (map second (filter (fn [e] (and
                                       ;; no infinite edges
                                       (not (any (= -1 (asarray (second e))))) ;; asarray nötig?
                                       ;; no edges crossing the path
                                       (or (and (> (abs (apply sub (first e))) 1)
                                                (not-in (tuple (sort (first e))) invalid-edges))
                                           (in (tuple (sort (first e))) inter-path-edges))
                                       ;; no edges with both points outside path
                                       (> 2 (len (.intersection invalid-vertices (second e))))))
                              (zip (.tolist vor.ridge_points)
                                   vor.ridge_vertices))))))

(defn filter-edges2 [vor paths] ;; check numpy stuff
  (let [[acc-lengths (list (reductions add (map (fn [p] (len (list p.vertices))) paths)))]

         ;; edge not crossing path with orthogonal edge having sequentiel indices
        [inter-path-edges (if (= 1 (len (list paths)))
                            (set [])
                            (set (zipwith (fn [x y] (tuple (sort [x (dec y)])))
                                          (cons 0 (butlast acc-lengths))
                                          (roll acc-lengths 1))))]

        ;; edge crossing path with orthogonal edge not having sequentiel indices
        [invalid-edges (set (map (fn [x y] (tuple (sort [x (dec y)])))
                                 (cons 0 (butlast acc-lengths))
                                 acc-lengths))]

        [new-edges (list (map second (filter (fn [e] (and
                                       ;; no infinite edges
                                       (not (any (= -1 (asarray (second e))))) ;; asarray nötig?
                                       ;; no edges crossing the path
                                       (or (and (> (abs (apply sub (first e))) 1)
                                                (not-in (tuple (sort (first e))) invalid-edges))
                                           (in (tuple (sort (first e))) inter-path-edges))))
                              (zip (.tolist vor.ridge_points)
                                   vor.ridge_vertices))))] ;; hier ansetzen fuer weights
        [G (zeros (len vor.vertices))]

        ;; vertices that are not within the compound path (in hole or outside any path)
        ;; (like raycasting: counts number of paths it is in; entire paths are needed)
        [invalid-vertices (set (map first
                                    (filter
                                     (fn [y] (even? (sum (second y))))
                                     (zip (range (len vor.vertices))
                                          (transpose
                                           (array
                                            (list (map (fn [p] (.contains_points p vor.vertices))
                                                       paths))))))))]]
    (list (map second (filter (fn [e] (and
                                       ;; no infinite edges
                                       (not (any (= -1 (asarray (second e))))) ;; asarray nötig?
                                       ;; no edges crossing the path
                                       (or (and (> (abs (apply sub (first e))) 1)
                                                (not-in (tuple (sort (first e))) invalid-edges))
                                           (in (tuple (sort (first e))) inter-path-edges))
                                       ;; no edges with both points outside path
                                       (> 2 (len (.intersection invalid-vertices (second e))))))
                              (zip (.tolist vor.ridge_points)
                                   vor.ridge_vertices))))))

(defn filter-edges3 [vor paths] ;; check numpy stuff
  (let [[acc-lengths (list (reductions add (map (fn [p] (len (list p.vertices))) paths)))]

         ;; edge not crossing path with orthogonal edge having sequentiel indices
        [inter-path-edges (if (= 1 (len (list paths)))
                            (set [])
                            (set (zipwith (fn [x y] (tuple (sort [x (dec y)])))
                                          (cons 0 (butlast acc-lengths))
                                          (roll acc-lengths 1))))]

        ;; edge crossing path with orthogonal edge not having sequentiel indices
        [invalid-edges (set (map (fn [x y] (tuple (sort [x (dec y)])))
                                 (cons 0 (butlast acc-lengths))
                                 acc-lengths))]

        [new-edges (list (map second (filter (fn [e] (and
                                       ;; no infinite edges
                                       (not (any (= -1 (asarray (second e))))) ;; asarray nötig?
                                       ;; no edges crossing the path
                                       (or (and (> (abs (apply sub (first e))) 1)
                                                (not-in (tuple (sort (first e))) invalid-edges))
                                           (in (tuple (sort (first e))) inter-path-edges))))
                              (zip (.tolist vor.ridge_points)
                                   vor.ridge_vertices))))] ;; hier ansetzen fuer weights
        [G (zeros (len vor.vertices))]

        ;; vertices that are not within the compound path (in hole or outside any path)
        ;; (like raycasting: counts number of paths it is in; entire paths are needed)
        [invalid-vertices (set (map first
                                    (filter
                                     (fn [y] (even? (sum (second y))))
                                     (zip (range (len vor.vertices))
                                          (transpose
                                           (array
                                            (list (map (fn [p] (.contains_points p vor.vertices))
                                                       paths))))))))]
        [filtered-edges (list (filter (fn [e] (and
                                               ;; no infinite edges
                                               (not (any (= -1 (asarray (second e))))) ;; asarray nötig?
                                               ;; no edges crossing the path
                                               (or (and (> (abs (apply sub (first e))) 1)
                                                        (not-in (tuple (sort (first e))) invalid-edges))
                                                   (in (tuple (sort (first e))) inter-path-edges))
                                               ;; no edges with both points outside path
                                               (> 2 (len (.intersection invalid-vertices (second e))))))
                                      (zip (.tolist vor.ridge_points)
                                           vor.ridge_vertices)))]
        ]
    (, (list (map second filtered-edges))
       (list (map first filtered-edges)))))

(defn divide-problem [mpl-paths]
  "Separate paths in smaller groups; paths are in the same group if one is included in the other "
  ;;(print "l" (len mpl-paths))
  (if (empty? mpl-paths)
    (, 0 [])
    (let [[p1-in-p2 (array (list (map (fn [p1] (list (map (fn [p2] (.contains_path p1 p2))
                                                          mpl-paths)))
                                      mpl-paths)))]
          [n (len mpl-paths)]]
      (setv cc (connected_components p1-in-p2 False))

      ;; calc-roots
      (setv depth (zeros (, n)))
      (setv depth2 (list (map (fn [ind] (get-tree-depth p1-in-p2
                                                        (first (where (= (second cc) ind))) ;;one group
                                                        depth)) ;; ind relative to view
                              (range (first cc)))))
      cc)))
