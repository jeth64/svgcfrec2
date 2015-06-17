(import [xml.etree.ElementTree [parse SubElement]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [numpy :as np]
        [operator [add sub]]
        [scipy.spatial [Voronoi]]
        [scipy.sparse.csgraph [connected_components depth_first_tree shortest_path]]
        [scipy.sparse [csr_matrix]]
        [matplotlib.path [Path]]
        [itertools :as it])

(defn get-edge-map [edges] ;test setitem to assoc
  (setv edge-map (dict []))
  (setv edge-set (make-set edges))
  (for [e edge-set]
    (if (in (first e) edge-map )
      (.add (get edge-map (first e)) (second e))
      (assoc edge-map (first e) (set [(second e)])))
    (if (in (second e) edge-map)
      (.add (get edge-map (second e)) (first e))
      (assoc edge-map (second e) (set [(first e)]))))
  edge-map)

(defn get-undir-adj-matrix [edges n-nodes &optional weights]
  "Returns adjacency matrix of undirected graph with given weights"
  (when (none? weights) (setv weights (np.ones (, (len edges) 1))))
  (csr_matrix (, (np.reshape (np.vstack [weights weights])
                             (, (* 2 (len edges))))
                 (if (empty? edges) (, [] [])
                   (np.transpose (np.vstack [edges (np.roll edges 1 1)]))))
              (, n-nodes n-nodes)))

(defn get-graph-leaves [edge-map]
  (set (map first (filter (fn [item] (= 1 (len (second item))))
                          (.iteritems edge-map)))))

(defn get-graph-isecs [edge-map]
  (set (map first (filter (fn [item] (> (len (second item)) 2))
                          (.iteritems edge-map))))            )


(defn simplify-graph-rec1 [node next-node edge-map not-2-conn visited new-edges vor paths]
  "produces double edges"
  (setv prev-node node)
  (while (not (in next-node not-2-conn))
    (setv candidate (first (.difference (get edge-map next-node) [prev-node])))
    (setv prev-node next-node)
    (setv line (Path (list (replace vor.vertices [node candidate])) [1 2]))
    (when (np.any (list (map (fn [p] (.intersects_path line p False)) paths)))
      (break))
    (setv next-node candidate))
  (.append new-edges [node next-node])
  (when (not (get visited next-node))
    (do
     (.__setitem__ visited next-node True)
     (for [n (.difference (get edge-map next-node) [prev-node])]
       (simplify-graph-rec1 next-node n edge-map not-2-conn visited new-edges vor paths))))
  new-edges)

(defn cond-divide-edge-rec [nodes vor paths]
  ;; check if connecting line crosses path
  (if (np.any (list (map (fn [p] (.intersects_path (Path [(get vor.vertices (first nodes))
                                                          (get vor.vertices (last nodes))]
                                                         [1 2])
                                                   p False))
                      paths)))
    (if (< (len nodes) 3)
      []
      (list (it.chain (cond-divide-edge-rec (list (take (np.ceil (/ (len nodes) 2)) nodes)) vor paths)
                      (cond-divide-edge-rec (list (drop (dec (np.ceil (/ (len nodes) 2))) nodes))
                                            vor paths))))
    [nodes]))

(defn simplify-graph-rec2 [node next-node edge-map not-2-conn visited new-edges vor paths]
  "produces too little edges in testfile 9"
  (setv prev-node node)
  (setv nodes [node next-node])
  ;; determine place of next intersection
  (while (not (in next-node not-2-conn))
    (setv candidate (first (.difference (get edge-map next-node) [prev-node])))
    (setv prev-node next-node)
    (setv next-node candidate)
    (.append nodes next-node))
  (.__setitem__ visited prev-node True)
  ;; divide edge if it intersect with a path
  (for [p (cond-divide-edge-rec nodes vor paths)]
    (.append new-edges [(first p) (last p)]))
  (when (not (get visited next-node))
    (do
     (.__setitem__ visited next-node True)
     (for [n (get edge-map next-node)]
       (when (not (get visited n))
         (simplify-graph-rec2 next-node n edge-map not-2-conn visited new-edges vor paths)))))
  new-edges)


(defn simplify-graph [edges vor paths]
  (if (> (len edges) 0)
    (let [[new-edges []]
          [edge-map (get-edge-map edges)]
          [leaves (get-graph-leaves edge-map)]
          [first-node (if (empty? leaves) (first (first edges)) (first leaves))]
          [next-node (first (get edge-map first-node))]
          [visited (.fromkeys (dict []) (.iterkeys edge-map) False)]
          [isecs (get-graph-isecs edge-map)]
          [not-2-conn (.union leaves isecs)]]
      (.__setitem__ visited first-node True)
      (simplify-graph-rec1 first-node next-node edge-map not-2-conn visited new-edges vor paths)
      ;;(array (list (make-set new-edges)))
new-edges
      ) ;; new: make-set
    edges))

(defn get-tree-depth-rec [adj-matrix node depths] ;; test
  (let [[neighbours (first (where (get adj-matrix node)))]]
     (if (empty? neighbours)
       depths
       (amax (array (list (map (fn [n] (get-tree-depth-rec adj-matrix
                                                           n
                                                           (do (assoc depths n (inc (get depths n)))
                                                               depths)))
                               neighbours))) 0))))

(defn get-tree-depth [adj-matrix idx depths] ;; test
  (get-tree-depth-rec adj-matrix
                      ;; get root:
                      (get idx (np.argmin (np.sum (np.array (list (replace (np.transpose adj-matrix)
                                                                           idx)))
                                                  1)))
                      depths))


(defn trace-back [pred start conn-limit]
  (setv node start)
  (setv path [start])
  (for [i (range (inc conn-limit))]
    (do (setv previous (get pred node))
        (if (np.isfinite previous)
          (do (.append path (int previous))
              (if (= previous start)
                (break)
                (setv node previous)))
          (do (print "Warning: Tree path tracing has reached dead end")
              (break))))
    (else (print "Warning: Tree path tracing has timed out!")))
  path)

(defn dls-rec [graph pred stop-node node pathlist cur-limit conn-limit length-limit]
  ;; side-effect: pathlist
  (when (> cur-limit 0)
    (for [v (second (.nonzero (.getrow graph node)))]
      (setv remaining-length (- length-limit (._get_single_element graph node v)))
      (when (not (or (= v (get pred node)) ;; avoid going directly back
                     (neg? remaining-length))) ;; do not consider node far along the path
        (if (= stop-node v)
          (do
           (assoc pred v node)
           (.append pathlist (trace-back pred v conn-limit)))
          (when (np.isinf (get pred v))
            (do (assoc pred v node)
                (dls-rec graph pred stop-node v pathlist
                         (dec cur-limit)
                         conn-limit
                         remaining-length)
                (assoc pred v inf))))))))

(defn get-leaves [graph]
  "Return graph leaves where graph is represented by sparse matrix"
  (first (np.where (= 1 (np.array (list (map (fn [i] (.getnnz (.getrow graph i)))
                                          (range (get (shape graph) 0)))))))))

(defn repulse-holes [holes]
  "Return list of unique holes"
  (list (map (fn [x] (first (second x)))
             (.iteritems (group-by (list (map (fn [x] (tuple (sorted (set x))))
                                              holes))
                                   holes)))))

(defn get-holes [edges vor conn-limit length-limit &optional [repulsed True]]
  "Returns lists of nodes describing paths around a hole"
  (if (empty? edges)
    (, [] [])
    (let [[N (len vor.vertices)]
          [dists (list (map (fn [e] (np.linalg.norm e)) (replace2d vor.vertices edges)))]
          [graph (get-undir-adj-matrix edges N dists)]
          [l (get-leaves graph)]
          [dfs-tree (depth_first_tree graph (first (first edges)) false)]
          [cycle-edges (.difference (make-set edges)
                                    (make-set (np.transpose (.nonzero dfs-tree))))]
          [paths []]]
      (list (map (fn [e] (dls-rec graph (* (ones (, N)) inf) (first e) (first e) paths
                                  conn-limit conn-limit length-limit))
                 cycle-edges))
      (setv res-paths (if repulsed (repulse-holes paths) paths))
      (, (list (map (fn [x] (list (butlast x))) res-paths)) cycle-edges))))

(defn fit-triangle [hole vor]
  (let [[entire-area (shoelace (list (replace vor.vertices hole)))]
        [tri-inds (list (it.combinations hole 3))]
        [tri-areas (list (map shoelace (replace2d vor.vertices tri-inds)))]
        [abs-errs (list (map (fn [A] (abs (- entire-area A))) tri-areas))]
        [min-ind (np.argmin abs-errs)]
        [rel-err (/ (get abs-errs min-ind) entire-area)]
        [best-triangle (list (get tri-inds min-ind))]]
    (dict [[:triangle best-triangle] [:rel-err rel-err] [:abs-err (get abs-errs min-ind)]]) ))

(defn get-wedge-holes-area [holes vor area-error] ;; test
  (if (empty? holes)
    (, [] [])
    (let [[entire-holes (copy holes)]
          [checked-holes (list (map (fn [hole]
                                      (if (> (len (list hole)) 3)
                                        (let [[best-fit (fit-triangle hole vor)]]
                                          (if (< (:rel-err best-fit) area-error)
                                            (:triangle best-fit)
                                            []))
                                        hole))
                                    entire-holes))]
          [filtered (list (remove (fn [hole] (empty? (first hole)))
                                  (zip checked-holes holes)))]]
      (, (list (map first filtered) (list (map second filtered)))))))

(defn angle-between [u v &optional otherwise]
  (let [[lu (np.linalg.norm u)]
        [lv (np.linalg.norm v)]
        [prod (np.dot (/ u lu) (/ v lv))]]
    (if (and (not (none? otherwise))
             (or (zero? lu) (zero? lv)))
      otherwise
      (np.degrees (np.arccos (cond [(< prod -1) -1]
                                   [(> prod 1) 1]
                                   [True prod]))))))

(defn closest-pt [hole &optional [mode :forward]]
  (let [[dists (list (map (fn [u v] (linalg.norm (- u v)))
                          hole (np.roll hole (if (= mode :forward) -1 1) 0)))]] ;; check
    (, (np.argmin dists) (np.amin dists))))

(defn widest-angle [hole]
  (let [[angles (list (map (fn [h1 h2 h3] (angle-between (- h2 h1) (- h2 h3)))
                           (np.roll hole 1 0) hole (np.roll hole -1 0)))]]
    (, (np.argmax angles) (np.amax angles))))

(defn remove-pts [hole vor max-d min-angle]
  (setv close-points True)
  (while (> (len hole) 3)
    (if close-points
      (do (setv (, i dist) (closest-pt (list (replace vor.vertices hole))))
          (if (< dist max-d)
            (.pop hole i)
            (setv close-points False)))
      (do (setv (, i angle) (widest-angle (list (replace vor.vertices hole))))
          (if (> angle min-angle)
            (.pop hole i)
            (break)))))
  hole)

(defn copy [list-of-lists]
  (list (map list list-of-lists)))

(defn get-wedge-holes-dist-angle [holes vor max-d min-angle]
  (if (pos? (len holes))
    (let [[entire-holes (copy holes)]
          [reduced-holes (list (map (fn [x] (remove-pts x vor max-d min-angle)) entire-holes))]
          [filtered (list (remove (fn [hole] (> (len (first hole)) 3))
                                  (zip reduced-holes holes)))]]
     ;; (list (remove (fn [hole] (> (len hole) 3)) (map (fn [x] (remove-pts x vor max-d min-angle)) holes)))
      (, (list (map first filtered)) (list (map second filtered))))
    (, [] [])))

(defn get-wedge-holes [holes vor config]
  (if (= (:hole-choser config) :area)
    (get-wedge-holes-area holes vor (:hole-choser-area-error config))
    (get-wedge-holes-dist-angle holes vor
                                (:hole-choser-max-d config)
                                (:hole-choser-min-angle config))))

;; (print (list (map fit-triangle true-holes)))
;; (print (list (map fit-triangle false-holes)))

;; (print (list (map (fn [x] (remove-pts (list (butlast x)))) false-holes)))

(def true-holes
  (array [[ [77.5176458238 44.6828993711]
            [82.8421137119 41.1319698207]
            [86.5771271012 46.4169850076]
            [86.0295802862 48.1089458024]
            [77.5176458238 44.6828993711]]
          [[115.351793478 12.3693285736]
           [117.222113814 18.9506303787]
           [124.256081903 11.8769560872]
           [118.293288624 10.4324921568]
           [115.351793478 12.3693285736]]]))

(def false-holes
  (array [[[177.048593727 38.8281980268]
           [175.183789555 37.0851391576]
           [180.033525113 28.8131659959]
           [185.829289552 25.9256470645]
           [202.789389693 26.1750889713]
           [201.644148337 39.030472472]
           [177.048593727 38.8281980268]]]))

(defn fit-triangle2 [hole vor] ;; compares areas of triangle and original hole, test!!!
  (if (> (len (list hole)) 3)
    (let [[entire-area (shoelace (list (replace vor.vertices hole)))]
          [tri-inds (list (it.combinations hole 3))]
          [tri-areas (list (map shoelace (replace2d vor.vertices tri-inds)))]
          [parts-list (list (map (fn [idx]
                              (list (map (fn [s e]
                                           (list (it.chain (it.takewhile
                                                            (fn [y] (not (= e y)))
                                                            (it.dropwhile (fn [x] (not (= s x)))
                                                                          (cycle hole)))
                                                           [e])))
                                         idx (roll idx -1))))
                                 tri-inds))]
          [limbs-list (list (map (fn [parts] ;; return chosen part?
                              (list (map (fn [nodes]
                                           (list (map (fn [n1 n2] [n1 n2])
                                                      (rest nodes) (rest (np.roll nodes -1)))))
                                         parts)))
                                 parts-list))]
          [abs-errs (list (map (fn [parts]
                                 (reduce add (map shoelace (replace2d vor.vertices parts))))
                               parts-list))]

          [min-ind (np.argmin abs-errs)]
          [rel-err (/ (get abs-errs min-ind)
                      (+ entire-area (get abs-errs min-ind)))]
          [best-triangle (list (get tri-inds min-ind))]]
      (dict [[:triangle best-triangle] [:original hole]
             [:rel-err rel-err] [:abs-err (get abs-errs min-ind)]]))
    (if (< (len (list hole)) 3)
      (dict [[:triangle None] [:original hole] [:rel-err 1] [:abs-err inf]])
      (dict [[:triangle hole] [:original hole] [:rel-err 0] [:abs-err 0]]))))

(defn get-path-edges [limb &optional [circular False]]
  (if circular
    (list (map (fn [i1 i2] (, i1 i2)) limb (roll limb -1)))
    (list (map (fn [i1 i2] (, i1 i2)) (butlast limb) (rest limb)))))

(defn get-wedge-holes2 [holes vor config] ;; iterative adding of holes, test!!!
  (if (empty? holes)
    (, [] [])
    (let [[entire-holes (copy holes)]
          [tri-fit (list (map (fn [hole] (fit-triangle2 hole vor)) entire-holes))]
          [used (get-undir-adj-matrix [] (len vor.vertices))]
          [chosen []]]
      (for [tri (sorted tri-fit None (fn [x] (get x :rel-err)))]
        (setv edges (get-path-edges (:original tri) True))
        (when (and (zero? (reduce add (map (fn [e] (get used (tuple e))) edges)))
                   (< (:rel-err tri) (:hole-choser-area-error config)))
          (list (map (fn [e] (assoc used (tuple e) 1) (assoc used (tuple (reversed e)) 1))
                     edges))
          (.append chosen tri)))
      (, (list (map (fn [x] (get x :triangle)) chosen)) (list (map (fn [x] (get x :original)) chosen))))))

(defn get-neighbours [graph node]
  (second (nonzero (get graph node))))

(defn get-incident-edges [graph node]
  "Operates on undirected graph as adjancency matrix"
  (map (fn [x] (get graph (, node x))) (get-neighbours graph node)))

(defn get-neighbours-except [graph node exceptions]
  (remove (fn [x] (in x exceptions)) (get-neighbours graph node)))



;; (angle-between [1 0] [-1 0.5]) ;; 180

;; (angle-between [1 1] [-1 -1] )

(defn wedge-skeleton [vertices edge-map vor paths center &optional whole-hole] ;; test with hole
  "Returns 3 limbs of wedge skeleton, center as coordinate"
  ;; center is mean of vertices if wedge has hole, isec when not"
  (list (map (fn [ind]
               (trace-skeleton-greedy1 edge-map vor ind center ;;urspr dir inst of center
                                      (if (none? whole-hole) paths
                                          (let [[h (Path (list (replace vor.vertices whole-hole)))]]
                                            (list (remove (fn [p] (.contains_path h p)) paths))))))
             vertices)))

(defn wedge-skeleton2 [vertices edge-map vor paths center &optional whole-hole] ;; test with hole
  "Returns 3 limbs of wedge skeleton, center as coordinate"
  ;; center is mean of vertices if wedge has hole, isec when not"
  (list (map (fn [ind]
               (trace-skeleton-greedy2 edge-map vor ind center ;;urspr dir inst of center
                                      (if (none? whole-hole) paths
                                          (let [[h (Path (list (replace vor.vertices whole-hole)))]]
                                            (list (remove (fn [p] (.contains_path h p)) paths))))))
             vertices)))

(defn trace-skeleton-greedy1 [edge-map vor start center paths &optional hole] ;; hole not used
  "Returns vertex list"
  (setv dir (- (get vor.vertices start) center))
  (setv path [start])
  (setv oldlength 0)
  (setv candidates (list (.difference (get edge-map start) path)))
  (while (not (empty? candidates))
    (setv angles (list (map (fn [x] (angle-between dir (- x center) 0) )
                            (replace vor.vertices candidates) )))
    (setv new-coordinate (get vor.vertices (get candidates (np.argmin angles))))
    (setv newlength (np.linalg.norm (-  new-coordinate (get vor.vertices start))))
    (setv cur-edge (Path [(get vor.vertices start) new-coordinate] [1 2]))
    ;; 1:moveto 2:lineto; problem not solvable by replacement with center
    (when (or (np.any (list (map (fn [p] (.intersects_path cur-edge p False)) paths)))
              (> oldlength newlength)
              (> (min angles) 70))
      (break))
    (setv oldlength newlength)
    (.append path (get candidates (np.argmin angles)))
    (setv candidates (list (.difference (get edge-map (last path)) path))))
  path)

(defn trace-skeleton-greedy2 [edge-map vor start center paths]
  "Returns vertex list"
  (setv dir (- (get vor.vertices start) center))
  (setv path [start])
  (setv oldlength 0)
  (setv candidates (list (.difference (get edge-map start) path)))
  (while (not (empty? candidates))
    (setv angles (list (map (fn [x] (angle-between dir (- x center) np.inf) )
                            (replace vor.vertices candidates) )))
    (setv new-coordinate (get vor.vertices (get candidates (np.argmin angles))))
    (setv newlength (np.linalg.norm (- new-coordinate (get vor.vertices start))))
    (setv cur-edge (Path [center ;; only change from ts1
                          new-coordinate] [1 2]))
    ;; 1:moveto 2:lineto; problem not solvable by replacement with center
    (when (or (np.any (list (map (fn [p] (.intersects_path cur-edge p False)) paths)))
              (> oldlength newlength)
              (> (min angles) 45))
      (break))
    (setv oldlength newlength)
    (.append path (get candidates (np.argmin angles)))
    (setv candidates (list (.difference (get edge-map (last path)) path))))
  path)

(defn trace-skeleton-greedy3 [edge-map vor start center paths]
  "Returns vertex list, based on 1, allows temporal path crossings"
  (setv dir (- (get vor.vertices start) center))
  (setv path [start])
  (setv oldlength 0)
  (setv candidates (list (.difference (get edge-map start) path)))
  (setv isec-steps 0)
  (while (not (empty? candidates))
    (setv angles (list (map (fn [x] (angle-between dir (- x center) np.inf) )
                            (replace vor.vertices candidates) )))
    (setv new-coordinate (get vor.vertices (get candidates (np.argmin angles))))
    (setv newlength (np.linalg.norm (- new-coordinate (get vor.vertices start))))
    (setv cur-edge (Path [(get vor.vertices start) new-coordinate] [1 2]))
    (when (or (> oldlength newlength)
              (> (min angles) 40)) ;; TODO: constrain both: the global angle and local angle
      (break))
    (setv oldlength newlength)
    (.append path (get candidates (np.argmin angles)))
    (if (np.any (list (map (fn [p] (.intersects_path cur-edge p False)) paths)))
      (do ;;(print isec-steps)
          (setv isec-steps (inc isec-steps))
          (when (> isec-steps 1) (break) ))
      (setv isec-steps 0))

    (setv candidates (list (.difference (get edge-map (last path)) path))))

  ;;(print path)
  ;;(print (if (zero? isec-steps) path (list (take (- (len path) isec-steps) path))))
  (if (zero? isec-steps) path (list (take (- (len path) isec-steps) path)))
  )



(defn get-used-edges [skeletons &optional centers]
  (when (none? centers) (setv centers (* [None] (len skeletons))))
  (if (empty? skeletons) []
      (np.vstack (list (map (fn [skeleton center]
                           (np.vstack (list (remove empty? (map (fn [limb]
                                                               (list (map (fn [i1 i2] [i1 i2])
                                                                          (rest (if (none? center) limb (cons center limb)))
                                                                          (rest (np.roll (if (none? center)
                                                                                        limb
                                                                                        (cons center limb)) 1)))))
                                                             skeleton)))))
                         skeletons
                         centers)))))

(defn get-path-edges [limb &optional [circular False]]
  (if circular
    (list (map (fn [i1 i2] [i1 i2]) limb (roll limb -1)))
    (list (map (fn [i1 i2] [i1 i2]) (butlast limb) (rest limb)))))

(defn get-used-edges-holes [skeletons holes] ;;test
  (if (empty? skeletons) []
      (np.vstack (list (map (fn [skeleton hole]
                           (np.vstack [(np.vstack (list (map (fn [limb]
                                                         (if (< (len limb) 2)
                                                           (np.array (, 0 2))
                                                           (get-path-edges limb)))
                                                       skeleton)))
                                     (get-path-edges hole True)]))
                         skeletons
                         holes)))))

(defn find-wedges2 [edges vor filtered-edges ridge-points wedge-skeletons holes paths t]
  (let [[edge-map (get-edge-map edges)]
        [isecs (get-graph-isecs edge-map)]
        [used (get-undir-adj-matrix (get-used-edges-holes wedge-skeletons holes)
                                    (len vor.vertices))]

        [edge-triples (list (map (fn [i] (list (map (fn [v] (, i v))
                                                    (get edge-map i))))
                                 isecs))]

        [isecpts (list (map (fn [i] (list (replace vor.points
                                                   (set (apply it.chain
                                                               (map second (filter (fn [e] (in i (set (first e))))
                                                                                   (zip filtered-edges ridge-points))))))))
                            isecs))]
        [dists (list (map (fn [x] (list (map (fn [y] (np.linalg.norm (- y (np.mean x 0)))) x)))
                          isecpts))]
        [a (print "d" dists)]
        [weights1 (list (map np.mean dists))] ;; außerdem: max, min, median
        [inner-angles (list (map (fn [edges]
                                   (let [[angles (list (map (fn [e1 e2]
                                                              (let [[u (apply sub e1)]
                                                                    [v (apply sub e2)]
                                                                    [angle (angle-between u v)]]
                                                                (if (neg? (.tolist (cross u v)))
                                                                  (* -1 angle) angle))) ;; urspr: 360-angle
                                                            (replace2d vor.vertices edges)
                                                            (roll (list (replace2d vor.vertices edges)) -1 0)))]]
                                     (list (filter (fn [x] (= (last x) (min angles)))
                                                   [angles (roll angles 1) (roll angles 2)]))))
                                 edge-triples))] ;; rightmost-angle = right angle = smallest angle
        [weights2 (list (map (fn [alphas] (int (and (or (all (negative (array alphas)))
                                                        (all (positive (array alphas))))
                                                    (all (> (abs (array alphas)) 90)))))
                             inner-angles))]
        [weights (zipwith mul weights1 weights2)]]
    (print weights)
    (setv ind-weights (sorted (remove (fn [x] (< (second x) t))
                                      (zip isecs weights)) None second True))
    (setv skeletons [])
    (setv centers [])
    (while (not (empty? ind-weights))
      (setv n-used (list (map (fn [i] (count-nonzero (list (map (fn [n] (get used (, (first i) n)))
                                                                (get edge-map (first i))))))
                              ind-weights)))
      (setv node-weight-used (first (sorted (zip ind-weights n-used) None second)))

      (when (= 3 (second node-weight-used)) (break))

      (setv node (first (first node-weight-used)))
      (.append skeletons (wedge-skeleton2 (get edge-map node) edge-map vor paths (get vor.vertices node)))
      (.append centers (get vor.vertices node))

      (list (map (fn [e] (assoc used (tuple e) 1) (assoc used (tuple (reversed e)) 1))
                 (get-used-edges [(last skeletons)] [node])))
      (setv ind-weights (list (remove (fn [x] (= node (first x))) ind-weights))))
    (, skeletons centers)))

(defn find-wedges [edges vor filtered-edges ridge-points wedge-skeletons paths]
  (let [[edge-map (get-edge-map edges)]
        [isecs (get-graph-isecs edge-map)]
        [used (get-undir-adj-matrix (get-used-edges wedge-skeletons)
                                    (len vor.vertices))]

        [edge-triples (list (map (fn [i] (list (map (fn [v] (, i v))
                                                    (get edge-map i))))
                                 isecs))]

        [isecpts (list (map (fn [i] (list (replace vor.points
                                                   (set (apply it.chain
                                                               (map second (filter (fn [e] (in i (set (first e))))
                                                                                   (zip filtered-edges ridge-points))))))))
                            isecs))]
        [dists (list (map (fn [x] (list (map (fn [y] (linalg.norm (- y (mean x 0)))) x)))
                          isecpts))]
        [a (print "d" dists)]
        [weights1 (list (map mean dists))] ;; außerdem: max, min, median
        [inner-angles (list (map (fn [edges]
                                   (let [[angles (list (map (fn [e1 e2]
                                                              (let [[u (apply sub e1)]
                                                                    [v (apply sub e2)]
                                                                    [angle (angle-between u v)]]
                                                                (if (neg? (.tolist (np.cross u v)))
                                                                  (* -1 angle) angle))) ;; urspr: 360-angle
                                                            (replace2d vor.vertices edges)
                                                            (np.roll (list (replace2d vor.vertices edges)) -1 0)))]]
                                     (list (filter (fn [x] (= (last x) (min angles)))
                                                   [angles (np.roll angles 1) (np.roll angles 2)]))))
                                 edge-triples))] ;; rightmost-angle = right angle = smallest angle
        [weights2 (list (map (fn [alphas] (int (and (or (np.all (negative (np.array alphas)))
                                                        (np.all (positive (np.array alphas))))
                                                    (np.all (> (abs (np.array alphas)) 70)))))
                             inner-angles))]
        [weights (zipwith mul weights1 weights2)]]

    (setv ind-weights (sorted (remove (fn [x] (zero? (second x))) (zip isecs weights)) None second True))
    (setv skeletons [])
    (setv centers [])
    (while (not (empty? ind-weights))
      (setv n-used (list (map (fn [i] (count-nonzero (list (map (fn [n] (get used (, (first i) n)))
                                                                (get edge-map (first i))))))
                              ind-weights)))
      (setv node-weight-used (first (sorted (zip ind-weights n-used) None second)))

      (when (= 3 (second node-weight-used)) (break))

      (setv node (first (first node-weight-used)))
      (.append skeletons (wedge-skeleton (get edge-map node) edge-map vor paths (get vor.vertices node)))
      (.append centers (get vor.vertices node))

      (list (map (fn [e] (assoc used (tuple e) 1) (assoc used (tuple (reversed e)) 1))
                 (get-used-edges [(last skeletons)] [node])))
      (setv ind-weights (list (remove (fn [x] (= node (first x))) ind-weights))))
    (, skeletons centers)))


(defn wedge-from-skeleton [limbs vor &optional center [mode :double]] ;; TODO: replace center by hole; hole with 1 point is center
  "Limbs consisting of vertex lists, modes :double or :quad :vert :vert-quad"
  ;; limbs may not start with same vertex if modes :vert or :vert-quad are used
  (when (none? center)
    (setv center (np.mean (list (replace vor.vertices (map first limbs))) 0)))
  (let [[corners (list (map last limbs))]
        [startpts (list (map first limbs))]]
    (cond [(= mode :double)
           (list (map (fn [c1 c2]  [(get vor.vertices c1) center
                                    center (get vor.vertices c2)])
                      corners (np.roll corners -1)))]
          [(= mode :quad)
           (list (map (fn [c1 c2] (elevate-bezier [(get vor.vertices c1)
                                                   center
                                                   (get vor.vertices c2)]))
                      corners (np.roll corners -1)))]
          [(= mode :vert)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                      (get vor.vertices s3)
                                      (get vor.vertices s3)
                                      (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :vert-quad)
           (list (map (fn [c1 c2 s3] (elevate-bezier [(get vor.vertices c1)
                                                      (get vor.vertices s3)
                                                      (get vor.vertices c2)]))
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :mid)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                   (np.mean [(get vor.vertices c1) (get vor.vertices s3)] 0)
                                   (np.mean [(get vor.vertices c2) (get vor.vertices s3)] 0)
                                   (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :mid-cross)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                   (np.mean [(get vor.vertices c2) (get vor.vertices s3)] 0)
                                   (np.mean [(get vor.vertices c1) (get vor.vertices s3)] 0)
                                   (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]

          [True (print "unknown reconstruction mode")])))

(defn connect-points-rec [p1 p2 max-dist]
  (if (> (np.linalg.norm (- p1 p2)) max-dist)
    (it.chain (connect-points-rec p1 (/ (+ p1 p2) 2) max-dist)
              (rest (connect-points-rec (/ (+ p1 p2) 2) p2 max-dist)))
    [p1 p2]))

(defn connect-points [ptlist max-dist]
  (list (apply it.chain (list (map (fn [p1 p2] (connect-points-rec p1 p2 max-dist))
                                   (butlast ptlist) (rest ptlist))))))


(defn wedge-from-skeleton2 [limbs vor &optional hole [mode :double]] ;; TODO: replace center by hole; hole with 1 point is center
  "Limbs consisting of vertex lists, modes :double or :quad :vert :vert-quad"
  ;; limbs may not start with same vertex if modes :vert or :vert-quad are used
  (if (= 1 (len hole))
    (setv center (get vor.vertices (first hole)))
    (setv center (np.mean (list (replace vor.vertices (map first limbs))) 0)))
  (let [[corners (list (map last limbs))]
        [startpts (list (map first limbs))]]
    (cond [(= mode :double)
           (list (map (fn [c1 c2]  [(get vor.vertices c1) center
                                    center (get vor.vertices c2)])
                      corners (np.roll corners -1)))]
          [(= mode :quad)
           (list (map (fn [c1 c2] (elevate-bezier [(get vor.vertices c1)
                                                   center
                                                   (get vor.vertices c2)]))
                      corners (np.roll corners -1)))]
          [(= mode :vert)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                      (get vor.vertices s3)
                                      (get vor.vertices s3)
                                      (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :vert-quad)
           (list (map (fn [c1 c2 s3] (elevate-bezier [(get vor.vertices c1)
                                                      (get vor.vertices s3)
                                                      (get vor.vertices c2)]))
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :mid)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                   (np.mean [(get vor.vertices c1) (get vor.vertices s3)] 0)
                                   (np.mean [(get vor.vertices c2) (get vor.vertices s3)] 0)
                                   (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :mid-cross)
           (list (map (fn [c1 c2 s3] [(get vor.vertices c1)
                                   (np.mean [(get vor.vertices c2) (get vor.vertices s3)] 0)
                                   (np.mean [(get vor.vertices c1) (get vor.vertices s3)] 0)
                                   (get vor.vertices c2)])
                      corners (np.roll corners -1) (np.roll startpts 1)))]
          [(= mode :bezier-fit) ;; test
           (let [[hole-verts (list (map first limbs))]
                 [parts (list (map (fn [s e]
                                     (list (rest (it.takewhile
                                                  (fn [y] (not (= e y)))
                                                  (it.dropwhile (fn [x] (not (= s x)))
                                                                (cycle hole))))))
                                   hole-verts (roll hole-verts -1)))]
                 [beziers (list (map (fn [l1 l2 part]
                                       (bezier-fit (connect-points (list (replace vor.vertices
                                                                                  (it.chain (reversed l1) part l2))) 5.0))) ;; allow adjusting
                                     limbs
                                     (roll limbs -1 0)
                                     parts))]]
             beziers)]
          [True (print "unknown reconstruction mode")])))
