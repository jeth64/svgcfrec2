(import [numpy [array transpose mean sum dot cross sin cos roll radians amin amax all
                isfinite isinf around vstack arccos sign sqrt]]
        [numpy.linalg [solve LinAlgError det norm]]
        [operator [sub]]
        [itertools :as it]
        [tools [*]])


(defn stretch [line factor]
  "Stretch given line in the 'end' direction"
  (let [[l (array (list line))]]
    (array [(first l) (+ (first l) (* -1 factor (apply sub l)))])))

(defn isec [lines &optional [mode :inf]] ; test modes
  "Intersection of lines (returns None if no intersection and [inf inf] if parallel"
  ;; modes:
  ;;  - :inf unbounded line to both sides
  ;;  - :inf-left unbounded in the direction of the first line point
  ;;  - :inf-right unbounded in the direction of the second line point
  ;;  - :segments if bounded in both directions
  (try (let [[l1 (array (first (list lines)))]
             [l2 (array (second (list lines)))]
             [r1 (- (last l1) (first l1))]
             [r2 (- (last l2) (first l2))]
             [b (- (first l2) (first l1))]
             [x (solve (transpose (array [r1 (- r2)])) b)]
             [point (+ (first l1) (* (first x) r1))]]
         (cond [(= mode :inf) point]
               [(= mode :inf-left) (if (all (< x 1)) point None)]
               [(= mode :inf-right) (if (all (< 0 x)) point None)]
               [(= mode :segments) (if (and (all (< 0 x)) (all (< x 1))) point None)]))
       (catch [e LinAlgError]
         (array [(float "inf") (float "inf")]))))

;;(isec (array [[[0 0] [1 1]] [[0 1] [1 0]]]))
;;(isec (array [[[0 0] [1 1]] [[0 1] [1 0]]]) :segments)
;;(isec (array [[[0 0] [0.4 0.4]] [[0 1] [1 0]]]))
;;(none? (isec (array [[[0 0] [0.4 0.4]] [[0 1] [1 0]]]) :segments)) ;; none
;;(none? (isec (array [[[0 0] [0.4 0.4]] [[0 1] [1 0]]]) :inf-right))
;;(none? (isec (array [[[0 0] [0.4 0.4]] [[0 1] [1 0]]]) :inf-left)) ;; none
;;(none? (isec (array [[[0 0] [0.4 0.4]] [[0 1] [1 0]]]) :inf))


(defn cross-line [line]
  "Returns line rotated around midpoint by 90 degree"
  (let [[l (array (list line))]
        [midpt (mean l 0) ]
        [rotated-dir (dot [[0 -1] [1 0]]
                          (apply sub l))]]
    (array [(+ midpt (/ rotated-dir 2))
            (- midpt (/ rotated-dir 2))])))

(defn rotation-matrix [angle]
  "Return rotation matrix of angle given in radian (counter-clockwise)"
  (let [[c (cos angle)]
        [s (sin angle)]]
    (array [[c (- s)] [s c]])))

(defn rotate [angle vectors]
  "Rotate all vectors of list counter-clockwise (angle in radian)"
  (dot vectors (rotation-matrix (- angle))))

(defn bounding-box-v [points dir-vec]
;; "Returns vertices of rotated rectangle containing given points.
;;   The direction is given by a normalized vector "
  (let [[angle (arccos (dot [0 1] dir-vec))]
        [p0 (mean points)]
        [normalized (- points p0)]
        [rotated (rotate angle normalized)]]
    (+ (rotate (- angle)
               (list-comp i [i (apply it.product (transpose [(amin rotated 0)
                                                             (amax rotated 0)]))]))
       p0)))

(defn bounding-box [points degree]
  "Returns vertices of rotated rectangle containing given points"
  (let [[p0 (mean points)]
        [whitened (- points p0)]
        [angle (radians degree)]
        [rotated (rotate angle whitened)]]
    (+ (rotate (- angle)
               (list-comp i [i (apply it.product (transpose [(amin rotated 0)
                                                             (amax rotated 0)]))]))
       p0)))

(defn points-in-triangle? [vertices points]
  "States if all given points are within the triangle defined by the vertices"
  (if (all (isfinite points))
    (let [[v-dirs (list (map sub vertices (roll vertices 1 0)))]
          [v2-dirs (list (map sub vertices (roll vertices 2 0)))]
          [p-dirs (list (map (fn [z] (list (map (fn [x] (sub z x)) points))) vertices))]]
      (all (list (map (fn [x y z] (list (map (fn [p] (= (pos? (first (array [(cross x y)])))
                                                        (pos? (first (array [(cross x p)])))))
                                             z)))
                      v-dirs v2-dirs p-dirs))))
    False))

(defn tri-area [vs]
  (abs (/ (+ (* (first (first vs))
                (- (second (second vs)) (second (last vs))))
             (* (first (second vs))
                (- (second (last vs)) (second (first vs))))
             (* (first (last vs))
                (- (second (first vs)) (second (second vs))))) 2)))

(defn heron [vertices]
  "Calculates the area of an triangle with Heron's formula"
  (if (> (len vertices) 2)
    (if (> (len vertices) 3)
      (throw (Exception "Vertex list must not have more than 3 entries"))
      (do (setv sides (list (map (fn [u v] (norm (- u v))) vertices (roll vertices 1 0))))
          (let [[s (/ (reduce add sides) 2)]]
            (sqrt (* s (- s (first sides)) (- s (second sides)) (- s (last sides)))))))
    0))

(defn shoelace [vertices]
  "Calculates area for simple polygons with shoelace formula"
  (if (> (len vertices) 2)
    (let [[XY (transpose vertices)]]
      (/ (abs (reduce add (map (fn [xyi xyi+1] (det (vstack [xyi xyi+1])))
                               vertices
                               (roll vertices 1 0))))
         2.0))
    0))

(defn area [vertices]
  "Area for simple self-intersecting polygons"
  ;; Splits path recursively at intersections; results similar to nonzero rule of path filling
  (setv segments (array (list (zip vertices (roll vertices -1 0)))))
  (setv intersecting? False)
  (for [ij (it.combinations (range (len segments)) 2)]
    (setv S (isec [(get segments (first ij))
                   (get segments (second ij))]
                  :segments))
    (when (not (or (none? S) (any (isinf S))))
      (setv outer (array (list (it.chain (take (inc (first ij)) vertices)
                                     [S]
                                     (drop (inc (second ij)) vertices)))))
      (setv middle (array (list (it.chain [S]
                                       (drop (inc (first ij)) (take (inc (second ij)) vertices))))))
      (setv intersecting? True)
      (break)))
  (if intersecting?
    (+ (area outer) (area middle))
    (shoelace vertices)))

(tri-area [[77 44] [82 41] [86 46]]) ; 18.5
(heron [[77 44] [82 41] [86 46]])
(area [[77 44] [82 41] [86 46]])
(shoelace [[77 44] [82 41] [86 46]])

;; (around (area (array [[0 0] [1 1] [2 1] [3 0]])) 2)
;; => 2.0 ???
;; (around (area (array [[0 0] [1 1] [2.5 -0.5] [3 -1] [4 0]])) 2)
;; => 2.0
;; (around (area (array [[0 0] [1 1] [2.5 -0.5] [3 -1]  [5 -1] [4 0]])) 2)
;; => 3.0


;;
;; algorthms for polygons with all but first and last vertex lying in the area
;; between the lines through the first and last vertex respectively placed orthogonal
;; to the line through both points
;;

  ;; only for given simple cases where all vertices are
(defn get-xy-segments-of-polygon [vertices]
  (let [[base (- (last vertices) (first vertices))]
       ;; [a (print base)]
        [XakkY (transpose (list (map (fn [v] (let [[diff-v (- v (first vertices))]
                                                   [alpha (arccos (dot (/ diff-v (norm diff-v))
                                                                       (/ base (norm base))))]
                                                   [lv (norm diff-v)
                                                    ]]
                                               [(* lv (cos alpha)) (* lv (sign (cross diff-v base))
                                                                      (sin alpha))]))
                                     (rest vertices))))]
        [X (list (it.chain [(first (first XakkY))] (map (fn [xi xi-1] (- xi xi-1))
                                         (rest (first XakkY)) (butlast (first XakkY)))))]]
    [X (second XakkY)]))

;; (around (get-xy-segments-of-polygon (array [[1 2] [2 2] [2 1.5] [2 1]])) 2)
;;=> [[0.70 0.35 0.35][0.70 0.35 0]]


(defn get-polygon-area [vertices &optional [normalized False] [absolute False]]
  ;; only for given simple cases where all vertices are in area between lines through first and last vertex respectively orthogonal to line through both points
  (let [[XY (get-xy-segments-of-polygon vertices)]
        [area (reduce add (map (fn [x y-1 y]
                                    ;; (print "input" [x y-1 y])
                                    ;; (print absolute)
                                    ;; (print (neg? (* y y-1)))
                                     (if absolute
                                       (if (neg? (* y y-1))
                                         (let [[S (isec [[[0 0] [x 0]] [[0 y-1] [x y]]] )]]
                                       ;;    (print "S" S)
                                       ;;    (print "x" x)
                                           (+ (abs (/ (* y-1 (first S)) 2))
                                              (abs (/ (* y (- x (first S))) 2))))
                                         (abs (* x (/ (+ y y-1) 2))))
                                       (* x (/ (+ y y-1) 2))))
                               (first XY)
                               (cons 0 (butlast (second XY)))
                               (second XY)))]]
    (if normalized
      (/ area (pow (sum (first XY)) 2))
      area)))

;;(get-polygon-area (array [[0 0] [1 1] [2 1] [3 0]]) True)

;; (around (get-polygon-area (array [[0 0] [1 1] [2 1] [3 0]]) False) 2)
;; => 2.0

;; (around (get-polygon-area (array [[0 0] [1 1] [2 0] [2.5 -0.5] [3 -1] [4 0]]) False True) 2)
;; => 2.0

;; (around (get-polygon-area (array [[0 0] [1 1] [2 1] [3 0]]) True) 1)
;; => 0.2



(defn get-mean-polygon-dists [vertices &optional [normalized False]] ;test, auch max
  (let [[XY (get-xy-segments-of-polygon vertices)]
        [meanY (mean (array (list (butlast (second XY)))))]]
    (if normalized
      (/ meanY (sum (first XY)))
      meanY)))

;; (get-mean-polygon-dists (array [[1 2] [2 2] [2 1.5] [2 1]]))
