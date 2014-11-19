(import [numpy [*]]
        [numpy.linalg [solve LinAlgError]]
        [operator [sub]]
        [itertools :as it]
        [tools [*]])


(defn stretch [line factor]
  "Stretch given line in the 'end' direction"
  (let [[l (array (list line))]]
    (array [(first l) (+ (first l) (* -1 factor (apply sub l)))])))


(defn isec [lines]
  "Intersection of lines that are infinite in the 'end' direction"
  (try (let [[l1 (array (first (list lines)))]
             [l2 (array (second (list lines)))]
             [r1 (- (last l1) (first l1))]
             [r2 (- (last l2) (first l2))]
             [b (- (first l2) (first l1))]
             [x (solve (transpose (array [r1 (- r2)])) b)]]
         (if (every? pos? x)
           (+ (first l1) (* (first x) r1))
       (array [(float "inf") (float "inf")])))
       (catch [e LinAlgError]
         (array [(float "inf") (float "inf")]))))


(defn cross-line [line]
  "Returns line rotated around midpoint by 90 degree"
  (let [[l (array (list line))]
        [midpt (mean l 0) ]
        [rotated-dir (dot [[0 -1] [1 0]]
                          (apply sub l))]]
    (array [(+ midpt (/ rotated-dir 2))
            (- midpt (/ rotated-dir 2))])))

(defn rotation-matrix [angle]
  "Return rotation matrix of angle given in radian."
  (let [[c (cos angle)]
        [s (sin angle)]]
    (array [[c s] [(- s) c]])))

(defn rotate [angle vectors]
  "Rotate all vectors of list clockwise (angle in radian)"
  (dot vectors (transpose (rotation-matrix angle))))

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
    (let [[v-dirs (list (map sub vertices (roll vertices 1)))]
          [v2-dirs (list (map sub vertices (roll vertices 2)))]
          [p-dirs (list (map (fn [z] (list (map (fn [x] (sub z x)) points))) vertices))]]
      (all (list (map (fn [x y z] (list (map (fn [p] (= (pos? (cross x y))
                                                        (pos? (cross x p))))
                                             z)))
                      v-dirs v2-dirs p-dirs))))
    False))
