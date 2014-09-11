(import [xml.etree.ElementTree [parse SubElement]]
        [svg [get-paths matrix2path add-path add-line add-circle svg?]]
        [bezier [merge-beziers dist de-casteljau-split de-casteljau]]
        [numpy [*]]
        [lisp-tools [*]]
        [operator [sub mul add]]
        [itertools [product groupby]]
        [scipy.spatial [ConvexHull]])


(defn angle-bisector [cubic-bezier]
  (let [[point (dict (zip [:p1 :p2 :p3 :p4] cubic-bezier))]]
    (/ (array [ (+ (get point :p1) (get point :p4))
                (+ (get point :p2) (get point :p3)) ]) 2)))


(defn choose-pts [coll xs ys];;working? necessary? ;; strange behaviour of get with array
  (get coll (.tolist (.transpose (array (list (product xs ys)))))))

(defn line-direction-vec [point-pair]
  (apply sub point-pair))

(defn convex? [part]
   (let [[a (line-direction-vec (angle-bisector part))]
         [b (line-direction-vec (slice part 0 4 3))]]
     (pos? (cross a b))))

(defn no-edge? [part1 part2]
  (let [[a (line-direction-vec [(first part1) (first part2)])]
        [b (line-direction-vec [(first part1) (second part2)])]]
    (= (pos? (cross a b)) (convex? part1))))

(defn transition-type [part1 part2]
  (if (convex? part1)
    (if (convex? part2)
      (if (no-edge? part1 part2) :convex-smooth :convex-edgy)
      :convex-concave)
    (if (convex? part2)
      :concave-convex
      (if (no-edge? part1 part2) :concave-smooth :concave-edgy))))


(defn partition-at-edge [parts]
  (let [[l [[(first parts)]]]]
    (for [consecutive (zip (butlast parts) (rest parts))]
      (if (apply no-edge? consecutive)
        (.append (last l) (last consecutive))
        (.append l [(last consecutive)])))
    l))

(defn try-merge [parts]
;  (print "merge")
  (if (= (len parts) 1)
    parts
    (map merge-beziers (partition-at-edge parts))))


;; "Splits segment such that in several parts
;;   'cut-off-attribute' can be
;;   - 'distance' : t is maximum distance between endpoints
;;   - 'pathlength' : t is maximum sum of point-to-point distances
;;   - 'curvature'y : t is maximum fraction of pathlength/distance "

(defn split-segment-to [cut-off-attribute t segment];;test
  (if (< (cond [(= cut-off-attribute "distance")
                (dist (first segment) (last segment))]
               [(= cut-off-attribute "pathlength")
                (reduce add (map dist (butlast segment) (rest segment)))]
               [(= cut-off-attribute "curvature")
                (/ (reduce add (map dist (butlast segment) (rest segment)))
                   (dist (first segment) (last segment)))])
         t)
    [segment]
    (let [[segments (de-casteljau-split 0.5 segment)]
          [new-segments (map (fn [s] (split-segment-to cut-off-attribute t s)) segments)]]
      (apply add new-segments))))


;;main
(do (setv tree (parse "more-strings-BC-0.3-2-1-0.2.svg"
      ;;                "testfiles/test1.svg"
                      ))
    (setv root (.getroot tree))
    (assert (svg? root))
    (setv outfile "out.svg")
    (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
    (setv paths (get-paths root namespace))
    (for [path paths]
      (let [[n (len path)]
          ;  [lines (map angle-bisector path)]
            [lends (map (fn [part] (slice part 0 4 3)) path)]
            [end-points (choose-pts path (range 0 n) [0])]
            [mid-points (choose-pts path (range 0 n) [1 2])]
            [c (reduce add (map (fn [x] (try-merge (list (second x)))) (filter first (groupby path convex?))) [])]
            [conv (map angle-bisector c)]
            [t 1.2]
            [new-path (array (reduce add (map (fn [x] (split-segment-to "curvature" t x)) path) []))]
        ;;    [new-path-ends (array (map first new-path))]
         ;   [b (print "new-path" (array (map first new-path)))]
         ;;   [ch (ConvexHull (array (map first new-path)))]
      ;      [y (print "y" (get new-path-ends (unique ch.simplices)))]
       ;     [a (print (unique ch.simplices))]
        ;    [x (print ch.ndim ch.nsimplex )]
            ]
        (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} path)
      ;;  (add-path root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} c)
        (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
        (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
       ; (add-path root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} path)
       ; (for [point (map first new-path)] (add-circle root {"fill" "orange" "r" "0.2"} point))
      ;  (for [point (get new-path-ends (unique ch.simplices))] (add-circle root {"fill" "black" "r" "0.2"} point))
      ;;  (for [point (choose-pts new-path (range 0 (len new-path)) [1 2])] (add-circle root {"fill" "black" "r" "0.2"} point))
       ; (for [line lends] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.3"} line))
        ;;(for [line lines] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
        (for [line conv] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
        ;;(for [line ] (add-line root {"fill" "none" "stroke" "black" "stroke-width" "0.3"} line))
        ))
     (.write tree outfile)
    None
    )



(setv part (map array [[ 11.904849     7.7520879 ] [ 22.912619    9.3858879 ] [ 23.272079  11.243688  ][ 12.456489 10.603088  ]]))


(let [[point (dict (zip [:p1 :p2 :p3 :p4] part))]
      [angle-bisector (/ (array [ (+ (get point :p1) (get point :p4))
                                  (+ (get point :p2) (get point :p3)) ]) 2) ]
      [bisector-dir-vec (apply sub angle-bisector)]
      [bisector-length (linalg.norm bisector-dir-vec)] ;; threshold
      [object-dir-vec (- (get point :p4) (get point :p1))]
      ]
  (print "angle-bisector" angle-bisector)
  (print "object-dir-vec" object-dir-vec)
  [angle-bisector (if (pos? (dot bisector-dir-vec object-dir-vec)) 1 0)])


(let [[angle-bisectors (map angle-bisector path)]
      [bisector-dir-vecs (map (fn [bis] (apply sub bis)) angle-bisectors)]
      [bisector-length (map linalg.norm bisector-dir-vecs)] ;; threshold
      [object-dir-vecs (map (fn [part] (apply sub (slice part 0 4 3))) path)]
      []])



(setv path (array [[[  0.16338942  11.341288  ]
                     [  0.49957942  10.465188  ]
                     [  0.99771942   8.6651879 ]
                     [  1.2703694    7.3412879 ]]
                    [[  1.2703694    7.3412879 ]
                     [  1.5430294    6.0173879 ]
                     [  1.9996794    3.8091879 ]
                     [  2.2851594    2.4341879 ]]
                    [[  2.2851594    2.4341879 ]
                     [  3.0451994   -1.2265121 ]
                     [  4.7661094   -0.61661206]
                     [  4.7661094    3.3133879 ]]
                    [[  4.7661094    3.3133879 ]
                     [  4.7661094    4.4397879 ]
                     [  4.7661094    5.5661879 ]
                     [  4.7661094    6.6925879 ]]
                    [[  4.7661094    6.6925879 ]
                     [  7.14568927   7.04575457]
                     [  9.52526913   7.39892123]
                     [ 11.904849     7.7520879 ]]
                    [[ 11.904849     7.7520879 ]
                     [ 22.912619     9.3858879 ]
                     [ 23.272079    11.243688  ]
                     [ 12.456489    10.603088  ]]
                    [[ 12.456489    10.603088  ]
                     [  5.3997394   10.185188  ]
                     [  3.5379594   10.365188  ]
                     [  3.0994994   11.507788  ]]
                    [[  3.0994994   11.507788  ]
                     [  2.7984494   12.292288  ]
                     [  1.8771294   12.934188  ]
                     [  1.0521294   12.934188  ]]
                    [[  1.0521294   12.934188  ]
                     [  0.09075942  12.934188  ]
                     [ -0.22838058  12.362188  ]
                     [  0.16338942  11.341288  ]]]))

(map angle-bisector (array v))
