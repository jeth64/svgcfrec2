(import [xml.etree.ElementTree [parse SubElement]]
        [svg [get-paths matrix2path]]
        [numpy [*]]
        [lisp-tools [*]]
        [operator [sub]]
        [itertools [product]])

(defn angle-bisector [cubic-bezier]
  (let [[point (dict (zip [:p1 :p2 :p3 :p4] cubic-bezier))]]
    (/ (array [ (+ (get point :p1) (get point :p4))
                (+ (get point :p2) (get point :p3)) ]) 2)))


(defn add-path [parent attribs path] ;; side-effects
  (.update attribs {"d" (matrix2path path)})
  (SubElement parent "ns0:path" attribs))

(defn add-line [parent attribs line] ;; side-effects
  (.update attribs (zip ["x1" "y1" "x2" "y2"] (map str (flatten line))))
  (SubElement parent "ns0:line" attribs))

(defn add-circle [parent attribs point]
  (.update attribs (zip ["cx" "cy"] (map str point)))
  (SubElement parent "ns0:circle" attribs))


(defn svg? [root]
  (= (.join "" (drop (- (len root.tag) 3) root.tag)) "svg"))

(defn choose-pts [coll xs ys];;working? necessary? ;; strange behaviour of get with array
  (get coll (.tolist (.transpose (array (list (product xs ys)))))))

(defn line-direction-vec [point-pair]
  (apply sub point-pair))

(defn convex? [part]
   (let [[a (line-direction-vec (angle-bisector part))]
         [b (line-direction-vec (slice part 0 4 3))]]
     (pos? (cross a b))))


(defn approx-points [points n]
  (rest (reductions (fn [P [i Q]]
                      (weighted-sum [Q P]
                                    [(/ n (- n i)) (/ (- i) (- n i))]))
                    [0 0]
                    (transpose [(range n) points]))))
(defn factor [i n]
  (* (Math/pow 2 (- 1 (* 2 n)))
     (reduce #(+ %1 (choose (* 2 n) (* 2 %2)))
             0 (range (inc i)))))

(defn reduce-degree [control-points]
  (let [n (dec (count control-points))
        Pr (approx-points control-points n)
        Pl (reverse (approx-points (reverse control-points) n))]
    (map #(let [lambda (factor %3 n)]
            (weighted-sum [%1 %2] [(- 1 lambda) lambda]))
         Pr Pl (range n) )))

;;main
(do (setv tree (parse ;"one-string-BC-0.45-2-1-0.2.svg"
                                "test3.svg"
                      ))
    (setv root (.getroot tree))
    (assert (svg? root))
    (setv outfile "out.svg")
    (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
    (setv paths (get-paths root namespace))
    (for [path paths]
      (let [[n (len path)]
            [lines (map angle-bisector path)]
            [lends (map (fn [part] (slice part 0 4 3)) path)]
            ;;[a (print lends)]
            [end-points (choose-pts path (range 0 n) [0])]
            [mid-points (choose-pts path (range 0 n) [1 2])]
            [conv (map angle-bisector (filter convex? path))]
       ;     [a (print conv)]
            ]
        (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} path)
       ; (for [line lends] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.3"} line))
        ;;(for [line lines] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
        (for [line conv] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
        (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
        (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
        ))
    (.write tree outfile)
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

;;choose convex segments: with point

;;[[11.9 5.1] ... [12.6 2.5]]
;;[[ 11.904849  ,   7.7520879 ],[ 22.912619  ,   9.3858879 ],
;; [ 23.272079  ,  11.243688  ],[ 12.456489  ,  10.603088  ]]

;;[[2.2 10.35] .. [4.5 9.6]]
;;


;; etwa 4 convex :2 concav : 3 gerade, angle bisector

;; (defn convex-parts [path])

(map (fn [x] (partition 2 2 x))
     path)

(map identity path)

(setv v [[[1 1] [2 2] [3 3] [4 4]]
         [[1 1] [2 2] [3 3] [4 4]]])



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
