(import [numpy [*]]
        [scipy.misc [comb]]
        [operator [sub mul add]]
        [lisp-tools [*]]
        [itertools [groupby]])


(defn dist [a b] (linalg.norm (- (array (list a))
                                 (array (list b)))))

(defn pt-coeffs [ts]
  (let [[t (array (list ts))]]
    (array (list-comp (* (comb 3 i)
                         (** t i)
                         (** (- 1 t) (- 3 i)))
                      [i (range 4)]))))

(defn eval-cubics [ts control-points]
  (dot (transpose (pt-coeffs ts)) (array (list control-points)) ))

(defn bezier-fit [ordered-pts]
  (let [[path-lengths (+ [0] (reductions add (map dist
                                                  (butlast ordered-pts)
                                                 (rest ordered-pts))))]
        [ts (array (list (map (fn [x] (/ x (last path-lengths))) path-lengths)))]
        [XT (pt-coeffs ts)]
        [control-points (dot (linalg.inv (dot XT (transpose XT)))
                             (dot XT ordered-pts))]]
    control-points))

(defn merge-beziers [cubic-beziers]
  (bezier-fit (vstack (list (map (fn [Cs] (eval-cubics (linspace 0 1 4 True) Cs))
                                 cubic-beziers)))))

;;(setv v (array [[0 1] [2 3] [4 5] [6 7]]))


;; Degree reduction: if needed needs to be corrected

(defn weighted-sum [points weights]
  (sum (array (map mul points weights))) 0)

(defn approx-points [points n]
  (rest (reductions (fn [P iQ]
                      (weighted-sum [(second iQ) P]
                                    [(/ n (- n (first iQ))) (/ (- (first iQ)) (- n (first iQ)))]))
                    [0 0]
                    (transpose [(range n) points]))))

(defn factor [i n]
  (* (Math/pow 2 (- 1 (* 2 n)))
     (reduce (fn [x y] (+ x (choose (* 2 n) (* 2 y))))
             0 (range (inc i)))))


(defn elevate-degree [p1 p2 p3] ;; delete?
  [p1 (/ (+ p1 p2 p2) 3) (/ (+ p2 p2 p3) 3) p3])

(defn elevate-bezier [control-points]
  "Elevates quadratic bezier curve to cubic"
  (let [[pts (array (list control-points))]]
    (array [(first pts)
            (/ (+ (first pts) (* 2 (second pts))) 3)
            (/ (+ (* 2 (second pts)) (last pts)) 3)
            (last pts)])))

(defn reduce-degree [control-points]
  (let [[n (dec (count control-points))]
        [Pr (approx-points control-points n)]
        [Pl (reverse (approx-points (reverse control-points) n))]]
    (map (fn [x y z] (let [[lmbd (factor z n)]]
                       (weighted-sum [x y] [(- 1 lmbd) lmbd])))
         Pr Pl (range n))))


;;
;; De-Casteljau algorithm
;;


(defn de-casteljau-rec [t control-points]
  (let [[coefs (array (last control-points))]]
    (if (< (len coefs) 2)
      control-points
      (+ control-points
         (de-casteljau-rec t [(list (map (fn [Pi Pinci] (+ (* Pi (- 1 t)) (* Pinci t)))
                                         (butlast coefs) (rest coefs)))])))))

(defn de-casteljau-split [t control-points];; test
  (let [[c-point-list (de-casteljau-rec t [control-points])]]
    [(array (map first c-point-list))
     (array (list (reversed (map last c-point-list))))]))

(defn de-casteljau [t control-points]
  (let [[c-point-list (de-casteljau-rec t [control-points])]]
    (array (first (last c-point-list)))))



;; helper functions

(defn line-direction-vec [point-pair]
  (apply sub point-pair))

(defn angle-bisector [cubic-bezier]
  (let [[point (dict (zip [:p1 :p2 :p3 :p4] cubic-bezier))]]
    (/ (array [ (+ (get point :p1) (get point :p4))
                (+ (get point :p2) (get point :p3)) ]) 2)))

(defn bezier-curvature [segment]
  (/ (reduce add (map dist (butlast segment) (rest segment)))
     (dist (first segment) (last segment))))

(defn bezier-curvature2 [segment1 segment2]
  (/ (dot (apply sub (take 2 segment1)) (apply sub (drop 2 segment2)))
     (* (linalg.norm segment1) (linalg.norm segment2))))


;;
;; high-level functions to manipulate poly-beziers
;;

(defn transition-type [part1 part2]
  (if (convex? part1)
    (if (convex? part2)
      (if (no-edge? part1 part2) :convex-smooth :convex-edgy)
      :convex-concave)
    (if (convex? part2)
      :concave-convex
      (if (no-edge? part1 part2) :concave-smooth :concave-edgy))))

(defn convex? [part]
   (let [[a (line-direction-vec (angle-bisector part))]
         [b (line-direction-vec (slice part 0 4 3))]]
     (pos? (cross a b))))

(defn no-edge? [part1 part2]
  (let [[a (line-direction-vec [(first part1) (first part2)])]
        [b (line-direction-vec [(first part1) (second part2)])]]
    (= (pos? (cross a b)) (convex? part1))))

(defn partition-at-edge [parts]
  (let [[l [[(first parts)]]]]
    (for [consecutive (zip (butlast parts) (rest parts))]
      (if (apply no-edge? consecutive)
        (.append (last l) (last consecutive))
        (.append l [(last consecutive)])))
    l))
;;
;; merge convex segments if they don't meet on an edge

(defn try-merge [parts]
  (if (= (len parts) 1)
    parts
    (map merge-beziers (partition-at-edge parts))))


(defn get-merged-convex [path]
  (let [[convex-segments (reduce add (map (fn [x] (try-merge (list (second x))))
                                          (filter first (groupby path convex?))) [])]]
    (if (< (dist (first (first convex-segments))
                 (last (last convex-segments)))
           0.5)
      (add (array (list (butlast (rest convex-segments))))
           (try-merge [(first convex-segments) (last convex-segments)]))
      convex-segments)
    ))

;;
;; merge convex segments if they don't meet on an edge

(defn merge-and-filter [parts]
  (print parts)
  (if (< (len parts) 2)
    parts
    (map (fn [x y] (do ;; (if (< (dist (last x) (first y)) 0.5) (print "no edge?"(no-edge? x y) "\ndot " (bezier-curvature2 x y)))
                    (if (and (< (dist (last x) (first y)) 0.5)
                             (no-edge? x y)
                             (< (bezier-curvature2 x y) 0))
                       (merge-beziers [x y])
                       [])))
         parts
         (+ (rest parts) [(first parts)]))))

(defn merge-and-filter2 [segs]
  (let [[grouped (groupby path (fn [x] (< (bezier-curvature2 x x)
                                          0)))]
        [parts (map (fn [x] (list (second x))) (remove first grouped))]
        [curved (map (fn [x] (list (second x))) (filter first grouped))]]
    (if (= (len parts) 1)
      parts
      (map (fn [x y] (do ;(if (< (dist (last x) (first y)) 0.5) (print "no edge?"(no-edge? x y) "dot " (bezier-curvature2 x y)))
                         (if (and (< (dist (last x) (first y)) 0.5)
                                  (no-edge? x y)
                                  (< (bezier-curvature2 x y) 0))
                           (merge-beziers [x y])
                           [])))
           parts
           (+ (rest parts) [(first parts)])))))

(defn try-merge2 [parts]
  (let [[last-part (first parts)]
        [l []]]
    (for [part (rest parts)]

      (if (no-edge? part last-part)
        (.append (last l) (last consecutive))
        (.append l [(last consecutive)])))
    l)
  )

(defn merge-consecutive [a b]
  (map (fn [x] (try-merge2 (list (second x))))
       (filter first (groupby path convex?))))


(defn get-merged-convex2 [path]
  (let [[convex-segments (map list (filter convex? path))]]
;    (print (merge-and-filter convex-segments))
 ;   (print (array (remove empty? (merge-and-filter convex-segments))))
    (remove empty? (merge-and-filter convex-segments))

    ))


                                ;array ( [[  5.42957003,  28.13219065],[  8.19975053,  16.10221783],[  2.30515485,  40.4036828 ],[-28.54816996,  31.52851723],[-27.50015975,  21.76728318],[-44.97547529,  35.59957164]])


;;
;; path manipulation routines
;;

;; TODO: put here: merge, split etc

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

(defn split-segments-to [cut-off-attribute t path]
  (array (reduce add (map (fn [x] (split-segment-to cut-off-attribute t x)) path))))
