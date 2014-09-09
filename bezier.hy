(import [numpy [*]]
        [scipy.misc [comb]]
        [lisp-tools [*]])


(defn dist [a b] (linalg.norm (- a b)))

(defn pt-coeffs [ts]
  (array (list-comp (* (comb 3 i) (** ts i) (** (- 1 ts) (- 3 i)))
                    [i (range 4)])))

(defn eval-cubics [ts control-points] ; test
  (dot (transpose (pt-coeffs ts)) (array control-points) ))

;;(eval-cubics (array [0 0.5 1]) [[0 1] [2 3] [4 5] [6 7]])

(defn bezier-fit [ordered-pts]
  (let [[path-lengths (+ [0] (reductions add (map dist
                                                 (butlast ordered-pts)
                                                 (rest ordered-pts))))]
        [ts (array (map (fn [x] (/ x (last path-lengths))) path-lengths))]
        [X (transpose (pt-coeffs ts))]
        [control-points (dot (linalg.inv (dot (transpose X) X))
                             (dot (transpose X) ordered-pts))]]
    control-points))

(defn merge-beziers [cubic-beziers]
  (bezier-fit (vstack (map (fn [Cs] (eval-cubics (linspace 0 1 4 True) Cs))
                           cubic-beziers))))

(merge-beziers [[[0 1] [7 8] [9 10] [6 7]] [[6 7] [12 13] [14 15] [11 12]]])

;;(around (bezier-fit (array [[0 1] [7 8] [9 10] [6 7]])))

;;(setv v (array [[0 1] [2 3] [4 5] [6 7]]))


;; Degree reduction: if needed needs to be corrected

(defn weighted-sum [points weights]
  (sum (array (map mul points weights))) 0)

(defn approx-points [points n]
  (rest (reductions (fn [P [i Q]]
                      (weighted-sum [Q P]
                                    [(/ n (- n i)) (/ (- i) (- n i))]))
                    [0 0]
                    (transpose [(range n) points]))))

(defn factor [i n]
  (* (Math/pow 2 (- 1 (* 2 n)))
     (reduce (fn [x y] (+ x (choose (* 2 n) (* 2 y))))
             0 (range (inc i)))))


(defn reduce-degree [control-points]
  (let [[n (dec (count control-points))]
        [Pr (approx-points control-points n)]
        [Pl (reverse (approx-points (reverse control-points) n))]]
    (map (fn [x y z] (let [[lambda (factor z n)]]
                       (weighted-sum [x y] [(- 1 lambda) lambda])))
         Pr Pl (range n))))
