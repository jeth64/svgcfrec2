(import [numpy [diag roll amin sort]]
        [itertools [groupby]]
        [operator [add]])


(defn last [coll]
  (let [[l (list coll)]]
    (nth l (dec (len l)))))

(defn zipmap [keys vals] (dict (zip keys vals)))

(defn partition [n step coll]
  (setv remaining (list coll))
  (setv res [])
  (while (>= (len remaining) n)
    (.append res (list (take n remaining)))
    (setv remaining (list (drop step remaining))))
  res)

(defn reductions [f coll &optional init]
  (let [[l (if (none? init) (list coll) (cons init (list coll)))]
        [res [(first l)]]]
    (for [item (rest l)]
      (.append res (f (last res) item)))
    res))

(defn concat [&rest colls]
  (for [coll colls]
    (cond [(= (type coll) dict)
           (for [x (.items coll)] (yield x))]
          [True
           (for [x coll] (yield x))])))


(defn replace [coll keys]
  "Replaces given keys by items in collection"
  (map (fn [x] (get coll x)) keys))

(defn replace2d [coll key-matrix]
  "Replaces keys in matrix by items in collection"
  (map (fn [key] (list (replace coll key))) key-matrix))

(defn group-by [labels coll] ;; not like clojure
  (dict (list (map (fn [x] [(first x) (list (map second (second x)))])
                   (groupby (sorted (zip labels coll) None first) first)))))

(defn partition-by [cmp coll] ;; not like clojure
  (setv prev (first coll))
  (setv l [[prev]])
  (for [item (rest coll)]
    (if (cmp prev item)
      (.append l [item])
      (.append (last l) item))
    (setv prev item))
  l)

(defn make-set [pairings]
  (set (list (map (fn [x] (tuple (sort x))) pairings))))


;;
;; distance functions
;;

(defn get-node2node [distances]
  "Returns for pairwise distance table: [0-1-dist 1-2-dist ... (n-2)-(n-1)-dist (n-1)-0-dist]"
  (list (concat (diag distances 1) [(first (last distances))])))

(defn get-tracedists1 [dists] ;; version 1, 4 s slower with 100000 repetitions than version 2
  "Returns the sum of node-to-node distances along the shortest path on the hull"
  (let [[N (len dists)]
        [tracedists-right
         (list (map (fn [i] (roll (cons 0 (reductions add (butlast (roll dists (- i)))))
                                  i))
                    (range N)))]
        [tracedists-left
         (list (map (fn [i] (roll (list (cons 0 (reversed
                                                 (reductions add
                                                             (reversed (list (rest (roll dists
                                                                                         (- i)))))))))
                                  i))
                    (range N)))]]
    (amin [tracedists-right tracedists-left] 0)))

(defn get-tracedists [dists]
  "Returns the sum of node-to-node distances along the shortest path on the hull"
  ;; dists: matrix of pairwise distances as returned by '(squareform (pdist v))'
  (list (map (fn [i] (let [[d (roll dists (- i))]]
                       (roll (list (cons 0 (amin [(reductions add (butlast d))
                                                  (list (reversed
                                                         (reductions add
                                                                     (reversed (list (rest d))))))]
                                                 0)))
                             i)))
             (range (len dists)))))
