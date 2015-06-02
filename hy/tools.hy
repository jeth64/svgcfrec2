(import [numpy [diag roll amin sort]]
        [itertools [groupby]]
        [operator [add]])


(defn last [coll] ;; in hy since 0.11.0
  (let [[l (list coll)]]
    (if (empty? l) None
        (nth l (dec (len l))))))

(defn zipmap [keys vals] (dict (zip keys vals)))

(defn partition [n step coll] ; clash with numpy
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

(defn concat [&rest colls] ; deprecated, use itertools chain method instead or numpy concatenate
  (for [coll colls]
    (cond [(= (type coll) dict)
           (for [x (.items coll)] (yield x))]
          [True
           (for [x coll] (yield x))])))

(defn catmap [f &rest colls]
  "Do map and concatenate lists returned by map"
  (list (apply concat (list (apply map (cons f colls))))))

(defn case [e &rest clauses]
  (setv res (if (odd? (len clauses)) (last clauses) None))
  (for [clause (partition 2 2 clauses)]
    (when (= e (first clause))
      (do (setv res (second clause))
          (break))))
  res)

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

(defn list-add [l1 l2] (list (map add l1 l2)))

(defn dict-min [coll &optional [f identity]]
  (setv (, mk mv) (first (.iteritems coll)))
  (for [(, k v) (.iteritems coll)]
    (when (< (f v) (f mv)) (do (setv mv v) (setv mk k))))
  mk)
