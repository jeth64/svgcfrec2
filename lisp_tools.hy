(import [numpy [around]]
        [itertools [groupby]])


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
  (let [[l (if (none? init) coll (cons init coll))]
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

(defn group-by [labels coll] ;; not like clojure
 (dict (map (fn [x] [(first x) (map second (second x))])
            (groupby (sorted (zip labels coll) None first) first))))

(defn partition-by [cmp coll] ;; not like clojure
  (setv prev (first coll))
  (setv l [[prev]])
  (for [item (rest coll)]
    (if (cmp prev item)
      (.append l [item])
      (.append (last l) item))
    (setv prev item))
  l)
