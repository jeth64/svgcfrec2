(import [hy.contrib [loop]])


(defn last [coll] (nth coll (dec (len coll))))

(defn butlast [coll] (list (take (dec (len coll)) coll)))

(defn zipmap [keys vals] (dict (zip keys vals)))

(defn partition [n step coll]
  (setv remaining coll)
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
