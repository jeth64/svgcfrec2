(import [tools [*]]
        [numpy [array roll diag amin]] ;; TODO: check behaviour of roll
        [scipy.spatial.distance [pdist cdist squareform]])


(defn get-pdists [vector-list]
  (squareform (pdist (array vector-list))))

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

(defn get-tracedists [dists] ;; check roll
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
