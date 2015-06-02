(import [tools [*]]
        [geometry [*]]
        [bezier [*]]
        [scipy.spatial.distance [cdist]]
        [itertools [combinations]]
        )

(defn reduce-wedges [wedge-list] ;; noch: Dartellung durch 3 Eckpunkte
  (setv wedges (list (map (fn [w] (tuple (map (fn [v] (tuple v)) w)))
                          wedge-list)))
  (setv vertex-dict {})
  (setv valid (set []))
  (setv invalid (set []))
  (setv unchecked (set wedges))

  (for [wedge wedges]
    (for [vertex wedge]
      (if (not (in vertex vertex-dict))
        (assoc vertex-dict vertex (set [wedge]))
        (.add (get vertex-dict vertex) wedge))))
  (while (> (len unchecked) 0)
    (setv min-count-vert (dict-min vertex-dict
                                   (fn [x] (len (.intersection unchecked x)))))
    (setv poss-new-wedges (.intersection unchecked (get vertex-dict min-count-vert)))
    (when (empty? poss-new-wedges)
      (.pop vertex-dict min-count-vert)
      (continue))

    (setv new-wedge (first poss-new-wedges))

    (.remove unchecked new-wedge)
    (.add valid new-wedge)

    (for [vertex new-wedge]
      (for [wedge (get vertex-dict vertex)]
        (when (in wedge unchecked)
          (.remove unchecked wedge)
          (.add invalid wedge)))
      (.pop vertex-dict vertex)))
  valid)

(defn calc-wedge [index-pairs inds isecs lines] ;; TODO: do only degree elevation here nad only for chosen wedges, rename to "reconstruct-wedges
  (let [[wlines (list (replace lines inds))]
        [vertices (list (apply concat wlines))]
        [gisecs (list (map (fn [x] (get isecs x)) index-pairs))]
        ;; determine wedge corners:
        [vs (array (list (map (fn [l] (let [[dists (cdist vertices l)]
                                            [cumdist (sum dists 0)]]
                                        (if (> (first cumdist) (second cumdist))
                                          (first l)
                                          (second l))))
                              wlines)))]]
    (if (points-in-triangle? vs gisecs)
      [(array (list (map (fn [start end] (elevate-bezier [start (mean gisecs 0) end]))
                         vs (roll vs -1 0))))
       gisecs
       vs]
      None)))

(defn get-wedges [lines]
  (let [[idx-pairs (combinations (range (len lines)) 2)]
        [line-pairs (combinations lines 2)]
        ;; the order should be the same:
        [isecs (dict (zip idx-pairs
                          (list (map isec line-pairs))))]
        [idx-triples (list (combinations (range (len lines)) 3))]
        [line-triples (list (combinations lines 3))]
        [pair-triples (list (map (fn [x] (list (reversed (list (combinations x 2)))))
                                 idx-triples))]
        [w (list (map (fn [x y] (calc-wedge x y isecs lines)) pair-triples idx-triples))]]
    (list (remove none? w))))


(defn find-wedges [edges] ;; TODO: implement check-wedges, reconstruct wedges
  (cond [(< (len edges) 3) []]
        [(and (= (len edges) 3) (check-wedge edges)) (reconstruct-wedge edges)]
        [(> (len edges) 3) (get-wedges edges)]
        [True []]))
