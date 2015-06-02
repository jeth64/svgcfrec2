(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [itertools [combinations]]
        [numpy [*]]
        [operator [add sub]]
        [scipy.sparse.csgraph [connected_components]]
)

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn get-bezier [file]
  (let [[a (apply fromfile [file] {"sep" " "})]
        [n (/ (len a) 8)]]
    (reshape a [n 4 2])))

(defn save-bezier [path file]
  (.tofile path file " "))

(def exclude
  {"test1-path.txt" [[9 11]]
   "test2-path.txt" [[19 21]]
   "test3-path.txt" [[7 9] [16 18] [101 103]]
   "test4-path.txt" [[62 64] [115 117]]})

(def not-exclude
  {"test1-path.txt" []
   "test2-path.txt" []
   "test3-path.txt" []
   "test4-path.txt" [[128 130]]})

(defn get-all-edges [node-count] ;; itertools object
  (combinations (range node-count) 2))


(defn get-simple-curvature [nodes]
  "Returns [0-2-curve 1-3-curve ... (n-3)-(n-3)-curve (n-2)-0-curve (n-1)-1-curve"
  (list (map (fn [x y z] (cross (- x y) (- x z)))
             nodes
             (roll nodes -1)
             (roll nodes -2))))

(defn get-curvature [curves bounds] ;; only one way -> problem?
  (let [[gap (dec (abs (apply sub bounds)))]
        [start (min bounds)]]
    ;(assert (> gap 0))
    (if (> gap 0)
      (reduce add (take gap (drop start curves)))
      0
      )))



(let [[e (get (.items exclude) 3)]
      [path (get-bezier "testfiles/test4-path.txt")]
      [nodes (list (map first path))]
     ; [a (print nodes)]
      [edges (get-all-edges (len path))]
      [curves (get-simple-curvature nodes)]
      [fedges (set (list (filter (fn [pair] (neg? (get-curvature curves pair)))
                                 edges)))
       ]
      ]
  (print (in fedges (, 62 64)))
  (print (in fedges (, 115 117)))
  (print (in fedges (, 128 130)))
  )
