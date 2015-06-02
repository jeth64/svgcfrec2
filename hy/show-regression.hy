(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]]
        [operator [add sub]]
        [itertools [combinations]]
        [sklearn [linear_model]]
        [matplotlib.pyplot [*]]
        [pylab [savefig]])


(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 5))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

(defn config [key]
  (case key
        :td 4.0
        :tf 1.3))

(defn show-scatter [x y name file]
  (figure)
  (title name)
  (ylabel "Y")
  (xlabel "X")
  (plot x y "ro")
  (show))

;;main
(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv filebasename (.join "" (list (concat (take (- (len infile) 4) infile)))))
      (setv outfile (+ filebasename "-hough.png"))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths (get-paths root namespace))
      (for [path paths]
        (let [[new-path (prepare path)]
              [end-points (get-end-points new-path)]
              [mid-points (cut-end-points new-path)]

              [nodes (list (map first new-path))]

              [n (len nodes)]

              [td (config :td)]
              [tf (config :td)]
              [dists (squareform (pdist (array nodes)))]
              [trace-dists (get-tracedists (get-node2node dists))]
              [graph (& (> td dists)
                        (< tf (/ trace-dists dists)))]

              [ind-edges (make-set (transpose (nonzero graph)))]
              [edges (list (replace2d nodes ind-edges))]

              [model (linear_model.LinearRegression)]
              [model_ransac (linear_model.RANSACRegressor (linear_model.LinearRegression))]

              [cc (second (connected_components graph False))]

              [conn-comp (list (filter (fn [y] (> (len y) 1))
                                       (list (.values (group-by cc (range n)))))) ]


              [grouped-ind-edges (list (map (fn [x] (.intersection ind-edges
                                                                   (make-set (combinations x 2))))
                                            conn-comp))]


              [grouped-edges (map (fn [g] (list (replace2d nodes g)))
                                  grouped-ind-edges)]

              [grouped-edge-means (map (fn [ge] (mean (array (list ge)) 1))
                                       grouped-edges)]
              ]


          (setv i 0)
          (for [gem grouped-edge-means]
            (do
             (let [[pts (transpose gem)]
                   [x (first pts)]
                   [y (second pts)]
                   [xres [(min x) (max x)]]]
               (model.fit (reshape x [(len x) 1]) (reshape y [(len y) 1]))
               ;;  (apply model_ransac.fit (transpose gem))

               (figure)

               (title (+ "Regression of " filebasename ", Component " (str i)))
               (ylabel "Distance from origin (r)")
               (xlabel "Angle (theta in °)")
               (plot x y "ro")
               (plot xres (.predict model (reshape xres [(len xres) 1])) "b-")

               (savefig (+ filebasename "-regression(comp" (str i) ").png")))
             (setv i (inc i)))))))))
