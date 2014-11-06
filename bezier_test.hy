(import [numpy [*]]
        [bezier [*]]
        [lisp-tools [*]])


(defn test_pt-coeffs []
  (assert (all (= [[1 0.125 0] [0 0.375 0] [0 0.375 0] [0 0.125 1]]
                  (around (pt-coeffs [0 0.5 1]) 3))))
  (assert (all (= [[1 0.125 0] [0 0.375 0] [0 0.375 0] [0 0.125 1]]
                  (around (pt-coeffs (map identity [0 0.5 1])) 3)))))

(defn test_eval-cubics []
  (assert (all (= [[0 1] [3 4] [6 7]]
                  (around (eval-cubics (array [0 0.5 1]) [[0 1] [2 3] [4 5] [6 7]]))))))

(defn test_bezier-fit []
  (assert (= [[0 1] [4 3] [17 18] [6 7]]
             (around (bezier-fit [[0 1] [7 8] [9 10] [6 7]]))))
  (for [points (take 1000 (repeatedly (fn [] (around (random.rand 4 2) 3))))]
    (let [[control-points (around (bezier-fit points) 3)]]
      (assert (all (= (first control-points) (first points))))
      (assert (all (= (last control-points) (last points)))))))

(defn test_merge-beziers []
  (assert (all (= [[ 0.158  1.158] [ 5.022  6.022] [ 9.191 10.191] [11.906  12.906]]
                  (around (merge-beziers [[[0 1] [7 8] [9 10] [6 7]]
                                          [[6 7] [12 13] [14 15] [11 12]]]) 3))))))

(defn test_de-casteljau []
  (assert (any (= [6.125 7]
                  (de-casteljau 0.5 [[0 0] [7 8] [8 9] [4 5]]))))
  (for [control-points (take 10 (repeatedly (fn [] (random.rand 4 2))))]
      (assert (all (= (de-casteljau 0 control-points)
                      (first control-points))))
      (assert (all (= (de-casteljau 1 control-points)
                      (last control-points))))))

(defn test_elevate-bezier []
  (assert (all (= (elevate-bezier (array [[1 2] [2 3] [3 4]]))
                   (elevate-bezier [[1 2] [2 3] [3 4]]))))
  (assert (all (= (around (elevate-bezier [[1 2] [2 3] [3 4]]) 2)
                   [[1.0 2.0] [1.67  2.67] [2.33 3.33] [3.0 4.0]]))))
