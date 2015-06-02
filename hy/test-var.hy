(import [var [*]]
        [numpy.random [rand]] ;; test for clashes: random/numpy
        [numpy [round array]])

(defn test_angle []
  (assert (= 0 (round (angle (sub (array [2 0]) [2 5])))))
  (assert (= 90 (round (angle (sub (array [0 2]) [5 2])))))
  (assert (= 135 (round (angle (sub (array [2 2]) [5 5]))))))

(defn test_dirvec []
  (for [x (range 100)]
    (let [[v (ext-normalize (rand 2))]]
      (assert (all (= (round v 5) (round (dirvec (angle v)) 5)))))))

(defn test_tiling []
  (let [[p (array (list (map identity (rand 5 2))))]
        [bins (get-tiling-bins p [2 2] [0.1 0.1])]
        [tiledp (array (tiling p bins))]]
    (assert (and (all (< tiledp 3))
                 (all (> tiledp -1))))))
