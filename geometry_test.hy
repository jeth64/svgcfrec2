(import [numpy [*]]
        [geometry [*]])

(defn test_stretch [line factor]
  (assert (all (= [[1 1] [3 3]]
                  (stretch [[1 1] [2 2]] 2)))))

(defn test_isec [lines]
  (assert (all (= [2 2]
                  (isec [[[1 1] [3 3]] [[1 3] [3 1]]])))))

(defn test_cross-line []
  (assert (all (= [[2 1] [1 2]]
                  (cross-line [[1 1] [2 2]])))))

(defn test_rotation-matrix []
  (assert (all (= [[0 1] [-1 0]]
                  (around (rotation-matrix (/ pi 2)))))))

(defn test_rotate []
  (assert (all (= [[1 0] [0 2]]
                  (around (rotate (/ pi 2) [[0 1] [-2 0]]))))))
