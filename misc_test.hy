(import [numpy [*]]
        [misc [*]])

(defn test_make-set []
  (assert (= (set [(, 1 2) (, 3 4)])
             (make-set [[1 2] [3 4] [2 1]]))))

(defn test_ext-normalize []
  (assert (= [1 0] (ext-normalize [1 0])))
  (assert (= [-1 0] (ext-normalize [-1 0]))))

(defn test_angle []
  (assert (= 0 (angle [1 0])))
  (assert (= 45 (angle [1 1])))
  (assert (= 90 (angle [0 1])))
  (assert (= 135 (angle [-1 1])))
  (assert (= 0 (angle [-1 0])))
  (assert (= 45 (angle [-1 -1]))))
