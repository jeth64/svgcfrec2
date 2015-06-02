(import [modulo [*]])

(defn test_mod-dist []
  (assert (= 2 (mod-dist 180 179 1)))
  (assert (= 1 (mod-dist 180 2 3))))

(defn test_mod-mean []
  (assert (= 138 (mod-mean 180 2 94)))
  (assert (= 91 (mod-mean 180 90 92)))
  (assert (= 1 (mod-mean 180 4 178)))
  (assert (= 48 (mod-mean 180 4 92))))
