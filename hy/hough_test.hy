(import [hough [*]]
        [numpy [shape array]])

(defn test_get-circle []
  (assert (= [1 2] (get-circle [1 2] 1 8 0)))
  (assert (= [17 2] (shape (array (get-circle [1 2] 2 8 2))))))

(defn test_hough-transform []
  (assert (= [340 2] (shape (hough-transform [[0 0] [1 1]]
                                             (repeat (range 0 10))
                                             {"radius" 1 "nphi" 8 "nrad" 2})))))
