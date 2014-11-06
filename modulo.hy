(import [numpy [abs min]])


(defn mod-dist [n x y]
  (min [(% (- x y) n) (% (- y x) n)]))

(defn mod-mean [n x y]
  (if (< (abs (- x y)) (/ n 2))
    (/ (+ x y) 2)
    (% (/ (+ n x y) 2) n)))
