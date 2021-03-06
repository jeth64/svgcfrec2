(import [tools [*]]
        [numpy [*]]
        [operator [add]])

(defn test_last []
  (assert (= 3 (last [1 2 3])))
  (assert (= 3 (last (map identity [1 2 3])))))

(defn test_zipmap []
  (assert (= {"a" 1 "b" 2}
             (zipmap ["a" "b"] [1 2] ))))

(defn test_partition []
  (assert (= [[1 2 3] [3 4 5]]
             (partition 3 2 [1 2 3 4 5])))
  (assert (= [[1 2 3] [3 4 5]]
             (partition 3 2 (map identity [1 2 3 4 5])))))

(defn test_reductions []
  (assert (= [1 3 6]
             (reductions add [1 2 3])))
  (assert (all (= (array [[1 1] [3 4] [6 8]])
                  (array (reductions add [(array [2 3]) (array [3 4])] (array [1 1]))))))
  (assert (all (= (array [[1 1] [3 4] [6 8]])
                  (array (reductions add (array [[2 3] [3 4]]) (array [1 1]))))))
  (assert (all (= (array [[1 1] [3 4] [6 8]])
                  (array (reductions add (array [[1 1] [2 3] [3 4]]) )))))
  (assert (= [[1 1] [3 4] [6 8]]
             (reductions (fn [x y] (list (map add x y))) [[1 1] [2 3] [3 4]]))))

(defn test_concat []
  (assert (= (list-comp x [x (concat [1 2 3] (array [4 5]) {"sechs" 6} (set [7]))])
             [1 2 3 4 5 (, "sechs" 6) 7])))

(defn test_case []
  (assert (= 1 (case "a" "b" 2 "a" 1 "c" 3 6)))
  (assert (= 2 (case "b" "b" 2 "a" 1 "c" 3 6)))
  (assert (= 6 (case "t" "b" 2 "a" 1 "c" 3 6)))
  (assert (= None (case "t" "b" 2 "a" 1 "c" 3 6))))


(defn test_partition-by []
  (assert (= (partition-by (fn [x y] (not (= (inc x) y))) [1 2 3 5 6 8 9])
             [[1 2 3] [5 6] [8 9]])))

(defn test_get-tracedists []
  (assert (all (= (get-tracedists [4 1 2 3])
                  (amin [[[0 4 5 7] [6 0 1 3] [5 9 0 2] [3 7 8 0]]
                         [[0 6 5 3] [4 0 9 7] [5 1 0 8] [7 3 2 0]]] 0))))
  (assert (all (= (get-tracedists [4 1 2 3])
                  [[0 4 5 3] [4 0 1 3] [5 1 0 2] [3 3 2 0]]))))
