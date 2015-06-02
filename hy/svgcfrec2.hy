
;;main
(do (setv tree (parse ;"more-strings-BC-0.3-2-1-0.2.svg"
                      "testfiles/test1.svg"
                      ))
    (setv root (.getroot tree))
    (assert (svg? root))
    (setv outfile "out.svg")
    (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
    (setv paths (get-paths root namespace))
    (for [path paths]
      (let [
      ;      [c (reduce add (map (fn [x] (try-merge (list (second x))) (filter first (groupby path convex?))) [])]

            [new-path (split-segments-to "pathlength" 3.0 path)]

            [a (print "new-path" new-path)]
