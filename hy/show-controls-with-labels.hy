(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
        [numpy [*]])

(def testfiles
  (list (map (fn [x] (.join "" ["testfiles/test" (.__str__ x) ".svg"]))
             (range 1 9))))

(defn prepare [path]
  (split-segments-to "pathlength" 3.0 path))

;;main

(for [infile testfiles]
  (do (setv tree (parse infile))
      (setv root (.getroot tree))
      (assert (svg? root))
      (setv outfile (.join "" (list (concat (take (- (len infile) 4) infile) "-labeled-controls.svg"))))
      (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
      (setv paths-list (list (get-paths2 root namespace)))
      (for [paths paths-list]
        (for [path paths]
          (let [[new-path (prepare path)]
                [end-points (get-end-points new-path)]
                [mid-points (cut-end-points new-path)]]
            (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
            (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
            (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
            (for [(, point text) (map (fn [i c]
                                      [(first c) (+ (str i) ": " (str (around (first c) 2)))])
                                  (range (len new-path))
                                  new-path)]
            (add-text root {"fill" "orange" "font-size" "0.3px" } point text)))))
      (.write tree outfile)))
