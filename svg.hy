(import [numpy [*]]
        [bezier [elevate-bezier]]
        [sys [exit]]
        [xml.etree.ElementTree [parse SubElement]]
        [itertools [imap]]
        [re [finditer split match]]
        [operator [add sub]]
        [warnings [warn]]
        [tools [*]])

(defn list-add [l1 l2] (list (map add l1 l2)))


(defn lines2beziers [point-pair-list]
  "Calculates two further points between start and end of line"
  (array (list (map (fn [pair] (column_stack (list (map (fn [x y] (linspace x y 4 True))
                                                        (first pair)  (second pair)))))
                    point-pair-list))))


(defn ptlist2matrix [points offset slicing absolute?]
; "Interpretation of the command point list
;   'slicing' determines the number of points for a single command
;   'offset' determines the initial brush position "
  (if (empty? (list points))
    []
    (let [[pts (partition slicing slicing points)]
          [first-pts (list (butlast (cons offset (list (map last pts)))))]
          [first-abs-pts (array (if absolute? first-pts (reductions list-add first-pts)))]
          [abs-pts (if absolute? (array pts) (array (list (zipwith add first-abs-pts (array pts)))))]]
      (hstack (, (reshape first-abs-pts (, (len first-pts) 1 2))
                 abs-pts)))))

(defn add-first-ctrl-pts [points pt-matrix]
  "Calculates first control points for 'smoothcurveto' command"
  (list (let [[first-ctrl-pts (map (fn [row] (+ (last row) (apply sub (take (int 2) (.reverse row)))))
                                   (cons (last points) pt-matrix))]]
          (zipwith (fn [l x] (.insert l (int 1) x)) pt-matrix first-ctrl-pts))))

(defn single-object-path2matrix [str] ;; archto command missing
  "For each command part calculates 4 points representing it"
  (setv points [])
  (for [it (finditer "([MmCcSsLl])([^A-DF-Za-df-z]+)" str)]
    (let [[cmd (.group it (int 1))]
          [nr-str (-> (.group it (int 2)) (.replace  "-" " -") (.replace  "e -" "e-") (.strip))]
          [numbers (map float (remove empty? (split " |," nr-str)))]
          [pts (partition 2 2 numbers)]
          [offset (if (empty? points) (first pts) (last (last points)))]
          [new-points (cond [(= (.upper cmd) "M")
                             (do (.extend points (array [[offset]]))
                                 (lines2beziers (ptlist2matrix (rest pts) offset 1 (.isupper cmd))))]
                            [(= (.upper cmd) "L")
                             (lines2beziers (ptlist2matrix pts offset 1 (.isupper cmd)))]
                            [(= (.upper cmd) "H")
                             (lines2beziers (ptlist2matrix (zip numbers (repeat (second offset)))
                                                           offset 1 (.isupper cmd)))]
                            [(= (.upper cmd) "V")
                             (lines2beziers (ptlist2matrix (zip (repeat (first offset)) numbers)
                                                           offset 1 (.isupper cmd)))]
                            [(= (.upper cmd) "C")
                             (ptlist2matrix pts offset (int 3) (.isupper cmd))]
                            [(= (.upper cmd) "S")
                             (add-first-ctrl-pts points (ptlist2matrix pts offset 2 (.isupper cmd)))]
                            [(= (.upper cmd) "Q")
                             (elevate-bezier (ptlist2matrix pts offset 2 (.isupper cmd)))]
                            [(= (.upper cmd) "T")
                             (elevate-bezier (add-first-ctrl-pts
                                              points (ptlist2matrix pts offset 1 (.isupper cmd))))]
                            [(= (.upper cmd) "Z")
                             (lines2beziers [[offset (first (first points))]])]
                            [True (raise (Exception "Unknown path command"))])]]
      (unless (empty? new-points)
        (.extend points new-points))))
  (array (list (rest points))))


(defn path2matrix [str]
  "Devides path string at 'moveto' commands (M/m) and returns a matrix for each substring"
  (list (imap (fn [m] (single-object-path2matrix (.group m (int 0))))
              (finditer "([Mm][^Mm]+)" str))))

(defn get-paths [node namespace]
  "Returns all paths which are subelements of given node as matrices"
  (setv paths [])
  (for [path (.iter node (+ namespace "path"))]
    (.extend paths (path2matrix (.get path "d"))))
  paths)


(defn matrix2path [matrix]
  "Return string for polybezier curve (input: nx4x2-matrix)"
  (let [[ptstr-matrix
         (list (map (fn [part] (list (map (fn [p] (.join "," (list-comp (.__str__ x) [x p])))
                                          part)))
                    (around matrix 3)))]
        [M (+ "M" (first (first ptstr-matrix)))]
        [Cs (.join " " (map (fn [part] (+ "C" (.join " " (rest part))))
                            ptstr-matrix))]]
    (+ M " " Cs)))

(defn add-path [parent attribs path] ;; side-effects
  (.update attribs {"d" (matrix2path path)})
  (SubElement parent "ns0:path" attribs))

(defn add-line [parent attribs line] ;; side-effects
  (.update attribs (zip ["x1" "y1" "x2" "y2"] (map str (flatten line))))
  (SubElement parent "ns0:line" attribs))

(defn add-circle [parent attribs center]
  (.update attribs (zip ["cx" "cy"] (map str center)))
  (SubElement parent "ns0:circle" attribs))

(defn add-rectangle [parent attribs center width height rotation]
  (let [[offset (- center (/ [width height] 2))]
        [transform (str "rotate("  (str rotation) " " (.join " " (map str center)) ")")]])
  (.update attribs (zip ["x" "y" "width" "height" "transform"]
                        [(first offset) (second offset) width height transform]))
  (SubElement parent "ns0:rect" attribs))

(defn add-polygon [parent attribs vertices]
  (.update attribs {"points" (.join " " (map (fn [x] (.join "," (map str x))) vertices))})
  (SubElement parent "ns0:polygon" attribs))

(defn svg? [root]
  (= (.join "" (drop (- (len root.tag) 3) root.tag)) "svg"))
