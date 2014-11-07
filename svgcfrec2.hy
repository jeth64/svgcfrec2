(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [lisp-tools [*]]
        [numpy [*]]
        [itertools [groupby combinations]]
        [itertools :as it]
        [operator [sub mul add]]
        [scipy.spatial [ConvexHull]]
        [scipy.spatial.distance [pdist squareform]]
        [scipy.sparse.csgraph [connected_components]]
        )


(def pt-list1 [
              [[ 11.597503   41.41253 ]  [ 11.36267784   40.20682687] [ 11.17037675   38.0891275 ]
               [ 11.04044378   35.35802187] [ 11.00506661   33.85575109]  [ 10.992723   32.3121  ]
               [ 10.97597583   29.98173047]  [ 10.90697066   28.10743125]  [ 10.49960426   25.3131 ]]

              [[  5.13649935   17.3531125 ]  [  3.1423731   15.2603   ]  [  6.1383231   14.2688   ]
               [  8.59075664   15.53787812]]

              ] )



(def dir-list1 [
               [[ 0.51310871   0.58636347]  [ 0.53599573   0.54316846]  [ 0.52171194   0.52319585]
                [ 0.51120641   0.51163272]  [ 0.50385948   0.50398208]  [ 0.50357717   0.50358024]
                [ 0.51838289   0.51805728]  [ 0.57179061   0.56689917]  [ 0.69281553   0.65715113]]

               [[ 0.72758603   0.70690073]  [ 0.67694864   0.70329069]  [ 0.63778636   0.67386013]
                [ 0.71780539   0.5968457 ]]

               ])
(def box-list1 [
               [[ 10.49960426  25.3131    ]
                [  4.02716805  29.38796891]
                [ 18.06993922  37.33766109]
                [ 11.597503    41.41253   ]]
               [[  5.87398678  12.67137425]
                [  3.1423731   15.2603    ]
                [  8.59075664  15.53787812]
                [  5.85914295  18.12680388]]

               ]
  )


(def pt-list2 [ [[ 16.29699788   5.63179749]
                 [ 17.969003     5.6604    ]
                 [ 20.15513488   5.69297156]
                 [ 22.04906542   5.72109746]
                 [ 25.03508175   5.77397969]
                 [ 25.73419675   7.9879875 ]
                 [ 23.47705659   8.43682656]
                 [ 20.912213     8.7935    ]
                 [ 19.18718872   9.00297324]
                 [ 17.77181191   9.22094844]]
                [[ 7.7427231   8.5195    ]
                 [ 6.3427231   7.37525781]
                 [ 4.9333481   6.2247625 ]
                 [ 2.7427231   4.4404    ]
                 [ 5.0304431   3.25      ]
                 [ 6.63983966  4.45506094]
                 [ 8.78123559  5.1557375 ]]
                [[  5.13649935  17.3531125 ]
                 [  3.1423731   15.2603    ]
                 [  6.1383231   14.2688    ]
                 [  8.59075664  15.53787812]]
                [[ 11.597503    41.41253   ]
                 [ 11.36267784  40.20682687]
                 [ 11.17037675  38.0891275 ]
                 [ 11.04044378  35.35802187]
                 [ 11.00506661  33.85575109]
                 [ 10.992723    32.3121    ]
                 [ 10.97597583  29.98173047]
                 [ 10.90697066  28.10743125]
                 [ 10.49960426  25.3131    ]]



               ])

(def dir-list2
  (map (fn [dirs] (map (fn [x] (rotate x [[0 1]])) dirs))
       [[1.1338232477877881 0.91682177324936731 0.85220383768266572 0.93270546389197762
         1.4191512382371059 164.44746282537938 170.41817705629273 172.57971787236528
         172.16064903816601 168.60210726187927]
        [41.221294460366508 39.242492158068188 39.194837243819549 46.344470521420895
         38.045059269647709 27.471495845328022 12.958231363406625]
        [48.160273895148642 48.669242948794555 25.077931233255043 19.966334530627051]
        [67.799852489568238 81.895161599371377 86.043794518442496 87.963588871677985
         89.096417432154226 89.56505149638366 88.739886439212569 84.798595949566192
         74.18581859736679]

        ])
)
(def box-list2 [[[ 25.78136603   8.55710463]
                 [ 25.47585946   4.87103972]
                 [ 16.60250444   9.31786239]
                 [ 16.29699788   5.63179749]]
                [[ 5.72318816  1.75794861]
                 [ 2.7427231   4.4404    ]
                 [ 9.98907186  6.49776137]
                 [ 7.0086068   9.18021276]]
                [[  2.77715652  16.78971238]
                 [  7.99432468  18.03554793]
                 [  3.5279884   13.64546419]
                 [  8.74515656  14.89129973]]
                [[ 11.597503    41.41253   ]
                 [ 16.75262094  39.0692038 ]
                 [  5.34448632  27.6564262 ]
                 [ 10.49960426  25.3131    ]]
                ])


;; Vorgehen:
;; - path preparation: merge or split segments
;; - find convex segments and calculate angle bisector (representative line for convex segment)
;; - calculate intersections of angle bisectors
;; - check if intersections are within triangle described by outer ends of lines
;; - recreate wedge by degree elevation of quadratic bezier with [line-start isec-avg line-end]

(setv part (map array [[ 11.904849     7.7520879 ] [ 22.912619    9.3858879 ] [ 23.272079  11.243688  ][ 12.456489 10.603088  ]]))


(let [[point (dict (zip [:p1 :p2 :p3 :p4] part))]
      [angle-bisector (/ (array [ (+ (get point :p1) (get point :p4))
                                  (+ (get point :p2) (get point :p3)) ]) 2) ]
      [bisector-dir-vec (apply sub angle-bisector)]
      [bisector-length (linalg.norm bisector-dir-vec)] ;; threshold
      [object-dir-vec (- (get point :p4) (get point :p1))]
      ]
  (print "angle-bisector" angle-bisector)
  (print "object-dir-vec" object-dir-vec)
  [angle-bisector (if (pos? (dot bisector-dir-vec object-dir-vec)) 1 0)])



(setv path (array [[[  0.16338942  11.341288  ]
                     [  0.49957942  10.465188  ]
                     [  0.99771942   8.6651879 ]
                     [  1.2703694    7.3412879 ]]
                    [[  1.2703694    7.3412879 ]
                     [  1.5430294    6.0173879 ]
                     [  1.9996794    3.8091879 ]
                     [  2.2851594    2.4341879 ]]
                    [[  2.2851594    2.4341879 ]
                     [  3.0451994   -1.2265121 ]
                     [  4.7661094   -0.61661206]
                     [  4.7661094    3.3133879 ]]
                    [[  4.7661094    3.3133879 ]
                     [  4.7661094    4.4397879 ]
                     [  4.7661094    5.5661879 ]
                     [  4.7661094    6.6925879 ]]
                    [[  4.7661094    6.6925879 ]
                     [  7.14568927   7.04575457]
                     [  9.52526913   7.39892123]
                     [ 11.904849     7.7520879 ]]
                    [[ 11.904849     7.7520879 ]
                     [ 22.912619     9.3858879 ]
                     [ 23.272079    11.243688  ]
                     [ 12.456489    10.603088  ]]
                    [[ 12.456489    10.603088  ]
                     [  5.3997394   10.185188  ]
                     [  3.5379594   10.365188  ]
                     [  3.0994994   11.507788  ]]
                    [[  3.0994994   11.507788  ]
                     [  2.7984494   12.292288  ]
                     [  1.8771294   12.934188  ]
                     [  1.0521294   12.934188  ]]
                    [[  1.0521294   12.934188  ]
                     [  0.09075942  12.934188  ]
                     [ -0.22838058  12.362188  ]
                     [  0.16338942  11.341288  ]]]))

(map angle-bisector (array v))

(def path (array [[[ 11.597503   41.41253  ]
                   [ 11.264873   40.54563  ]
                   [ 10.992723   36.4505   ]
                   [ 10.992723   32.3121   ]]

                  [[ 10.992723   32.3121   ]
                   [ 10.992723   25.3977   ]
                   [ 10.735523   24.4519   ]
                   [  7.8193731  20.6439   ]]

                  [[  7.8193731  20.6439   ]
                   [  6.0740331  18.3648   ]
                   [  3.9693831  15.9421   ]
                   [  3.1423731  15.2603   ]]

                  [[  3.1423731  15.2603   ]
                   [  0.6236831  13.1838   ]
                   [  2.9223231  12.423    ]
                   [  6.1383231  14.2688   ]]

                  [[  6.1383231  14.2688   ]
                   [ 10.153573   16.5732   ]
                   [ 10.992723   16.4557   ]
                   [ 10.992723   13.5889   ]]

                  [[ 10.992723   13.5889   ]
                   [ 10.992723   12.0806   ]
                   [  9.7757631  10.1824   ]
                   [  7.7427231   8.5195   ]]

                  [[  7.7427231   8.5195   ]
                   [  5.9552231   7.0575   ]
                   [  3.7052231   5.2219   ]
                   [  2.7427231   4.4404   ]]

                  [[  2.7427231   4.4404   ]
                   [  0.9265631   2.9658   ]
                   [  0.3902631   1.       ]
                   [  1.8041631   1.       ]]

                  [[  1.8041631   1.       ]
                   [  2.2504531   1.       ]
                   [  3.7022831   2.0125   ]
                   [  5.0304431   3.25     ]]

                  [[  5.0304431   3.25     ]
                   [  7.1893031   5.2615   ]
                   [  8.5608431   5.517    ]
                   [ 17.969003    5.6604   ]]

                  [[ 17.969003    5.6604   ]
                   [ 30.440803    5.8504   ]
                   [ 29.647913    5.7444   ]
                   [ 28.349253    7.0435   ]]

                  [[ 28.349253    7.0435   ]
                   [ 27.775353    7.6174   ]
                   [ 24.428683    8.4049   ]
                   [ 20.912213    8.7935   ]]

                  [[ 20.912213    8.7935   ]
                   [ 15.873713    9.3503   ]
                   [ 14.444983    9.8833   ]
                   [ 14.171113   11.3084   ]]

                  [[ 14.171113   11.3084   ]
                   [ 13.447953   15.0714   ]
                   [ 14.535043   15.7787   ]
                   [ 19.659103   14.8792   ]]

                  [[ 19.659103   14.8792   ]
                   [ 27.727203   13.4627   ]
                   [ 28.261623   13.4551   ]
                   [ 28.718423   14.75     ]]

                  [[ 28.718423   14.75     ]
                   [ 28.981303   15.4952   ]
                   [ 28.304583   16.       ]
                   [ 27.042813   16.       ]]

                  [[ 27.042813   16.       ]
                   [ 23.835553   16.       ]
                   [ 15.219013   20.7086   ]
                   [ 13.893303   23.1858   ]]

                  [[ 13.893303   23.1858   ]
                   [ 13.117933   24.6346   ]
                   [ 12.919693   28.1745   ]
                   [ 13.306763   33.6594   ]]

                  [[ 13.306763   33.6594   ]
                   [ 13.849003   41.34313  ]
                   [ 13.017763   45.11363  ]
                   [ 11.597503   41.41253  ]]]))

(find-thinnings (split-segments-to "pathlength" 5.0  path))



;;main
(do (setv tree (parse ;"more-strings-BC-0.3-2-1-0.2.svg"
                      "testfiles/test2.svg"
                      ))
    (setv root (.getroot tree))
    (assert (svg? root))
    (setv outfile "out.svg")
    (setv namespace (.join "" (take (- (len root.tag) 3) root.tag)))
    (setv paths (get-paths root namespace))
    (for [path paths]
      (let [
          ;  [lines (map angle-bisector path)]
            [lends (list (map (fn [part] (list (slice part 0 4 3))) path))]
      ;      [c (reduce add (map (fn [x] (try-merge (list (second x))) (filter first (groupby path convex?))) [])]
            ;[t 1.2]
            [t 3.0]
            [new-path (array (list (reduce add (list (map (fn [x] (split-segment-to "pathlength" t x))
                                                          path)) [])))]

            [n (len new-path)]
            [end-points (choose-pts new-path (range 0 n) [(int 0)])]
            [mid-points (choose-pts new-path (range 0 n) [1 2])]

      ;;     [c (reduce add (map (fn [x] (try-merge (list (second x)))) (filter first (groupby new-path convex?))) [])]
            [c (get-merged-convex2 new-path)]
         ;   [a (print c (shape c))blueblue]
            [lines (list (map flipud (list (map angle-bisector c))))]
            [thinnings (find-thinnings new-path 4.0 3)]
           ; [lines (vstack [thinnings (map flipud thinnings) (map flipud (map angle-bisector c))])]

            [thinnings (find-thinnings new-path 4.0 3)]

            ;; zum testen
            [data (find-thinnings2 new-path 4.0 3)]
            [pt-list (list (map first data))]
        ;    [a (print "pt-list" pt-list)]
            [dir-list (list (map (fn [dirs] (squeeze (array (list (map (fn [x] (rotate (radians x) [[0 1]])) dirs)))))
                                 (map second data)))]
       ;     [a (print "dir-list" dir-list)]
            [box-list (list (map last data))]

      ;      [a (print "box-list" box-list)]

            ;[a (print "b" bi-dir-thinnings)]
            [lines2 (list (map (fn [x] (stretch x 7)) lines))]
      ;      [d (print "l" lines "l2" lines2)]


        ;;    [new-path-ends (array (map first new-path))]
         ;   [b (print "new-path" (array (map first new-path)))]
         ;;   [ch (ConvexHull (array (map first new-path)))]
      ;      [y (print "y" (get new-path-ends (unique ch.simplices)))]
       ;     [a (print (unique ch.simplices))]
                                ;    [x (print ch.ndim ch.nsimplex )]

            ;; new part:
            [vertices (list (map first lines))]

            [idx-pairs (combinations (range (len lines)) 2)]
            [line-pairs (combinations lines 2)]

            ;; the order should be the same:
            [isecs (dict (zip idx-pairs
                              (list (map isec line-pairs))))]

            [idx-triples (list (combinations (range (len lines)) 3))]
            [line-triples (list (combinations lines 3))]
            [pair-triples (list (map (fn [x] (list (reversed (list (combinations x 2)))))
                                     idx-triples))]

            [wedges2 (list (remove none? (map (fn [x y] (get-wedge x y isecs vertices))
                                              pair-triples idx-triples)))]
            [wedges (list (map first wedges2))]
            [centers (list (map second wedges2))]
     ;       [a (print wedges)]
      ;      [b (for [wedge (remove none? wedges)] (print "w" wedge "\n" ))]
      ;      [a (print centers)]

            ]
        (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
     ;   (add-path root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} c)
      ;  (for [wedge wedges] (add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.3"} wedge))
        (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
        (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
       ; (add-path root {"fill" "none" "stroke" "red" "stroke-width" "0.2"} path)
       ; (for [point (map first new-path)] (add-circle root {"fill" "orange" "r" "0.2"} point))
      ;  (for [point (get new-path-ends (unique ch.simplices))] (add-circle root {"fill" "black" "r" "0.2"} point))
      ;;  (for [point (choose-pts new-path (range 0 (len new-path)) [1 2])] (add-circle root {"fill" "black" "r" "0.2"} point))
       ; (for [line lends] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.3"} line))
        ;;(for [line lines] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
      ;  (for [line lines2] (add-line root {"fill" "none" "stroke" "red" "stroke-width" "0.3"} line))
        (for [line thinnings] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.3"} line))
       ; (for [point (first pt-list1)] (add-circle root {"fill" "red" "r" "0.2"} point))
        (for [rec box-list] (add-polygon root {"fill" "none" "stroke" "red" "stroke-width" "0.1"} rec))
            ;; (for [line (array (zip (reduce add pt-list) (list (map (fn [x y] (add (array x) (array y))) (reduce add pt-list) (reduce add dir-list)))))]

       ; (print "pt-list" pt-list)
       ; (print "dirs" dir-list)
       ; (print "lines" (+ (vstack pt-list) (vstack dir-list)))
        (for [line (list (zip (vstack pt-list)
                              (+ (vstack pt-list) (vstack dir-list))))]
      ;    [a (print "line"line)]
          (add-line root {"fill" "none" "stroke" "magenta" "stroke-width" "0.2"} line))
                                ;   (for [point (.values isecs)] (add-circle root {"fill" "black" "r" "0.2"} point))
                                ;   (for [point centers] (add-circle root {"fill" "orange" "r" "0.2"} point))
              ;;(for [line ] (add-line root {"fill" "none" "stroke" "black" "stroke-width" "path0.3"} line))
        ))
    (.write tree outfile)
    None
    )
