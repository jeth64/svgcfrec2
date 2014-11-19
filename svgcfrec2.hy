(import [xml.etree.ElementTree [parse SubElement]]
        [misc [*]]
        [svg [*]]
        [bezier [*]]
        [geometry [*]]
        [tools [*]]
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
(def new-path [[[ 11.597503    41.41253   ]
  [ 11.5143455   41.195805  ]
  [ 11.434968    40.77731562]
  [ 11.36267784  40.20682687]]

 [[ 11.36267784  40.20682687]
  [ 11.29038769  39.63633812]
  [ 11.22518488  38.91385   ]
  [ 11.17037675  38.0891275 ]]

 [[ 11.17037675  38.0891275 ]
  [ 11.11556862  37.264405  ]
  [ 11.07115519  36.33744812]
  [ 11.04044378  35.35802187]]

 [[ 11.04044378  35.35802187]
  [ 11.02508808  34.86830875]
  [ 11.01315788  34.36547828]
  [ 11.00506661  33.85575109]]

 [[ 11.00506661  33.85575109]
  [ 10.99697534  33.34602391]
  [ 10.992723    32.8294    ]
  [ 10.992723    32.3121    ]]

 [[ 10.992723    32.3121    ]
  [ 10.992723    31.4478    ]
  [ 10.98870425  30.67675938]
  [ 10.97597583  29.98173047]]

 [[ 10.97597583  29.98173047]
  [ 10.96324741  29.28670156]
  [ 10.94180933  28.66768438]
  [ 10.90697066  28.10743125]]

 [[ 10.90697066  28.10743125]
  [ 10.83729332  26.986925  ]
  [ 10.71401363  26.101475  ]
  [ 10.49960426  25.3131    ]]

 [[ 10.49960426  25.3131    ]
  [ 10.28519489  24.524725  ]
  [  9.97965584  23.833425  ]
  [  9.54545976  23.10121875]]

 [[  9.54545976  23.10121875]
  [  9.32836172  22.73511562]
  [  9.07909942  22.35878594]
  [  8.79298195  21.95498203]]

 [[  8.79298195  21.95498203]
  [  8.50686448  21.55117812]
  [  8.18389184  21.1199    ]
  [  7.8193731   20.6439    ]]

 [[  7.8193731   20.6439    ]
  [  7.3830381   20.074125  ]
  [  6.92424623  19.495375  ]
  [  6.46857482  18.93709531]]

 [[  6.46857482  18.93709531]
  [  6.01290341  18.37881563]
  [  5.56035248  17.84100625]
  [  5.13649935  17.3531125 ]]

 [[  5.13649935  17.3531125 ]
  [  4.2887931   16.377325  ]
  [  3.5558781   15.6012    ]
  [  3.1423731   15.2603    ]]

 [[  3.1423731   15.2603    ]
  [  2.5127006   14.741175  ]
  [  2.18411123  14.30428125]
  [  2.09566794  13.96978906]]

 [[  2.09566794  13.96978906]
  [  2.00722466  13.63529688]
  [  2.15892748  13.40320625]
  [  2.48983935  13.2936875 ]]

 [[  2.48983935  13.2936875 ]
  [  2.82075123  13.18416875]
  [  3.33087216  13.19722188]
  [  3.95926513  13.35301719]]

 [[  3.95926513  13.35301719]
  [  4.5876581   13.5088125 ]
  [  5.3343231   13.80735   ]
  [  6.1383231   14.2688    ]]

 [[  6.1383231   14.2688    ]
  [  7.14213558  14.8449    ]
  [  7.94744181  15.26963125]
  [  8.59075664  15.53787812]]

 [[  8.59075664  15.53787812]
  [  9.23407147  15.806125  ]
  [  9.71539489  15.9178875 ]
  [ 10.07124176  15.86805   ]]

 [[ 10.07124176  15.86805   ]
  [ 10.7829355   15.768375  ]
  [ 10.992723    15.0223    ]
  [ 10.992723    13.5889    ]]

 [[ 10.992723    13.5889    ]
  [ 10.992723    12.83475   ]
  [ 10.68848303  11.983125  ]
  [ 10.13011305  11.112175  ]]

 [[ 10.13011305  11.112175  ]
  [  9.85092806  10.6767    ]
  [  9.50821058  10.23639375]
  [  9.10822433   9.801025  ]]

 [[  9.10822433   9.801025  ]
  [  8.70823809   9.36565625]
  [  8.2509831    8.935225  ]
  [  7.7427231    8.5195    ]]

 [[  7.7427231    8.5195    ]
  [  7.2958481    8.154     ]
  [  6.82006685   7.76515   ]
  [  6.3427231    7.37525781]]

 [[  6.3427231    7.37525781]
  [  5.86537935   6.98536563]
  [  5.3864731    6.59443125]
  [  4.9333481    6.2247625 ]]

 [[  4.9333481    6.2247625 ]
  [  4.0270981    5.485425  ]
  [  3.2239731    4.83115   ]
  [  2.7427231    4.4404    ]]

 [[  2.7427231    4.4404    ]
  [  1.8346431    3.7031    ]
  [  1.2465281    2.843     ]
  [  1.0621706    2.167225  ]]

 [[  1.0621706    2.167225  ]
  [  0.8778131    1.49145   ]
  [  1.0972131    1.        ]
  [  1.8041631    1.        ]]

 [[  1.8041631    1.        ]
  [  2.0273081    1.        ]
  [  2.5018381    1.253125  ]
  [  3.08660185   1.6609375 ]]

 [[  3.08660185   1.6609375 ]
  [  3.6713656    2.06875   ]
  [  4.3663631    2.63125   ]
  [  5.0304431    3.25      ]]

 [[  5.0304431    3.25      ]
  [  5.5701581    3.752875  ]
  [  6.0606656    4.146     ]
  [  6.63983966   4.45506094]]

 [[  6.63983966   4.45506094]
  [  7.21901372   4.76412187]
  [  7.88685434   4.98911875]
  [  8.78123559   5.1557375 ]]

 [[  8.78123559   5.1557375 ]
  [  9.22842621   5.23904687]
  [  9.73225199   5.30776172]
  [ 10.30994718   5.36509277]]

 [[ 10.30994718   5.36509277]
  [ 10.88764237   5.42242383]
  [ 11.53920697   5.46837109]
  [ 12.28187525   5.50614531]]

 [[ 12.28187525   5.50614531]
  [ 13.02454352   5.54391953]
  [ 13.85831546   5.5735207 ]
  [ 14.80042534   5.59815957]]

 [[ 14.80042534   5.59815957]
  [ 15.27148027   5.610479  ]
  [ 15.76961969   5.62155786]
  [ 16.29699788   5.63179749]]

 [[ 16.29699788   5.63179749]
  [ 16.82437606   5.64203711]
  [ 17.38099301   5.6514375 ]
  [ 17.969003     5.6604    ]]

 [[ 17.969003     5.6604    ]
  [ 18.7484905    5.672275  ]
  [ 19.4761628    5.68299375]
  [ 20.15513488   5.69297156]]

 [[ 20.15513488   5.69297156]
  [ 20.83410696   5.70294937]
  [ 21.46437882   5.71218623]
  [ 22.04906542   5.72109746]]

 [[ 22.04906542   5.72109746]
  [ 23.21843863   5.73891992]
  [ 24.20547081   5.75543984]
  [ 25.03508175   5.77397969]]

 [[ 25.03508175   5.77397969]
  [ 25.86469269   5.79251953]
  [ 26.53688238   5.8130793 ]
  [ 27.07657058   5.83898145]]

 [[ 27.07657058   5.83898145]
  [ 27.61625878   5.86488359]
  [ 28.0234455    5.89612812]
  [ 28.3230505    5.9360375 ]]

 [[ 28.3230505    5.9360375 ]
  [ 29.5214705    6.095675  ]
  [ 28.998583     6.39395   ]
  [ 28.349253     7.0435    ]]

 [[ 28.349253     7.0435    ]
  [ 28.062303     7.33045   ]
  [ 27.0821605    7.6708    ]
  [ 25.73419675   7.9879875 ]]

 [[ 25.73419675   7.9879875 ]
  [ 25.06021487   8.14658125]
  [ 24.29427769   8.29938438]
  [ 23.47705659   8.43682656]]

 [[ 23.47705659   8.43682656]
  [ 22.6598355    8.57426875]
  [ 21.7913305    8.69635   ]
  [ 20.912213     8.7935    ]]

 [[ 20.912213     8.7935    ]
  [ 20.2824005    8.8631    ]
  [ 19.70899066   8.93232812]
  [ 19.18718872   9.00297324]]

 [[ 19.18718872   9.00297324]
  [ 18.66538679   9.07361836]
  [ 18.19519277   9.14568047]
  [ 17.77181191   9.22094844]]

 [[ 17.77181191   9.22094844]
  [ 16.92505019   9.37148437]
  [ 16.26554113   9.53484375]
  [ 15.75492675   9.7253375 ]]

 [[ 15.75492675   9.7253375 ]
  [ 14.733698    10.106325  ]
  [ 14.308048    10.59585   ]
  [ 14.171113    11.3084    ]]

 [[ 14.171113    11.3084    ]
  [ 13.990323    12.24915   ]
  [ 13.92267363  12.99891875]
  [ 14.00295738  13.58034531]]

 [[ 14.00295738  13.58034531]
  [ 14.08324113  14.16177187]
  [ 14.311458    14.57485625]
  [ 14.7224005   14.8422375 ]]

 [[ 14.7224005   14.8422375 ]
  [ 15.133343    15.10961875]
  [ 15.72701113  15.23129687]
  [ 16.53819738  15.22991094]]

 [[ 16.53819738  15.22991094]
  [ 16.9437905   15.22921797]
  [ 17.40376316  15.19775898]
  [ 17.92246441  15.13836387]]

 [[ 17.92246441  15.13836387]
  [ 18.44116566  15.07896875]
  [ 19.0185955   14.9916375 ]
  [ 19.659103    14.8792    ]]

 [[ 19.659103    14.8792    ]
  [ 20.6676155   14.7021375 ]
  [ 21.55841425  14.54708906]
  [ 22.34606187  14.41384687]]

 [[ 22.34606187  14.41384687]
  [ 23.13370948  14.28060469]
  [ 23.81820597  14.16916875]
  [ 24.41411394  14.07933125]]

 [[ 24.41411394  14.07933125]
  [ 25.60592988  13.89965625]
  [ 26.44339175  13.806375  ]
  [ 27.0430005   13.797825  ]]

 [[ 27.0430005   13.797825  ]
  [ 28.242218    13.780725  ]
  [ 28.490023    14.10255   ]
  [ 28.718423    14.75      ]]

 [[ 28.718423    14.75      ]
  [ 28.981303    15.4952    ]
  [ 28.304583    16.        ]
  [ 27.042813    16.        ]]

 [[ 27.042813    16.        ]
  [ 26.6419055   16.        ]
  [ 26.156478    16.07357187]
  [ 25.6113354   16.20716094]]

 [[ 25.6113354   16.20716094]
  [ 25.0661928   16.34075   ]
  [ 24.46133511  16.53435625]
  [ 23.82156722  16.774425  ]]

 [[ 23.82156722  16.774425  ]
  [ 23.18179933  17.01449375]
  [ 22.50712124  17.301025  ]
  [ 21.82233786  17.62046406]]

 [[ 21.82233786  17.62046406]
  [ 21.13755448  17.93990313]
  [ 20.44266581  18.29225   ]
  [ 19.76247675  18.66395   ]]

 [[ 19.76247675  18.66395   ]
  [ 19.08228769  19.03565   ]
  [ 18.41679823  19.42670313]
  [ 17.79081329  19.82355469]]

 [[ 17.79081329  19.82355469]
  [ 17.16482835  20.22040625]
  [ 16.57834792  20.62305625]
  [ 16.05617691  21.01795   ]]

 [[ 16.05617691  21.01795   ]
  [ 15.53400589  21.41284375]
  [ 15.07614429  21.79998125]
  [ 14.707397    22.16580781]]

 [[ 14.707397    22.16580781]
  [ 14.33864972  22.53163437]
  [ 14.05901675  22.87615   ]
  [ 13.893303    23.1858    ]]

 [[ 13.893303    23.1858    ]
  [ 13.6994605   23.548     ]
  [ 13.54168863  24.04089375]
  [ 13.42011519  24.66219844]]

 [[ 13.42011519  24.66219844]
  [ 13.29854175  25.28350313]
  [ 13.21316675  26.03321875]
  [ 13.164118    26.9090625 ]]

 [[ 13.164118    26.9090625 ]
  [ 13.13959363  27.34698437]
  [ 13.12415081  27.81643828]
  [ 13.11780554  28.31713887]]

 [[ 13.11780554  28.31713887]
  [ 13.11146027  28.81783945]
  [ 13.11421253  29.34978672]
  [ 13.12607831  29.91269531]]

 [[ 13.12607831  29.91269531]
  [ 13.13794409  30.47560391]
  [ 13.15892339  31.06947383]
  [ 13.18903218  31.69401973]]

 [[ 13.18903218  31.69401973]
  [ 13.21914097  32.31856562]
  [ 13.25837925  32.9737875 ]
  [ 13.306763    33.6594    ]]

 [[ 13.306763    33.6594    ]
  [ 13.374543    34.61986625]
  [ 13.42086238  35.51918828]
  [ 13.44725327  36.35041615]]

 [[ 13.44725327  36.35041615]
  [ 13.47364417  37.18164402]
  [ 13.48010659  37.94477773]
  [ 13.46817269  38.63286734]]

 [[ 13.46817269  38.63286734]
  [ 13.45623878  39.32095695]
  [ 13.42590855  39.93400246]
  [ 13.37871413  40.46505393]]

 [[ 13.37871413  40.46505393]
  [ 13.33151972  40.99610539]
  [ 13.26746113  41.44516281]
  [ 13.1880705   41.80527625]]

 [[ 13.1880705   41.80527625]
  [ 13.02928925  42.52550313]
  [ 12.80917988  42.88995406]
  [ 12.53999956  42.84302953]]

 [[ 12.53999956  42.84302953]
  [ 12.27081925  42.796105  ]
  [ 11.952568    42.337805  ]
  [ 11.597503    41.41253   ]]]
  )

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
            [n (len new-path)]

            [end-points (get-end-points new-path)]
            [mid-points (cut-end-points new-path)]

      ;;     [c (reduce add (map (fn [x] (try-merge (list (second x)))) (filter first (groupby new-path convex?))) [])]
            [c (get-merged-convex2 new-path)]
         ;   [a (print c (shape c))]
            [lines (list (map flipud (list (map angle-bisector c))))]
          ;  [thinnings (find-thinnings new-path 4.0 3)]
           ; [lines (vstack [thinnings (map flipud thinnings) (map flipud (map angle-bisector c))])]


            ;; zum testen
            [data (find-thinnings2 new-path 4.0 3)]
            [pt-list (list (map first data))]
        ;    [a (print "pt-list" pt-list)]
            [dir-list (list (map (fn [dirs] (squeeze (array (list (map (fn [x] (rotate (radians x) [[0 1]])) dirs)))))
                                 (map second data)))]
       ;     [a (print "dir-list" dir-list)]
            [box-list (list (map last data))]

         ;   [a (print "new-path" new-path) ]
            [thinnings (find-thinnings-new new-path 6)]

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

            ]
        (add-path root {"fill" "none" "stroke" "yellow" "stroke-width" "0.3"} new-path)
      ;  (for [wedge wedges] (add-path root {"fill" "none" "stroke" "blue" "stroke-width" "0.3"} wedge))
        (for [point end-points] (add-circle root {"fill" "blue" "r" "0.2"} point))
        (for [point mid-points] (add-circle root {"fill" "green" "r" "0.2"} point))
        (for [line thinnings] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.3"} line))
       ; (for [line lines] (add-line root {"fill" "none" "stroke" "orange" "stroke-width" "0.1"} line))
        (for [rec box-list] (add-polygon root {"fill" "none" "stroke" "red" "stroke-width" "0.1"} rec))

        (for [line (list (zip (vstack pt-list)
                              (+ (vstack pt-list) (vstack dir-list))))]
          (add-line root {"fill" "none" "stroke" "magenta" "stroke-width" "0.2"} line))
        ))
    (.write tree outfile)
    None
    )
