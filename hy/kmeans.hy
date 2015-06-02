(import [numpy [dot transpose ones argmin all array mean arctan ceil amax amin squeeze reshape cos
                arctan rad2deg]]
        [numpy.linalg [inv norm LinAlgError]]
        [geometry [cross-line]]
        [tools [*]])


(defn least-squares [X y]
  "Returns None if (XTX) singular"
  (try (dot (inv (dot (transpose X) X)) (dot (transpose X) y))
       (catch [e LinAlgError]
         (print (len X))
         None)))


(defn ind-k-means-rec [points edges inds]
  (if (> (len inds) 1)
    (least-squares (array (list (replace points inds))) (ones [(len inds) 1]))
    (least-squares (cross-line (get edges (first inds))) (ones [2 1]))))


(defn ind-k-means [P edges k maxit] ;; for lines
  (setv n (len P))
  (print P)
  (setv groups (list (map (fn [i] (list (range (int (ceil (* i (/ n k))))
                                               (int (ceil (* (inc i) (/ n k)))))))
                          (range k))))

  (setv betas (reshape (array (list (map (fn [g] (ind-k-means-rec P edges g))
                                         groups)))
                       [(len groups) 2]))
  (for [t (range maxit)]
    (setv old-betas betas)
    (setv difference (abs (- (dot P (transpose betas)) 1)))
    (setv labels (argmin difference 1))
    (setv groups (.values (group-by labels (range n))))
   ; (unless (= k (len groups)) (print "Warning: number of groups < k"))
    (setv betas (reshape (array (list (map (fn [g] (ind-k-means-rec P edges g))
                                           groups))) [(len groups) 2]))
    (when (all (= old-betas betas))
      (break)));was wenn kein break?
  (, betas groups (max (amin difference 1))))


(defn ind-get-thetas [edges maxerror maxit]
  (setv points (mean (array (list edges)) 1))
  (setv e 10000)
  (setv k 1)
  (while (> e maxerror) ;; try for different k until distance to lines small enough
    (setv (, betas groups e) (ind-k-means points edges k maxit))
    ;(print k (shape betas) betas)
    (setv k (inc k)))
  (, betas groups e))


(defn k-means [points k maxit] ;; for points in hough-space
  (setv n (len (list points)))
  (setv groups (list (map (fn [i] (list (range (int (ceil (/ (* i n) k)))
                                               (int (ceil (/ (* (inc i) n) k))))))
                          (range k))))
  (setv centroids (list (map (fn [g] (mean (list (replace points g)) 0))
                             groups)))
  (for [t (range maxit)]
    (setv old-centroids centroids)
    (setv difference (list (map (fn [c] (norm (- points c) 2 1)) centroids)))
    (setv labels (argmin difference 0))
    (setv groups (.values (group-by labels (range n))))
    (setv centroids (list (map (fn [g] (mean (list (replace points g)) 0))
                               groups)))
    (when (all (= (array old-centroids) (array centroids)))
      (break)))
  (, centroids groups (max (amin difference 1))))

(defn get-thetas [hough-pts maxerror maxit] ;; for hough-space: not working yet
  (setv pts (array (list hough-pts)))
  (setv k 1)
  (setv e (amax (norm (- pts (mean pts 0)))))
  (while (> e maxerror) ;; try for different k until distance to lines small enough
    (setv k (inc k))
    (setv (, centroids groups e) (k-means pts k maxit)))
  centroids)

; new try...

(defn get-differences [points line-params]
  (let [[rs (list (map (fn [p] (/ (cos (arctan (/ (second p) (first p))))
                                  (first p)))
                       line-params))]]
    (abs (* (- (dot points (transpose line-params)) 1) rs))))

(defn k-lines-rec [points edges inds]; angle node-dirs instead of edges?
  (if (> (len inds) 1)
    (least-squares (array (list (replace points inds))) (ones [(len inds) 1]))
    (least-squares (cross-line (get edges (first inds))) (ones [2 1]))))

(defn k-lines [P edges k maxit]
  (setv n (len P))
  (setv groups (list (map (fn [i] (list (range (int (ceil (* i (/ n k))))
                                               (int (ceil (* (inc i) (/ n k)))))))
                          (range k))))

  (setv parameters (reshape (array (list (map (fn [g] (k-lines-rec P edges g))
                                         groups)))
                            [(len groups) 2]))
  ;;(when (= (len P) 15) (do (print "pts" P) (print "params " parameters)))
  (setv difference (get-differences P parameters))

  (for [t (range maxit)]
    (setv old-parameters parameters)
    (setv labels (argmin difference 1))
    (setv groups (.values (group-by labels (range n))))
   ; (unless (= k (len groups)) (print "Warning: number of groups < k"))
    (setv parameters (reshape (array (list (map (fn [g] (k-lines-rec P edges g))
                                                groups)))
                              [(len groups) 2]))
    (setv difference (get-differences P parameters))

    (when (all (= old-parameters parameters))
      (break)));was wenn kein break?

  (setv e (max (amin difference 1)))
  (, parameters groups e))


(defn get-lines [edges maxerror maxit]
  (setv points (mean (array (list edges)) 1))
   ;;bei 3) 8 + 15 zu 3,5 und 3,12
  (setv e 10000)
  (setv k 1)
  (while (> e maxerror) ;; try for different k until distance to lines small enough
    (setv (, parameters groups e) (k-lines points edges k maxit))
    (setv k (inc k)))
  (setv alphas (list (map (fn [p] (- 90 (rad2deg (arctan (/ (second p) (first p))))))
                          parameters)))
  (, alphas groups e))
