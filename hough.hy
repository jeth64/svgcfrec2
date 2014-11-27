(import [numpy [cos sin deg2rad linspace pi transpose]]
        [tools [concat catmap]])

(defn get-circle [midpoint &optional [radius 1] [nphi 8] [nrad 1]]
  "Returns list of points in circle around midpoint"
  ;; n-phi specifies number of angles, that points should be calculated for
  ;; n-rad specifies number of points in one angle direction (midpoint excluded)
  (let [[phis (linspace 0 (* 2 pi) nphi False)]
        [unit-circles (transpose [(cos phis) (sin phis)])]]
    (cons midpoint (catmap (fn [rad] (+ (* unit-circles rad) midpoint))
                           (rest (linspace 0 radius (inc nrad)))))))

(defn full-hough-range [hough-pts]
  "Mirrors points to extend range from 0-180 to 0-360"
  (list (concat (map (fn [p] (, (+ (first p) 180)
                                (- (second p))))
                     hough-pts)
                hough-pts)))

(defn hough-transform [points &optional
                       [theta-list (repeat (range 0 180))]
                       [circle-ops {"radius" 0 "nphi" 8 "nrad" 0}]
                       [r-range [0]]] ;test
  "Performs hough-transform of points"
  ;; theta-list can be used to set certain thetas for each point,
  ;;   e.g. theta-list = (repeat (range 70 90))
  ;; r-range can be used to produce additional points for each theta in r-direction,
  ;;   e.g. r-range = (arange -1 1 0.02)
  ;; for calculating hough-points for points around each given point, give circle-ops,
  ;;   e.g.: {"radius" 1 "n-phi" 8 "n-rad" 2} for 16 additional points in 8 different
  ;;         directions with radius 1 and 0.5
  (catmap (fn [pt thetas]
            (let [[pts (apply get-circle [pt] circle-ops)]]
              (catmap (fn [theta]
                        (catmap (fn [p]
                                  (list (map (fn [offset]
                                               [theta
                                                (+ (* (first p) (cos (deg2rad theta)))
                                                   (* (second p) (sin (deg2rad theta)))
                                                   offset)])
                                             r-range)))
                                pts))
                      (list thetas))))
            points
            theta-list))
