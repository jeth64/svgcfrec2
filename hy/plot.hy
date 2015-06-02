(import [numpy [*]]
        [tools [*]]
        [matplotlib.pyplot [*]]
        )

(defn show-scatter [points file]
  (figure)
  (title (+ "Hough transform of " file))
  (ylabel "Distance from origin (r)")
  (xlabel "Angle (theta in °)")
  (plot (first pts) (second pts) "ro")
  (show))

(defn show-hough [points file]
  (let [[pts (transpose (list (concat (map (fn [p] (, (+ (first p) 180)
                                                      (- (second p))))
                                           points)
                                      points)))]]
    (figure)
    (title (+ "Hough transform of " file))
    (ylabel "Distance from origin (r)")
    (xlabel "Angle (theta in °)")
    (hist2d (first pts) (second pts) [180 200])
    (show)))
