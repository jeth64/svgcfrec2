(import [numpy [*]]
        [time [clock]])

(defn compare-performance [fn1 fn2 arg-list reps]
  "Example: (compare-performance get-tracedists get-tracedists1 [(random.rand 4)] 10000)"
  (do (setv start-time1 (clock))
      (setv x (list (take reps (repeatedly (fn [] (apply fn1 arg-list))))))
      (setv time1 (- (clock) start-time1))
      (setv start-time2 (clock))
      (setv x (list (take reps (repeatedly (fn [] (apply fn2 arg-list))))))
      (setv time2 (- (clock) start-time2))
      (setv diff (abs (- time1 time2)))
      (print "\nTesting functions with" reps "repetitions...\n")
      (print "  " fn1 "takes" time1 "seconds\n")
      (print "  " fn2 "takes" time2 "seconds\n")
      (print "Conclusion: \n")
      (cond [(= time1 time2) (print "Similar performance")]
            [(< time1 time2) (print "  " fn1 "performs by" diff "s better (with this input)")]
            [(> time1 time2) (print "  " fn2 "performs by" diff "s better (with this input)")])))
