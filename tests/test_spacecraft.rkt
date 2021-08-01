#lang racket/base

(require racket/flonum math/flonum)
(require "spacecraft.dau")

(define (runge-kutta-4 x0 get-deriv delta)
  (define k1 (flvector-scale (get-deriv x0) delta))
  (define k2 (flvector-scale (get-deriv (flvector+ x0 (flvector-scale k1 0.5))) delta))
  (define k3 (flvector-scale (get-deriv (flvector+ x0 (flvector-scale k2 0.5))) delta))
  (define k4 (flvector-scale (get-deriv (flvector+ x0 k3)) delta))

  (for/flvector
      ((x (in-flvector x0))
       (a (in-flvector k1))
       (b (in-flvector k2))
       (c (in-flvector k3))
       (d (in-flvector k4)))
    (+ x (/ a 6.0) (/ b 3.0) (/ c 3.0) (/ d 6.0))))

(define (init-state pos vel)
  (let ((state (make-flvector 48)))
    ;; initial position and velocity
    (for ((i (in-range 3)))
      (flvector-set! state i (flvector-ref pos i))
      (flvector-set! state (+ 3 i) (flvector-ref vel i)))
    ;; identity matrix of partials w.r.t initial state
    (for ((i (in-range 6)))
      (flvector-set! state (+ 6 (* i 7)) 1.0))
    ;; partials w.r.t. GM
    (for ((i (in-range 6)))
      (flvector-set! state (+ i 42) 0.0))
    state))

(let* ((pos0 (flvector 1.989508835731483E+00
                       1.241709323225251E+00
                       5.080233571221722E-01))
       (vel0 (flvector 3.060543697816455E-03
                       1.176041674085188E-02
                       4.933250034401462E-03))
       (state (init-state pos0 vel0))
       (pos-delta (flvector 0.0001 0.0002 0.0003))
       (vel-delta (flvector 0.000001 0.000002 0.000003))
       (state-mod-pos (init-state (flvector+ pos0 pos-delta)
                                  vel0))
       (state-mod-vel (init-state pos0
                                  (flvector+ vel0 vel-delta)))
       (gm 0.0002959122082898587)
       (gm-delta 0.0000001))
  (printf "initial state = ~a\n" state)
  (printf (format "xdot = ~a\n" (xdot state gm)))
  (let-values
      (((state-final state-mod-pos-final state-mod-vel-final state-mod-gm-final)
        (let loop ((i 0)
                   (state state)
                   (state-mod-pos state-mod-pos)
                   (state-mod-vel state-mod-vel)
                   (state-mod-gm state))
          (if (= i 1000)
              (values state state-mod-pos state-mod-vel state-mod-gm)
              (loop (add1 i)
                (runge-kutta-4 state (lambda (x) (xdot x gm)) 0.25)
                (runge-kutta-4 state-mod-pos (lambda (x) (xdot x gm)) 0.25)
                (runge-kutta-4 state-mod-vel (lambda (x) (xdot x gm)) 0.25)
                (runge-kutta-4 state-mod-gm  (lambda (x) (xdot x (fl+ gm gm-delta))) 0.25))))))
    (printf "state: ~a\n" state)
    (printf "difference from pos: ~a\n"
            (flvector- (flvector-copy state-mod-pos-final 0 3)
                       (flvector-copy state-final 0 3)))
    (printf "exp. diff. from pos: ~a\n"
            (for/flvector ((i 3))
              (for/sum ((a pos-delta)
                        (b (flvector-copy state-final (+ 6 (* i 6)) (+ 9 (* i 6)))))
                (fl* a b))))
    (printf "difference from vel: ~a\n"
            (flvector- (flvector-copy state-mod-vel-final 0 3)
                       (flvector-copy state-final 0 3)))
    (printf "exp. diff. from vel: ~a\n"
            (for/flvector ((i 3))
              (for/sum ((a vel-delta)
                        (b (flvector-copy state-final (+ 9 (* i 6)) (+ 12 (* i 6)))))
                (fl* a b))))
    (printf "difference from gm: ~a\n"
            (flvector- (flvector-copy state-mod-gm-final 0 3)
                       (flvector-copy state-final 0 3)))
    (printf "exp. diff. from gm: ~a\n"
            (flvector-scale (flvector-copy state-final (+ 6 36) (+ 9 36))
                            gm-delta))))
