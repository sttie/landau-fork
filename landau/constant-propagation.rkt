#lang racket/base

(require
   racket/extflonum
   racket/fixnum
   racket/flonum
   racket/function
   racket/list
   racket/vector
   "target-config.rkt"
   "environment.rkt")

(provide rl0.0 rl* rl/ rl+ rl- rl-neg rl-sin 
         rl-cos rl-tan ->rl inexact->rl rl-sqrt 
         rl-sqr rl-expt rl-vector-ref
         list->rl-vector make-rl-vector any-number?)

(define rl0.0
  (if (target-extfloat? TARGET)
    0.0t0
    0.0))

(define rl*
  (if (target-extfloat? TARGET)
    extfl*
    fl*))

(define rl/
  (if (target-extfloat? TARGET)
    extfl/
    fl/))

(define rl+
  (if (target-extfloat? TARGET)
    extfl+
    fl+))

(define rl-
  (if (target-extfloat? TARGET)
    extfl-
    fl-))

(define (rl-neg n)
  (if (target-extfloat? TARGET)
    (extfl- 0.0t0 n)
    (fl- n)))

(define rl-sin
  (if (target-extfloat? TARGET)
    extflsin
    flsin))

(define rl-cos
  (if (target-extfloat? TARGET)
    extflcos
    flcos))

(define rl-tan
  (if (target-extfloat? TARGET)
    extfltan
    fltan))

(define ->rl
  (if (target-extfloat? TARGET)
    ->extfl
    ->fl))

(define inexact->rl
  (if (target-extfloat? TARGET)
    real->extfl
    identity))

(define (rl-sqr x)
  (if (target-extfloat? TARGET)
    (lambda (x) (extfl* x x))
    (lambda (x) (fl* x x))))

(define rl-sqrt
  (if (target-extfloat? TARGET)
    extflsqrt
    flsqrt))

(define rl-expt
  (if (target-extfloat? TARGET)
    extflexpt
    flexpt))

(define rl-vector-ref
  (if (target-extfloat? TARGET)
    extflvector-ref
    flvector-ref))

(define make-rl-vector
  (if (target-extfloat? TARGET)
    make-extflvector
    make-flvector))

(define (list->rl-vector value-list)
  (if (target-extfloat? TARGET)
    (for/extflvector ((x (in-list value-list)))
                     (begin
                      (atom-number x)))
    (for/flvector ((x (in-list value-list)))
                  (begin
                   (atom-number x)))))

(define (normilize-rl rl)
  (if (and (equal? 'racket (target-lang TARGET)) 
           (target-extfloat? TARGET))
    #'(format "~a" rl)
    #'rl))

(define (any-number? n)
  (if (target-extfloat? TARGET)
    (or (extflonum? n) (number? n))
    (number? n)))

