#lang racket
(require (for-syntax racket/base
                     racket/syntax
                     racket/match
                     racket/flonum
                     racket/fixnum
                     racket/extflonum
                     "constant-propagation.rkt"
                     "target-config.rkt"
                     "environment.rkt"))

(provide (for-syntax (all-defined-out)))


(define-for-syntax (instantiate type-stx)
    (let* ((type (syntax->datum type-stx))
           (res (match type
      ((or (list 'real (list size))
           (list 'dual-r (list size))
           (list 'dual-l (list size)))
       (with-syntax ((size-stx (datum->syntax type-stx size)))
          #`(#,make-rl-vector #,#'size-stx)))
      ((list 'int (list size))
       (with-syntax ((size-stx (datum->syntax type-stx size)))
          #`(#,make-rl-vector #,#'size)))
      ((list 'real '())
       #`#,rl0.0)
      ((list 'int '())
       #'0))))
     res))

(define-for-syntax (expand-type type-stx)
  (let* 
    ((type (syntax->datum type-stx))
     (res 
       (match type
         ((list base-type (list size))
          (with-syntax* 
            ((size-stx (local-expand (datum->syntax type-stx size) 'expression '())))
            #`(list (quote #,base-type) (list #,#'size-stx))))
         ((list base-type '())
          #`(list (quote #,base-type) (list)))
         )))
    res))


(define-for-syntax (expand-type-to-datum type-stx)
  (define ns (make-base-namespace))
  (eval (syntax->datum (expand-type type-stx)) ns))


(define-for-syntax (type-range type-stx)
  (let* ((type (syntax->datum type-stx)))
     (match type
      ((list 'list base-type (list 'list size))
       (datum->syntax type-stx size))
      ((list 'list base-type (list 'list))
       (datum->syntax type-stx #'f))
      )))

(define-for-syntax (type-base type-stx)
  (let* ((type (syntax->datum type-stx)))
     (match type
      ((list 'list base-type _)
       base-type)
      )))
