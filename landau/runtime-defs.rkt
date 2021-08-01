#lang racket
(require 
 racket/base
 racket/string
 racket/fixnum
 "target-config.rkt")

(provide (all-defined-out))

;; NOTE: dirty hacks to break hygiene 
(define func_slice_idx 'bad_runtime_func_slice_idx)
(define df_slice_idx 'bad_runtime_df_slice_idx)
(define dx_slice_idx 'bad_runtime_dx_slice_idx)
(define slice_idx 'bad-runtime)
(define al_index_name_symbol "al_index_name_symbol")
(define mappings_full_idx_symbol 'mappings_full_idx_not_set)
(define mapped_idx 'mapped_idx_not_set)

(define debug-runtime #f)
(define (log-debug msg)
  (unless (not debug-runtime)
    (displayln msg)))

(define (not-void? it) (not (equal? (void) it)))

(define BREAK #f)

(define (loop i end body)
  (when (fx< i end)
    (unless (equal? BREAK (body i))
      (loop (add1 i) end body))))

(define (loop-forever i body)
  (let ((ret-code (body i)))
    (unless (equal? BREAK ret-code)
      (loop-forever (add1 i) body))))
