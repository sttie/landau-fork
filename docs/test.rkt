#lang racket

(define-syntax (something stx)
    "(1 2 3)")

(something 1 2 3)