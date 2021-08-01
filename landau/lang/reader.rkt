#lang s-exp syntax/module-reader
#:language (collection-file-path "semantics.rkt" "landau")
#:read my-read
#:read-syntax my-read-syntax
#:whole-body-readers? #t

(require "../grammar.rkt")
(require "../lexer.rkt")

(define (my-read in)
  (syntax->datum (my-read-syntax #f in)))

(define (my-read-syntax src ip)
  (list (parse src (tokenize ip))))

(define (elem v lst)
  (ormap (lambda (x) (equal? v x)) lst))
