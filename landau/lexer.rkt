#lang racket

(require parser-tools/lex
         (prefix-in : parser-tools/lex-sre))
(require ragg/support)
(require (for-syntax racket/base syntax/parse))

(provide tokenize)

(define-lex-abbrevs
  (lex:letter (:or (:/ #\a #\z) (:/ #\A #\Z) #\_))
  (lex:digit (:/ #\0 #\9))
  (lex:comment (:: #\# (:* (:~ #\newline)) #\newline))
  (lex:whitespace (:or #\newline #\return #\tab #\space #\vtab))
  (lex:identifier (:: lex:letter (:* (:or lex:letter lex:digit))))
  (lex:integer (:+ lex:digit))
  (lex:float (:or
               (:: (:+ lex:digit) #\. (:* lex:digit)
                   (:? (:: (:or "e" "E") (:? (:or #\- #\+)) (:+ lex:digit))))
               (:: #\. (:+ lex:digit)
                   (:? (:: (:or "e" "E") (:? (:or #\- #\+)) (:+ lex:digit))))
               (:: (:+ lex:digit)
                   (:: (:or "e" "E") (:or #\- #\+) (:+ lex:digit)))))
  )

(define-syntax (lexer-src-pos-with-tokens stx)
  (syntax-parse stx
    ((_ tokens other ...)
     #`(lexer-src-pos
        #,@(for/list ((t (syntax->datum #'tokens)))
             (list t (list 'token t 'lexeme)))
        other ...))))

(define (tokenize ip)
  (port-count-lines! ip)
  (define my-lexer
    (lexer-src-pos-with-tokens
     ("+" "-" "*" "/" "^" "(" ")" "," "[" "]"
          "{" "}" "real" "int"
          "const" "=" "+=" "-=" "*=" "/=" 
          "if" "else" "True" 
          ">" "<" ">=" "<=" "=="
          "!=" "and" "or" "not" "for" ":"
          "'" "discard" "parameter" "print")
     [lex:integer
      (token 'INTEGER (string->number lexeme))]
     [lex:float
      (token 'FLOAT (string->number lexeme))]
     [lex:identifier (token 'IDENTIFIER (string->symbol lexeme))]
     [lex:comment
           (token 'COMMENT lexeme #:skip? #t)]
     [whitespace
      (token 'WHITESPACE lexeme #:skip? #t)]
     [(eof)
      (void)]))
  (define (next-token) (my-lexer ip))
  next-token)
