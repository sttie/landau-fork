#lang racket/base

(require racket/match
         racket/fixnum
         racket/extflonum
         racket/vector
         racket/contract
         racket/string
         "environment.rkt"
         "target-config.rkt")

(provide (all-defined-out))

(define indent-step 2)
(define offset (make-parameter 0))
(define (offset-string level) (make-string (fx* level indent-step) #\ ))

(define (c+ . x) (format "(~a)" (string-join (map (lambda (y) (format "~a" y)) x) " + ")))
(define (c- x y) (format "(~a - ~a)" x y))
(define (c* . x) (format "(~a)" (string-join (map (lambda (y) (format "~a" y)) x) " * ")))
(define (c/ x y) (format "(~a / ~a)" x y))
(define (c< x y) (format "(~a < ~a)" x y))
(define (c> x y) (format "(~a > ~a)" x y))
(define (c== x y) (format "(~a == ~a)" x y))
(define (c>= x y) (format "(~a >= ~a)" x y))
(define (c<= x y) (format "(~a <= ~a)" x y))
(define (c-pow x y) (format "pow(~a, ~a)" x y))
(define (c-sqrt x) (format "sqrt(~a)" x))
(define (c-powl x y) (format "powl(~a, ~a)" x y))
(define (c-sqrtl x) (format "sqrtl(~a)" x))
(define (c-neg x) (format "-~a" x))
(define (c-not x) (format "!~a" x))
(define (c-or x y) (format "(~a || ~a)" x y))
(define (c-and x y) (format "(~a && ~a)" x y))
(define (c-cos x) (format "cos(~a)" x))
(define (c-cosl x) (format "cosl(~a)" x))
(define (c-sin x) (format "sin(~a)" x))
(define (c-sinl x) (format "sinl(~a)" x))
(define (c-tan x) (format "tan(~a)" x))
(define (c-tanl x) (format "tanl(~a)" x))

(define (c-for index start end body (pragma #false))
  (let* ((indentation (offset-string (offset)))
         (pragma-str (match pragma
                       ("vectorize"  "#pragma vector always\n")
                       (else "")))
         (for-header (if pragma
                       (string-append indentation pragma-str indentation)
                       indentation)))
    (string-append for-header
                   (format "for (int ~a = ~a; ~a < ~a; ~a++) {\n~a~a}\n" index start index end index
                           (parameterize
                             ((offset (+ (offset) 1)))
                             (call-with-parameterization (current-parameterization) body)) indentation))))

(define (c-forever index start body)
  (let ((indentation (offset-string (offset))))
    (string-append indentation (format "for (int ~a = ~a;; ~a++) {\n~a~a}\n" index start index
                                       (parameterize ((offset (+ (offset) 1))) (call-with-parameterization (current-parameterization) body)) indentation))))

(define (c-break)
  (let ((indentation (offset-string (offset))))
    (format "~abreak;\n" indentation)))

(define (c-with-brackets body)
  (let ((indentation (offset-string (offset))))
    (format "~a{\n~a~a}\n"
            indentation
            (parameterize ((offset (+ (offset) 1))) (call-with-parameterization (current-parameterization) body))
            indentation)))

(define (c-if-expr pred true-body false-body)
  (format "(~a ? ~a : ~a)" pred true-body false-body))

;; NOTE: for experiment with branchless code.
#| (define (c-if-expr pred true-body false-body) |#
#|   (c-or (c-and pred true-body) (c-and (c-not pred) false-body))) |#

(define (c-if pred true-body (false-body #f))
  (let ((indentation (offset-string (offset))))
    (if false-body
      (format "~aif (~a) {\n~a~a} else {\n~a~a}\n"
              indentation
              pred
              (parameterize ((offset (+ (offset) 1))) (call-with-parameterization (current-parameterization) true-body))
              indentation
              (parameterize ((offset (+ (offset) 1))) (call-with-parameterization (current-parameterization) false-body))
              indentation)
      (format "~aif (~a) {\n~a~a}\n"
              indentation
              pred
              (parameterize ((offset (+ (offset) 1))) (call-with-parameterization (current-parameterization) true-body))
              indentation))))

(define (c-indent value)
  (let ((indentation (offset-string (offset))))
    (format "~a~a" indentation value)))

(define (c-return value)
  (let ((indentation (offset-string (offset))))
    (format "~areturn ~a;\n" indentation value)))

(define (c-define type name value (modifier ""))
  (let ((indentation (offset-string (offset))))
    (if (equal? modifier "")
      (format "~a~a ~a = ~a;\n" indentation type name value)
      (format "~a~a ~a ~a = ~a;\n" indentation modifier type name value))))

(define (c-define-array type name size value (modifier ""))
  (let ((indentation (offset-string (offset))))
    (if (equal? modifier "")
      (format "~a~a ~a[~a] = ~a;\n" indentation type name size value)
      (format "~a~a ~a ~a[~a] = ~a;\n" indentation modifier type name size value))))

(define (c-declare-array type name size (modifier ""))
  (let ((indentation (offset-string (offset))))
    (if (equal? modifier "")
      (format "~a~a ~a[~a];\n" indentation type name size)
      (format "~a~a ~a ~a[~a];\n" indentation modifier type name size))))


(define (c-func-decl return-type name return-value args body)
  (format "~a ~a(~a~a) {\n~a}\n\n"
          return-type
          name
          return-value
          args 
          (parameterize ((offset (+ (offset) 1)))
                        (call-with-parameterization (current-parameterization) body))))

(define (c-func-call func-name func-ret-symb args-list)
  (let ((c-args-list (string-join (cons func-ret-symb args-list) ", ")))
    (format "~a(~a)" func-name c-args-list)))

(define (c-pure-func-call func-name args-list)
  (let ((c-args-list (string-join (map to-string args-list) ", ")))
    (format "~a(~a)" func-name c-args-list)))

(define (c-new-line line)
  (format "~a\n" line))

(define (c-line-end line)
  (format "~a;\n" line))

(define (c-array-ref name idx)
  (format "~a[~a]" name idx))

(define (to-string x)
  (if (extflonum? x)
    ;; NOTE: print extflonum as double. Need for C backend 
    (string-replace (format "~a" x) "t" "e") 
    (format "~a" x)))

(define (c-make-array values-list)
  (format "{ ~a }" (string-join (map to-string values-list) ", ")))

(define (vector->carray vec)
  (format "{ ~a }"
          (string-join
           (vector->list (vector-map! number->string vec)) ", ")))

(define (c-set-array symb idx value)
  (let ((indentation (offset-string (offset))))
    (format "~a~a[~a] = ~a;\n" indentation symb idx value)))

(define (c-set symb value)
  (let ((indentation (offset-string (offset))))
    (format "~a~a = ~a;\n" indentation symb value)))

(define (c-dereference pointer-symb)
  (format "*~a" pointer-symb))

(define c-real-type 
  (if (target-extfloat? TARGET)
    "long double"
    "double"))

(define c-zero-filled-array "{ 0.0 }")

(define (c-exact->inexact value)
  (if (target-extfloat? TARGET)
    (format "((long double) ~a)" value)
    (format "((double) ~a)" value)))

;; FIXME use runtime
(define (to-c-decl landau-parsed-type name-str [value #'#f])
  (let ((value (if (syntax? value) (syntax->datum value) value)))
    (with-syntax ((decl-str (match landau-parsed-type
      [(list 'real '()) (if value
                          #`(format "double ~a = ~a;\n" #,name-str #,value)
                          #`(format "double ~a;\n" #,name-str))]
      [(list 'int '()) (if value
                         #`(format "int ~a = ~a;\n" #,name-str #,value)
                         #`(format "int ~a;\n" #,name-str))]
      [(list 'real (list size)) #`(format "static double ~a[~a];\n" #,name-str #,size)]
      [(list 'dual-l (list size)) #`(format "static double ~a[~a];\n" #,name-str #,size)]
      [(list 'dual-r (list size)) #`(format "static double ~a[~a];\n" #,name-str #,size)]
      [(list 'dual-l '()) (if value
                            #`(format "double ~a = ~a;\n" #,name-str #,value)
                            #`(format "double ~a;\n" #,name-str))]
    
      [else (error (format "bug: unsupported type: ~a" landau-parsed-type))])))
#'(string-append (offset-string (offset)) decl-str))))

(define (to-c-der-decl df-type dx-size name-str)  
  (with-syntax ((decl-str (match df-type
         [(list 'dual-l (list size)) #`(format "static double ~a[~a];\n" #,name-str (fx* #,size #,dx-size))]
         [(list 'real (list size)) #`(format "static double ~a[~a];\n" #,name-str (fx* #,size #,dx-size))]
         [(list 'real '()) #`(format "static double ~a[~a];\n" #,name-str #,dx-size)]
    [else (error (format "bug: to-c-der-decl:varaible of type: ~a can not have derivatives" df-type))])))
#'(string-append (offset-string (offset)) decl-str)))

(define (prepand-headers src-str)
  (format "#include <math.h>\n\n~a" src-str))

(define/contract
  (c-declare-var name-str type (modifier-pragma 'on-stack) (value #f) (const? #f))
  (->* (string? landau-type/c) (any/c any/c boolean?)
       string?)
  (let* ((target_ TARGET)
         (c-real-type (if (target-extfloat? target_) "long double" "double"))
         (indentation (offset-string (offset)))
         ;; NOTE: all real arrays are on stack
         (modifier (match modifier-pragma
                     ('static "static ") 
                     ('on-stack "")))
         (const-modifier (if const?
                           "const "
                           ""))
         (decl (match type
                 ((or (list 'real (list size))
                      (list 'dual-r (list size))
                      (list 'dual-l (list size)))
                  (let ((arr-value (if value
                                     value
                                     "{ 0.0 }")))
                    (format "~a~a~a ~a[~a] = ~a;\n" modifier const-modifier c-real-type name-str size arr-value)))

                 ((list 'int (list size))
                  (let ((arr-value (if value
                                     (format " = ~a" value)
                                     "")))
                    (format "~a~a~a ~a[~a]~a;\n" modifier const-modifier "int" name-str size arr-value)))

                 ((list 'real '())
                  (let ((var-value (if value
                                     (format " = ~a" value)
                                     "")))
                    (format "~a~a~a ~a~a;\n" modifier const-modifier c-real-type name-str var-value)))

                 ((list 'int '())
                  (let ((var-value (if value
                                     (format " = ~a" value)
                                     "")))
                    (format "~a~a~a ~a~a;\n" modifier const-modifier "int" name-str var-value))))))
        (string-append indentation decl)))


(define (c-ignore-void value)
  (if (equal? value (void))
    "\n"
    value))

(define (to-c-func-param landau-parsed-type name-str target (is-return-variable? #f))
  (let ((target-real (if (target-extfloat? target) "long double" "double")))
    (match landau-parsed-type
    [(list 'real '()) (if is-return-variable? ;; Function returns using mutation
                        (format "~a* ~a" target-real name-str)
                        (format "~a ~a" target-real name-str))]
    [(list 'int '()) (if is-return-variable?
                       (format "int* ~a" name-str)
                       (format "int ~a" name-str))]
    [(list 'real (list size)) (format "~a *restrict ~a" target-real name-str)]
    [(list 'int (list size)) (format "~a *restrict ~a" "int" name-str)]
    [else (error (format "bug: unsupported type: ~a" landau-parsed-type))])))


(define/contract 
  (c-func-arg-decl argnames args-list)
  (-> syntax? (listof any/c) 
      string?)
  (if (null? (syntax-e argnames))
    ""
    (string-append ", " (string-join args-list ", "))))


(define (c-print-assignation arr-name idx value)
  (let ((indentation (offset-string (offset)))
        (cast-expr (if (target-extfloat? TARGET)
                     "(long double)"
                     "(double)")))
    (format "~aprintf(\"~a[%d] = %.15Le\\n\", ~a, ~a~a);\n" indentation (to-string arr-name) idx cast-expr (to-string value))))
