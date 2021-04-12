#lang racket/base
(require 
 ;(for-syntax (for-syntax racket/base racket/syntax syntax/parse))
 (for-syntax racket/base
             racket/contract
             racket/extflonum
             racket/fixnum
             racket/flonum
             racket/function
             racket/list
             racket/match
             racket/set
             racket/string
             racket/stxparam
             racket/syntax
             racket/pretty
             racket/vector
             syntax/parse
             "constant-propagation.rkt" 
             "environment.rkt"
             "target-config.rkt"
             "metalang.rkt")
         "constant-propagation.rkt" 
         "environment.rkt"
         "combinators.rkt"
         "metalang.rkt"
         racket/stxparam
         racket/list
         racket/flonum
         racket/extflonum
         racket/fixnum)

(provide (all-defined-out))
(provide (for-syntax (all-from-out "constant-propagation.rkt")))
(provide (for-syntax (except-out (all-defined-out) binary-op-cast)))

(define-syntax-rule (define/contract-for-syntax definition type body ...)
  (begin-for-syntax (define/contract definition type (begin body ...))))

(define-for-syntax TIME-TABLE (make-hash (list (cons 'head-time 0.0)
                                               (cons 'tail-time 0.0)
                                               (cons 'local-expand-2 0.0)
                                               (cons 'to-landau-type 0.0)
                                               (cons 'search-name 0.0)
                                               (cons 'local-expand 0.0))))
(define-for-syntax (fmap f maybe-false)
  (if maybe-false
      (f maybe-false)
      #f))

(define-for-syntax (zip proc l1 l2)
  (for/list ((lr (in-list l1)) (ll (in-list l2)))
    (proc lr ll)))


(define-for-syntax (with-syntax-property property-name property-value stx)
  (syntax-property stx property-name property-value #t))

(define-for-syntax (read-syntax-property property-name stx)
  (syntax-property stx property-name))

(define-for-syntax (is-type stx type)
  (syntax-property stx 'landau-type type #t))

(define-for-syntax (is-type_ type stx)
  (syntax-property stx 'landau-type type #t))

(define-for-syntax (is-int stx)
  (is-type stx 'int))

(define-for-syntax (is-real stx)
  (is-type stx 'real))

(define-for-syntax (is-dual stx)
  (is-type stx 'dual))

(define-for-syntax (has-dual stx)
  (syntax-property stx 'has-dual #t))

(define-for-syntax (is-array-type? type)
  (not (equal? (car (cdr type)) '())))

(define-for-syntax (get-array-type type)
  (car type))

(define-for-syntax (fxvector->vector fxvec)
  (for/vector
      ((i (in-range (fxvector-length fxvec))))
    (fxvector-ref fxvec i)))

(define-for-syntax (get-value-stx dual-b-expanded-stx)
  (caddr (syntax->list dual-b-expanded-stx)))

(define-for-syntax (get-derivative-stx dual-b-expanded-stx)
  (cadddr (syntax->list dual-b-expanded-stx)))

(define-for-syntax (is-slice? type)
  (if (list? type)
      ;; NOTE: Cant check that second value is nonegative-integer because it is not known in compile time
      (equal? (length type) 2)
      #f))

(define-for-syntax (get-slice-range type)
  (match type
    ((list base (list rng)) rng)
    (base #f)))

(define-for-syntax (get-slice-type type)
  (match type
    ((list base _) base)
    (base #f)))

(define-for-syntax (is-slice-of-type basic-type type)
  (and (is-slice? type) (equal? (get-slice-type type) basic-type)))

(define-for-syntax (landau-type base-type (slice-range #f))
  (if slice-range
      (list base-type (list slice-range))
      base-type))

(define-for-syntax (check-result stx err-str func-res)
  (if func-res
      func-res
      (raise-syntax-error #f err-str stx)))

(define-for-syntax (binary-op-cast op1 op2 stx) (binary-op-cast-helper op1 op2 stx))


(define-for-syntax (binary-op-cast-helper op1 op2 stx)
  (let ((type1 (syntax-property op1 'landau-type))
        (type2 (syntax-property op2 'landau-type))
        (op1-atom (atom-number op1))
        (op2-atom (atom-number op2)))
    (cond
      ((and (equal? type1 'int) (equal? type2 'int))
       (values op1 op2 'int))
      ((and (equal? type1 'real) (equal? type2 'real))
       (values op1 op2 'real))
      ((and (equal? type1 'int) (equal? type2 'real))
       (values (if op1-atom
                 (datum->syntax stx `,(->rl op1-atom))
                 (datum->syntax stx `(_exact->inexact ,op1)))
               op2 'real))
      ((and (equal? type1 'real) (equal? type2 'int))
       (values op1 (if op2-atom
                     (datum->syntax stx `,(->rl op2-atom))
                     (datum->syntax stx `(_exact->inexact ,op2)))
               'real))

      ((and (equal? type1 'dual-b) (equal? type2 'dual-b))
       (values op1 op2 'dual-b))
      ((and (equal? type1 'dual-b) (equal? type2 'real))
       (values op1 op2 'dual-b))
      ((and (equal? type1 'real) (equal? type2 'dual-b))
       (values op1 op2 'dual-b))
      ((and (equal? type1 'dual-b) (equal? type2 'int))
       (values op1 (if op2-atom
                     (datum->syntax stx `,(->rl op2-atom))
                     (datum->syntax stx `(_exact->inexact ,op2)))
               'dual-b))
      ((and (equal? type1 'int) (equal? type2 'dual-b))
       (values (if op1-atom
                 (datum->syntax stx `,(->rl op1-atom))
                 (datum->syntax stx `(_exact->inexact ,op1)))
               op2 'dual-b))

      ((and (is-slice? type1) (is-slice? type2))
       (let ((range1 (get-slice-range type1))
             (range2 (get-slice-range type2))
             (slice-type1 (get-slice-type type1))
             (slice-type2 (get-slice-type type2)))
         ;; NOTE: ranges should be checked in backran, ranges can have loop variable in tnem, thus they are unknown in compile time
         #;(unless (equal? range1 range2)
             (raise-syntax-error #f (format "cannot cast ranges ~v and ~v" range1 range2) stx))
         (match (list slice-type1 slice-type2)
           ((list 'dual-b 'dual-b) (values op1 op2 type1))
           ((list 'dual-b 'real) (values op1 op2 type1))
           ((list 'real 'dual-b) (values op1 op2 type2))
           ((list 'real 'real) (values op1 op2 type1))
           (else
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type1 slice-type2) stx)))))
      ((and (is-slice? type1) (not (is-slice? type2)))
       (let ((range1 (get-slice-range type1))
             (slice-type1 (get-slice-type type1)))
         (match (list slice-type1 type2)
           ((list 'dual-b 'dual-b) (values op1 op2 type1))
           ((list 'dual-b 'real) (values op1 op2 type1))
           ((list 'dual-b 'int) (values op1 (if op2-atom
                                              (datum->syntax stx `,(->rl op2-atom))
                                              (datum->syntax stx `(_exact->inexact ,op2)))
                                        type1))
           ((list 'real 'dual-b) (values op1 op2 (landau-type 'dual-b (get-slice-range type1))))
           ((list 'real 'real) (values op1 op2 type1))
           ((list 'real 'int) (values op1 (if op2-atom
                                            (datum->syntax stx `,(->rl op2-atom))
                                            (datum->syntax stx `(_exact->inexact ,op2)))
                                      type1))
           (else
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type1 type2) stx)))))

      ((and (is-slice? type2) (not (is-slice? type1)))
       (let ((range2 (get-slice-range type2))
             (slice-type2 (get-slice-type type2)))
         (match (list slice-type2 type1)
           ((list 'dual-b 'dual-b) (values op1 op2 type2))
           ((list 'dual-b 'real) (values op1 op2 type2))
           ((list 'dual-b 'int) 
            (values (if op2-atom
                      (datum->syntax stx `,(->rl op1-atom))
                      (datum->syntax stx `(_exact->inexact ,op1)))
                    op2
                    type2))
           ((list 'real 'dual-b) (values op1 op2 (landau-type 'dual-b (get-slice-range type2))))
           ((list 'real 'real) (values op1 op2 type2))
           ((list 'real 'int) (values (if op1-atom
                                        (datum->syntax stx `,(->rl op1-atom))
                                        (datum->syntax stx `(_exact->inexact ,op1)))
                                      op2
                                      type2))
           (else
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type2 type1) stx)))))
      (else
       (begin
        (raise-syntax-error #f (format "cannot cast ~v and ~v" type1 type2) stx))))))


(define-for-syntax (xor e1 e2)
  (let ((b1 (not (not e1)))
        (b2 (not (not e2))))
    (cond
      [(equal? b1 b2) #f]
      [else #t])))
      

(define/contract-for-syntax
  (make-need-derivatives-key arr-name-str)
  (-> string? string?)
  arr-name-str)

(define/contract-for-syntax
  (get-der-variable-dx-range-helper df-name dx-name need-derivatives-table stx)
  (-> var-symbol/c (or/c string? false/c) need-derivatives-table/c (syntax/c any/c)
      (or/c integer? false/c))
  (cond
    ((not dx-name) #f)
    (else
      (fmap mapping-sizes-.mapping-period 
            (ndt-get-mapping-sizes need-derivatives-table df-name dx-name)))))





;; NOTE: iterate through args, if arg need derivatives then
;; - add -der variable
;; - change type 'real -> 'dual-l
;; - allocate vector
;; - dx is always dual-l
(define/contract-for-syntax
  (make-der-decl-list-helper!
    args
    dx-names-hash-set 
    need-derivatives-table
    current-variables
    current-arguments
    lang
    stx)
  (-> current-arguments/c
      dx-names-hash-set/c
      need-derivatives-table/c
      current-variables/c
      current-arguments/c
      lang/c
      (syntax/c any/c)

      (listof any/c))

  (let* ((dx-names dx-names-hash-set)
         (fake-src-pos 0))
    (for/list ((arg-kv (in-list (hash->list args))))
      (let* ((arg-name (car arg-kv))
             (arg-vs (var-symbol (symbol->string arg-name) fake-src-pos))
             (arg-struct (cdr arg-kv))
             (arg-basic-type (car (argument-type arg-struct))))
        (when (equal? arg-basic-type 'real)
          (let*
            ((df-type (get-arg-type-helper (datum->syntax stx arg-name) current-arguments)))
            (when (hash-has-key? dx-names (symbol->string arg-name))
              (change-arg-type!
                current-arguments
                arg-name
                (make-dual-type df-type 'dual-l)))
            (unless (var-symbol/c arg-vs)
              (error "436"))
            (when (ndt-member? need-derivatives-table arg-vs)
              (with-syntax
                ((der-vec-decl-list-forall-dx
                   (let ((dx-keys (ndt-get-dx-names need-derivatives-table arg-vs)))
                     (for/list ((dx-name-str (in-list dx-keys)))
                       (let* ((df-type (get-arg-type-helper (datum->syntax stx arg-name) current-arguments))

                              (arg-name-str (symbol->string arg-name))
                              (dx-mapped-size (check-result
                                                stx 
                                                (format 
                                                  "bug: get-der-variable-dx-range-helper returned #f for a key: ~a ~a" 
                                                  arg-vs
                                                  dx-name-str) 
                                                (get-der-variable-dx-range-helper 
                                                  arg-vs
                                                  dx-name-str 
                                                  need-derivatives-table 
                                                  stx)))
                              (df-size (car (if (empty? (cadr df-type)) (list 1) (cadr df-type))))
                              (dx-mapped-type (make-landau-type 'real dx-mapped-size))
                              (dual-r-var (add-variable!
                                            current-variables
                                            (datum->syntax #f (make-var-der-name arg-vs dx-name-str))
                                            ;; NOTE: it is a period number (dx range), not the vector capasity
                                            (make-dual-type dx-mapped-type 'dual-r))))
                         ;(log-debug (format "debug: add derivative. key: ~a" (make-var-der-name arg-name dx-name-str)))
                         ;; NOTE: After dual-r-var been added to variables
                         ;;       arg type is changed from 'real to 'dual-l
                         (change-arg-type!
                           current-arguments
                           arg-name
                           (make-dual-type df-type 'dual-l))
                         (match lang 
                           ('racket
                            (with-syntax ((der-vec (datum->syntax stx dual-r-var)))
                              (quasisyntax/loc stx
                                               (define der-vec #,(instantiate-dual-var df-type dx-mapped-size stx)))))
                           ('ansi-c
                            (with-syntax* ((der-vec (symbol->string dual-r-var))
                                           (df-size (local-expand (datum->syntax stx df-size) 'expression '()))
                                           (arr-size (fx* (atom-number #'df-size) dx-mapped-size)))
                                          (quasisyntax/loc stx  
                                                           (c-define-array c-real-type
                                                                           der-vec arr-size c-zero-filled-array))))))))))
                (if (equal? lang 'racket)
                  (cons 'begin (syntax-e #'der-vec-decl-list-forall-dx))
                  (with-syntax ((der-vec-decl-list-forall-dx_ (cons 'list (syntax-e #'der-vec-decl-list-forall-dx))))
                    #'(flatten der-vec-decl-list-forall-dx_)))))))))))


(define-for-syntax (get-arg-type-helper d-name current-arguments)
  (let
      ((d-name-arg (search-argument
                    (syntax->datum d-name) current-arguments)))
    (cond
      (d-name-arg
       (let ((base-type (car (argument-type d-name-arg))))
         (if (or (equal? 'real base-type) (equal? 'dual-l base-type))
             (argument-type d-name-arg)
             (raise-syntax-error #f "derivative of non-real in prohibited" d-name))))
      (else
       (raise-syntax-error #f "name not found; only function parameters can be used in annotation" d-name)))))

;; NOTE: it is used close to the declaration where types are just parsed,
;; except the func-value case, where it's type changed from real to dual-l and then der-vector is allocated
(define-for-syntax (instantiate-dual-var df-type dx-mapped-size stx)
  (let ((msg "bug: instantiate-dual-var: df-type can be only 'real array"))
    (cond
      ((exact-positive-integer? dx-mapped-size)
       (match (list df-type dx-mapped-size)
         ((list (list type (list df-size)) dx-mapped-size)
          (if (or (equal? type 'real) (equal? type 'dual-l))
              (with-syntax ((expanded-df-size (cadr (syntax->datum (local-expand (datum->syntax stx df-size) 'expression '())))))
                (quasisyntax/loc stx 
                  (#,make-rl-vector (fx* #,#'expanded-df-size #,dx-mapped-size))))
              (error msg)))
         ((list (list type (list)) dx-mapped-size)
          (if (or (equal? type 'real) (equal? type 'dual-l))
              (quasisyntax/loc stx
                (begin
                  (#,make-rl-vector #,dx-mapped-size)))
              (error msg)))
         (else "bug: instantiate-dual-var: match type failed")))
         
      (else (error (format "bug: instantiate-dual-var: exact-positive-integer? dx-mapped-size constraint violation. Given: ~a" dx-mapped-size))))))

(define-for-syntax (set-need-only-value-set! grouped-keys-table der-table need-only-value-set)
  (for/list ((var-name-key (in-list (hash-keys grouped-keys-table))))
    (let* ((is-basic-variable (equal? (hash-ref grouped-keys-table var-name-key) (list #f))))
      (cond
        (is-basic-variable
         ;; NOTE: non-array variable
         (let ((df-table (hash-ref der-table (ref-to-key (list 'basic-ref var-name-key)))))  
           (when (hash-empty? df-table)
             (set-add! need-only-value-set var-name-key))))
        (else
         ;; NOTE: array variable
         (let* ((df-tables-are-empty
                 (for/and ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
                   (let ((df-table (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))
                     (hash-empty? df-table)))))
           (when df-tables-are-empty
             (set-add! need-only-value-set var-name-key))))))))

;; NOTE: dual-type: dual-l | dual-r
;; TODO: change parsed-type to df-type and maybe-dx-type and multiply sizes
(define-for-syntax (make-dual-type parsed-type dual-type)
  (match parsed-type
    [(list basic (list size)) (list dual-type (list size))]
    [(list 'dual-l '()) (list dual-type '())]
    [(list 'real '()) (list dual-type '())]
    [else (error (format "unsupported type: ~v" parsed-type))]))
    

(define/contract-for-syntax
  (check-duplicate-variable-name-helper name stx-name current-variables current-arguments function-name)
  (-> symbol? (syntax/c any/c) current-variables/c current-arguments/c symbol? void?)
  #| (when (search-variable name current-variables) |#
       #|   (pretty-print current-variables) |#
       #|   (raise-syntax-error #f "duplicate variable declaration" stx-name)) |#
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "variable name shadows constant" stx-name))
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "variable name shadows parameter" stx-name))
  (when (hash-has-key? current-arguments name)
    (raise-syntax-error #f "variable name shadows argument" stx-name))
  (when (equal? function-name name)
    (raise-syntax-error #f "variable name shadows function" stx-name)))

(define/contract-for-syntax
  (check-duplicate-argument-name args funcname name stx-name)
  (-> current-arguments/c symbol? symbol? (syntax/c any/c) void?)
  (when (equal? funcname name)
    (raise-syntax-error #f "argument name shadows function" stx-name))
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "variable name shadows parameter" stx-name))
  (when (hash-has-key? args name)
    (raise-syntax-error #f "duplicate argument" stx-name))
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "argument name shadows constant" stx-name)))

;; NOTE: extract syntax after syntax-parameterize
(define-for-syntax (extract stx)
  (syntax-parse stx
    (((~literal let-values) _ 
                            ((~literal #%expression) 
                             ((~literal let-values) _ expr))) #'expr)
    (((~literal let-values) _ 
                            ((~literal #%expression) 
                             ((~literal let-values) _ _ _ expr))) #'expr)))

(define-for-syntax (fxvec->vec fxvec)
  (for/vector ((i (fxvector-length fxvec)))
    (fxvector-ref fxvec i)))

(define-for-syntax (vec->fxvec vec)
  (for/fxvector ((i (vector-length vec)))
    (vector-ref vec i)))

(define (fxvec->vec fxvec)
  (for/vector ((i (fxvector-length fxvec)))
    (fxvector-ref fxvec i)))

(define (vec->fxvec vec)
  (for/fxvector ((i (vector-length vec)))
    (vector-ref vec i)))

(define-for-syntax debug #f)

(define/contract-for-syntax
  (get-symbol-helper stx name dx-name-str make-name current-variables err-msg [throw-error? #t])
  (->* ((syntax/c any/c) 
        var-symbol/c 
        string? 
        (-> var-symbol/c string? symbol?)
        current-variables/c
        string?)
       (boolean?)
       (or/c (syntax/c symbol?) #f))
  (let ([var
          (search-variable
            (make-name name dx-name-str)
            current-variables)])
    (cond
      [var (datum->syntax stx (variable-symbol var))]
      [else (if throw-error?
              (raise-syntax-error #f (format "bug: ~a array not found. key: ~a"
                                             err-msg
                                             (make-name name dx-name-str)) stx)
              #f)])))


(define (is-basic-ref? ref) (equal? (car ref) 'basic-ref))

(define (is-array-ref? ref) (equal? (car ref) 'array-ref))

(define/contract-for-syntax
  (_key<? vs-1 vs-2)
  (-> var-symbol/c var-symbol/c 
      boolean?)
  (string<? (var-symbol->string vs-1)
            (var-symbol->string vs-2)))

(define/contract-for-syntax
  (group-by-names keys)
  (-> (listof df-name/c) 
      grouped-keys-table/c)
  (let* ((sorted (sort keys _key<? #:key car))
         (ht (make-hash)))
    (for ((p (in-list sorted)))
      (let ((name (car p))
            (idx (if (equal? (length p) 2)
                   (cadr p)
                   #f)))
        (if (hash-has-key? ht name)
          (let* ((l (hash-ref ht name)))
            (hash-set! ht name (append (list idx) l)))
          (hash-set! ht name (list idx)))))
    ht))

(begin-for-syntax
  (define-syntax-class plus-or-minus
    (pattern "+")
    (pattern "-"))
  (define-syntax-class mul-or-div
    (pattern "*")
    (pattern "/"))
  (define-syntax-class colon
    (pattern ":")))

(begin-for-syntax
  (define-syntax-class par-spec
    #:attributes (par-value)
    (pattern
     ({~literal par} par-value))

    (pattern
     ({~literal other-par} "," par-value))))

(begin-for-syntax
 (define-splicing-syntax-class expr-body-cls
   #:attributes (body)
   (pattern
    (~seq "{" body "}"))
   (pattern
    (~seq body))))


(define-for-syntax (get-slice-start-and-range stx slice-colon index-start index-end array-range)
  (if (and (syntax->datum array-range)
           slice-colon)
    (with-syntax* 
      ((index-start-expanded (if (syntax->datum index-start)
                               (local-expand index-start 'expression '())
                               0))
       (index-end-expanded (if (syntax->datum index-end)
                             (local-expand index-end 'expression '())
                             (car (get-type-range (to-landau-type stx (list 'int array-range))))))
       (index-start-expanded_ (if (atom-number #'index-start-expanded)
                                (atom-number #'index-start-expanded)
                                #'index-start-expanded))
       (index-end-expanded_ (if (atom-number #'index-end-expanded)
                              (atom-number #'index-end-expanded)
                              #'index-end-expanded))
       (slice-range (if (and (number? (syntax->datum #'index-end-expanded_)) 
                             (number? (syntax->datum #'index-start-expanded_)))
                      (fx- (syntax->datum #'index-end-expanded_) (syntax->datum #'index-start-expanded_))
                      (datum->syntax stx `(_int- ,#'index-end-expanded_ ,#'index-start-expanded_)))))
      (values #'index-start-expanded_ #'slice-range))
    (values #f #f)))


(define/contract-for-syntax
  (check-func stx func-vs func-args funcs-info current-variables current-arguments)
  (-> (syntax/c any/c) var-symbol/c (listof any/c) funcs-info/c current-variables/c current-arguments/c
      void?)
  (unless (hash-has-key? funcs-info func-vs)
    (raise-syntax-error #f (format "function ~a is not defined" (var-symbol-.name func-vs)) stx))
  (let* ((func-info (hash-ref funcs-info func-vs))
         (arity (func-info-arity func-info))
         (args-type (func-info-args-type func-info))
         (output-type (func-info-output-base-type func-info))
         (output-range (func-info-output-range func-info)))
    (unless (equal? arity (length func-args))
      (raise-syntax-error
        #f
        (format "function ~a arity mismatch: expects ~a, but given ~a" (var-symbol-.name func-vs) arity (length func-args))
        stx))
    (for ((expected-arg-type (in-list args-type))
          (provided-arg (in-list func-args))
          (arg-number (in-naturals 1)))
      ;; provided arg could be:
      ;; expr resoled to a slice expresion (not implemented)
      ;; expr, resolved to a number
      ;; expr, resolved to a variable
      (let* ((expanded-arg (local-expand provided-arg 'expression '()))
             (expanded-arg-type (syntax-property expanded-arg 'landau-type))
             (var-name (syntax-property expanded-arg 'get-value-name))
             (var-is-slice? (is-slice? expanded-arg-type))
             (var-type (if var-name 
                         (let ((var (search-variable var-name current-variables)))
                           (if var 
                             (variable-type var)
                             (let ((arg (search-argument var-name current-arguments)))
                               (if arg
                                 (argument-type arg)
                                 #f))))
                         #f))
             (var-base-type (if var-type (car var-type) expanded-arg-type))
             (var-range (if var-type (cadr var-type) '()))
             (var-range-expanded (if (empty? var-range)
                                   #f 
                                   (atom-number (local-expand (datum->syntax stx (car var-range)) 'expression '()))))

             (expected-arg-base-type (car expected-arg-type))
             (expected-arg-range (cdr expected-arg-type)))

        (when (and expected-arg-range (not var-range-expanded))
          (raise-syntax-error
            #f
            (format "argument ~a should be an array" arg-number)
            stx))

        (when (and (not expected-arg-range) var-range-expanded)
          (raise-syntax-error
            #f
            (format "argument ~a should not be an array" arg-number)
            stx))

        (when var-is-slice?
          (raise-syntax-error
            #f
            (format "slices can not be passed to funcitons. Function ~a was given a slice as a ~a argument" (var-symbol-.name func-vs) arg-number)
            stx))


        (unless (equal? (format "~a" expected-arg-base-type) (format "~a" var-base-type))
          (unless (equal? var-base-type 'dual-l) 
            (raise-syntax-error
              #f
              (format "argument ~a type mismatch. Expects ~a, but ~a given."
                      arg-number
                      expected-arg-base-type
                      var-base-type)
              stx)))

        (unless (equal? expected-arg-range var-range-expanded)
          (raise-syntax-error
            #f
            (format "argument ~a range mismatch. Expects array of length ~a, but the provided argumnt has length of ~a"
                    arg-number
                    expected-arg-range
                    var-range-expanded)
            stx))))))

; (begin-for-syntax
;   (define/contract
;     (get-dx-names need-derivatives-table df-name)
;     (-> need-derivatives-table/c symbol? 
;         (listof string?))
   
;     (let ((key (make-need-derivatives-key (symbol->string df-name))))
;       (if (hash-has-key? need-derivatives-table key)
;           (hash-keys (hash-ref need-derivatives-table key))
;           #f))))



(define-for-syntax (not-void? it) (not (equal? (void) it)))

(begin-for-syntax
 (define-splicing-syntax-class getter-cls
  ;  #:attributes (name)
   (pattern
    (~seq "[" (~optional _index-start #:defaults ((_index-start #'#f)))
                                   _slice-colon:colon
                                   (~optional _index-end #:defaults ((_index-end #'#f))) "]")
    #:attr index #'#f
    #:attr slice-colon #'_slice-colon
    #:attr index-start #'_index-start
    #:attr index-end #'_index-end)
   ;; NOTE: pattern matching order does matter
   (pattern 
    (~seq "[" _index "]")
    #:attr index #'_index
    #:attr slice-colon #'#f
    #:attr index-start #'#f
    #:attr index-end #'#f)

   (pattern
    (~seq)
    #:attr index #'#f
    #:attr slice-colon #'#f
    #:attr index-start #'#f
    #:attr index-end #'#f)
   )

 (define-syntax-class get-value-cls
   (pattern
    (_ _name:id pat:getter-cls)
    #:attr name #'_name
    #:attr index #'pat.index
    #:attr slice-colon #'pat.slice-colon
    #:attr index-start #'pat.index-start
    #:attr index-end #'pat.index-end)
   ))

(begin-for-syntax
 (define (timeit! time-table time-key proc)
   (define tik (current-inexact-milliseconds))
   (define r (proc))
   (define tok (current-inexact-milliseconds))
   (hash-update! time-table time-key (lambda (old-time) (fl+ old-time (fl- tok tik))) 0.0)
   (hash-update! time-table (string->symbol (format "~a-counts" time-key))
                 (lambda (x) (add1 x)) 0)
   r)

 (define (timeit/values! time-table time-key proc)
   (define tik (current-inexact-milliseconds))
   (define vals (call-with-values proc list))
   (define tok (current-inexact-milliseconds))
   (hash-update! time-table time-key (lambda (old-time) (fl+ old-time (fl- tok tik))) 0.0)
   (hash-update! time-table (string->symbol (format "~a-counts" time-key))
                 (lambda (x) (add1 x)) 0)
   (apply values vals))


 (define
   (search-name stx ctx name name-stx dx-name-str-in-current-al get-name-mode-on idx)
   #| (-> (syntax/c any/c) |# 
          #|     (or/c func-context/c #f) |# 
          #|     symbol? |# 
          #|     (syntax/c any/c) |# 
          #|     (or/c string? false/c) |#
          #|     boolean? |#
          #|     (syntax/c any/c) |#
          #|     ;; NOTE: (resolved-value type name-is-dx src-pos) |#
          #|     (values (or/c any-number? (syntax/c any/c)) type/c boolean? integer?)) |#

   (let 
     ;; NOTE: src-pos is used to separate variables with the same name
     ((fake-src-pos 0)
      (not-dx-name #f))
     (cond
       ((timeit! TIME-TABLE 'search-constant (thunk (hash-has-key? constants name)))
        (let ((c (hash-ref constants name)))
          (let ((type (constant-type c)))
            (if (and (is-array-type? type) (atom-number idx))
              (let ((idx-atom-number (atom-number idx)))
                (values ;; Where to store const& in compile time or runtime? not all indexes could be resolved
                  (rl-vector-ref (constant-array c) idx-atom-number)
                  (constant-type c)
                  not-dx-name
                  fake-src-pos))
              (values
                (datum->syntax stx (constant-value c))
                (constant-type c)
                not-dx-name
                fake-src-pos)))))
       (else
         (begin
           (when (equal? ctx #f)
             (raise-syntax-error #f "name not found" name-stx))
           (let ((func-name (func-context-.function-name ctx))
                 (var (timeit! TIME-TABLE 'search-variable (thunk (search-variable
                                                         name (func-context-.current-variables ctx))))))
             (cond
               (var
                 (values (datum->syntax stx (variable-symbol var))
                         (variable-type var)
                         (equal? (symbol->string name) dx-name-str-in-current-al)
                         (variable-src-pos var)))
               ((equal? name func-name)
                (let ((func-return-value (func-context-.function-return-value ctx)))
                  (values
                    (datum->syntax stx func-return-value)
                    (func-context-.function-return-type ctx)
                    not-dx-name
                    fake-src-pos)))
               (else
                 (let ((arg (timeit! TIME-TABLE 'search-argument (thunk (search-argument
                                                               name (func-context-.current-arguments ctx))))))
                   (cond
                     (arg
                       (values
                         (datum->syntax stx (argument-symbol arg))
                         (argument-type arg)
                         (equal? (symbol->string name) dx-name-str-in-current-al)
                         fake-src-pos))
                     ((and get-name-mode-on (hash-has-key? parameters name))
                      (values
                        name-stx
                        (list 'real (hash-ref parameters name))
                        not-dx-name
                        fake-src-pos))
                     (else
                       (begin
                         (pretty-print (func-context-.current-variables ctx))
                         (raise-syntax-error #f "name not found" name-stx))))))))))))))
