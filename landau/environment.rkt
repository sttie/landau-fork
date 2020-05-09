#lang racket/base

(require 
 racket/syntax
 racket/match
 racket/base
 racket/contract/base
 racket/contract/region
 racket/contract/combinator
 racket/flonum
 racket/extflonum
 racket/function
 "target-config.rkt")

(provide (all-defined-out))

(define constants (make-hash))
(define parameters (make-hash))

(define (syntax->string stx)
  (symbol->string (syntax->datum stx)))

(define (hash->string ht)
  (for/fold ((str ""))
            (((k v) (in-hash ht)))
             (string-append str (format "~a: ~a\n" k v))))

(struct constant
  (value
   type
   array)
  #:prefab)

(struct variable
  (symbol
   type
   src-pos)
  #:prefab)

(struct variables
  (current-level
   next-level)
  #:prefab)

(struct argument
  (symbol
   type)
  #:prefab)

;; NOTE: To resolve ambiguity between local vars with the same name
(struct var-symbol
        (.name
         .src-pos)
        #:prefab)

;; NOTE: types

(define var-symbol/c
  (struct/c var-symbol
            string?
            integer?))
(define/contract
  (var-symbol->string vs)
  (-> var-symbol/c string?)
  (format "~a~a" (var-symbol-.name vs) (var-symbol-.src-pos vs)))

(struct mapping-sizes 
        (.mapping-period 
         .inv-mapping-period) 
        #:prefab)
(define mapping-sizes/c (struct/c mapping-sizes integer? integer?))
(define (mutable-set/c a) any/c)
(define df-name/c (list/c var-symbol/c integer?))
(define dx-name/c string?)
(define df-der-table/c (hash/c dx-name/c (list/c (symbols 'dx-idxs) mutable-set/c)))

(begin
 (define der-table/c
   (hash/c
    df-name/c
    df-der-table/c))
 (define/contract (dtbl-new)
   (-> der-table/c)
   (make-hash))
 
 (define/contract (dtbl-get-df-table dtbl df-name)
   (-> der-table/c df-name/c
       df-der-table/c)
   (hash-ref dtbl df-name (thunk (error (format "bug: der-table has no ~a key" df-name))))))


(begin
 (define real-vars-table/c (hash/c var-symbol/c (or/c integer? false/c)))
 (define/contract (rvt-new)
   (-> real-vars-table/c)
   (make-hash))

 (define/contract (rvt-get-var-size rvt var)
   (-> real-vars-table/c var-symbol/c 
       (or/c integer? false/c))
   (hash-ref rvt var (thunk (error (format "bug: real-vars-table has no key ~a\n~a" var (hash->string rvt))))))
 
 (define/contract (rvt-var-is-array rvt var)
   (-> real-vars-table/c var-symbol/c
       boolean?)
   (integer? (rvt-get-var-size rvt var))))


(define dx-names-set/c (hash/c string? boolean?))
(define grouped-keys-table/c (hash/c var-symbol/c (listof integer?)))


(begin
 (define dx-mapping-sizes-ht/c (hash/c string? mapping-sizes/c))
 (define/contract (dms-member? dms dx-name)
   (-> dx-mapping-sizes-ht/c string?
       boolean?)
   (hash-has-key? dms dx-name))
 
 (define/contract (dms-get-mapping-sizes dms dx-name)
   (-> dx-mapping-sizes-ht/c string?
       mapping-sizes/c)
   (hash-ref dms dx-name (thunk (error (format "bug: dx-mapping-sizes-ht has no ~a key" dx-name))))))


(begin
 (define need-derivatives-table/c (hash/c var-symbol/c dx-mapping-sizes-ht/c))
 (define/contract (ndt-new)
   (-> need-derivatives-table/c)
   (make-hash))

 (define/contract (ndt-member? ndt df-name)
   (-> need-derivatives-table/c var-symbol/c
       boolean?)
   (hash-has-key? ndt df-name))

 (define/contract (ndt-set! ndt df-name value-ht)
   (-> need-derivatives-table/c var-symbol/c dx-mapping-sizes-ht/c
       boolean?)
   (hash-set! ndt df-name value-ht))

 (define/contract (ndt-get-dx-mapping-sizes-ht ndt df-name)
   (-> need-derivatives-table/c var-symbol/c 
       dx-mapping-sizes-ht/c)
   (hash-ref! ndt df-name (thunk (error (format "bug: ndt-get-dx-mapping-sizes-ht has no key ~a" df-name)))))

 (define/contract (ndt-get-dx-names ndt df-name)
   (-> need-derivatives-table/c var-symbol/c
       (listof string?))
   (define df-table (ndt-get-dx-mapping-sizes-ht ndt df-name))
   (hash-keys df-table))

 (define/contract (ndt-get-mapping-sizes ndt df-name dx-name)
   (-> need-derivatives-table/c var-symbol/c string?
       (or/c mapping-sizes/c false/c))
   (cond
     [(ndt-member? ndt df-name)
      (let ((dx-name->dx-sizes (ndt-get-dx-mapping-sizes-ht ndt df-name)))
         (cond
           ((dms-member? dx-name->dx-sizes dx-name)
            (dms-get-mapping-sizes dx-name->dx-sizes dx-name))
           (else #f)))]
     [else #f])))


(define base-type/c (one-of/c 'real 'int 'dual-l 'dual-r))
(define landau-type/c (list/c base-type/c (or/c (list/c integer?) (list/c))))
;; FIXME: landau-type/c with unexpanded `size` part
; (define type/c (list/c base-type/c any/c))
(define type/c any/c)
(define variable/c (struct/c variable 
                             symbol?
                             type/c 
                             (or/c integer? false/c)))
(define current-variables/c
  (struct/c variables
            (hash/c symbol? variable/c)
            (or/c false/c (recursive-contract current-variables/c #:chaperone))))
(define argument/c (struct/c argument symbol? type/c))
(define current-arguments/c (hash/c symbol? argument/c))
(define need-only-value-set/c (mutable-set/c var-symbol/c))
(define dx-names-hash-set/c (hash/c string? boolean?))
;; NOTE: Reference to a array or a func-return-value, variables are mapped to singleton arrays
(define backrun-ref/c (list/c (or/c 'array-ref 'func-ref) var-symbol/c integer?))
(define functions-symbols/c (hash/c string? symbol?))

(struct func-context
        (.current-variables
         .function-name
         .dx-names-hash-set
         .der-table
         .real-vars-table
         .need-only-value-set
         .need-derivatives-table
         .function-return-value
         (.function-return-type #:mutable)
         .current-arguments)
        #:prefab)

(define func-context/c
  (struct/c func-context
            current-variables/c
            symbol?
            dx-names-hash-set/c
            der-table/c
            real-vars-table/c
            need-only-value-set/c
            need-derivatives-table/c
            symbol?
            type/c
            current-arguments/c))

(struct func-info
    (func-return-symbol
     (arity #:mutable)
     args-type 
     output-base-type 
     output-range)
    #:prefab)

(define func-info/c 
  (struct/c func-info
            symbol?
            integer?
            (listof (cons/c base-type/c (or/c integer? boolean?)))
            base-type/c
            (or/c integer? boolean?)))

(define funcs-info/c (hash/c var-symbol/c func-info/c))

(struct func-call-info
    (.name
     .type
     .arg-list)
    #:prefab)

(define func-call-info/c
  (struct/c func-call-info
            string?
            landau-type/c
            (listof (syntax/c any/c))))

(define func-call-ht/c (hash/c symbol? func-call-info/c))
(define binding-type/c (or/c 'constant 'parameter 'variable 'function))
(define prohibition-list/c (listof binding-type/c))
;; der-table keys, with grouped array cells 
; (define grouped-by-names/c (hash/c var-symbol/c (or/c (listof integer?) (list/c boolean?))))

(struct derivatives-info 
  (.der-table
    .need-only-value-set
    .need-derivatives-table)
  #:prefab)

(define derivatives-info/c
  (struct/c derivatives-info
            der-table/c
            need-only-value-set/c
            need-derivatives-table/c))

;; NOTE: functions

(define (new-variables-nesting vars)
  (variables (make-hash) vars))

(define (search-variable name vars)
  (cond
    ((not vars) #f)
    ((hash-has-key? (variables-current-level vars) name)
     (hash-ref (variables-current-level vars) name))
    ((variables-next-level vars)
     (search-variable name (variables-next-level vars)))
    (else #f)))

(define (search-argument name args)
  (hash-ref args name (lambda () #f)))

(define/contract
  (add-variable! vars name type)
  (-> current-variables/c (syntax/c symbol?) type/c
      symbol?)
  (let* ((name_ (syntax->datum name))
         (sym (gensym name_))
         (_src-pos (syntax-position name))
         (src-pos (if _src-pos _src-pos 0)))
        ; (println (format "add-variable! ~a ~a" name src-pos))
        (hash-set!
         (variables-current-level vars) name_ (variable sym type src-pos))
        sym))

(define (add-argument! args name type)
  (let ((sym (gensym name)))
    (hash-set!
     args name (argument sym type))
    sym))

(define (parse-type-to-syntax t)
  (match t
    ((list 'type (list 'basic-type basic))
     #`(list #,(string->symbol basic) '()))
    ((list 'type (list 'array-type (list 'basic-type basic)
                             "[" size "]"))
     (if (equal? basic "real")
         #`(list #,(string->symbol basic) #,size)
         (raise-syntax-error #f "int arrays are not supported yet" t)))))

(define (parse-type t)
  (match t
    ((list 'type (list 'basic-type basic))
     (list (string->symbol basic) '()))
    ((list 'type (list 'array-type (list 'basic-type basic)
                             "[" size "]"))
     (list (string->symbol basic) (list size)))))

;; NOTE: usage: 'real is changed to 'dual-l after der-annot
(define (change-arg-type! args name new-type)
  (let* ((arg (search-argument name args))
         (arg-sym
          (cond
            (arg (argument-symbol arg))
            (error (format "Bug: no argument with name: ~v" name)))))
   (hash-set! args name (argument arg-sym new-type))))

(define (throw-if-not-int stx [annotation ""])
  (let ([type (syntax-property stx 'landau-type)])
    (unless
        (equal? 'int type)
      (raise-syntax-error #f (format "expect 'int, given: ~a. ~a" type annotation) stx)))
    )

(define (throw-if-not-type expected-type stx [annotation ""])
  (let ([type (syntax-property stx 'landau-type)])
    (unless
        (equal? expected-type type)
      (raise-syntax-error #f (format "expect ~a, given: ~a. ~a" expected-type type annotation) stx))))

(define/contract
  (make-landau-type base-type size)
  (-> base-type/c (or/c integer? boolean?)
      landau-type/c)
  (if size
    (list base-type (list size))
    (list base-type (list))))


(define (any-number? n)
  (if (target-extfloat? TARGET)
    (or (extflonum? n) (number? n))
    (number? n)))

(define/contract
  BUILT-IN-FUNCTIONS
  (hash/c string? (listof base-type/c))
  (make-hash
   (list
    (cons "sqrt" (list 'real))
    (cons "sqr"  (list 'real))
    (cons "sin"  (list 'real))
    (cons "cos"  (list 'real))
    (cons "pow"  (list 'real 'real)))))

(define/contract
  (add-real-var! name type stx real-vars-table)
  (-> (syntax/c symbol?) type/c (syntax/c any/c) real-vars-table/c
      void?)
  (let* ((src-pos_ (syntax-position name))
         (src-pos (if src-pos_ src-pos_ 0))
         (basic-type (car type))
         (name-vs (var-symbol
                   (symbol->string (syntax->datum name))
                   src-pos))
         (range (cadr type))
         (expanded-range (if (equal? range '())
                           #f
                             ;; NOTE: var-size is something like ''5
                             ;; TODO: refactor to use `expand-type`
                            (cadr (cadr (syntax->datum (local-expand (datum->syntax stx range) 'expression '())))))))
        
        (when (equal? basic-type 'real)
          (hash-set! real-vars-table
                     name-vs
                     expanded-range))))


(define/contract
  (ref-to-key ref)
  (-> backrun-ref/c df-name/c)
  (match ref
    ((list 'array-ref var-symb idx) (list var-symb idx))
    ((list 'basic-ref var-symb) var-symb)
    (else (error "bug: unknown ref"))))


(define/contract
  (make-var-mappings-name var-symb dx-name-str)
  (-> var-symbol/c string?
    symbol?)
  (string->symbol (format "d~ad~a_mpg" (var-symbol->string var-symb) dx-name-str)))

(define/contract
  (make-var-der-name var-der-symb dx-name-str)
  (-> var-symbol/c string?
    symbol?)
  (string->symbol (format "d~ad~a_der" (var-symbol->string var-der-symb) dx-name-str)))

(define/contract
  (make-dx-idxs-mappings-name var-symb dx-name-str)
  (-> var-symbol? string?
    symbol?)
  (string->symbol (format "d~ad~a_inv_mpg" (var-symbol->string var-symb) dx-name-str)))
