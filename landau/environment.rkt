#lang racket/base

(require 
  racket/base
  file/md5
  racket/contract/base
  racket/contract/combinator
  racket/contract/region
  racket/extflonum
  racket/flonum
  racket/fixnum
  racket/function
  racket/list
  racket/match
  racket/serialize
  racket/set
  racket/syntax
  "target-config.rkt")

(provide (all-defined-out))

(define constants (make-hash))
(define parameters (make-hash))
(define parameters/c (hash/c symbol? list?))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NOTE: Utils
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (syntax->string stx)
  (symbol->string (syntax->datum stx)))

(define (hash->string ht)
  (for/fold ((str ""))
            (((k v) (in-hash ht)))
             (string-append str (format "~a: ~a\n" k v))))

(define (deepcopy x)
  (deserialize (serialize x)))

(define (raise-error stx msg)
  (raise-syntax-error #f msg stx))

(define (atom-number stx)
  (let ((v (syntax->datum stx)))
    (match v
      ((? number?) v)
      ((? extflonum?) v)
      ((list 'quote (? number?)) (cadr v))
      ((list 'quote (? extflonum?)) (cadr v))
      (else #f))))

(define gensym_
  (let ([counter 0])
    (lambda ([x 'g])
      (if (number? x)
        (set! counter x)
        (begin0 (string->unreadable-symbol
                 (format "~a~a" x counter))
          (set! counter (add1 counter)))))))

(define (syntax->hash-key stx)
 (md5 (format "~a" (syntax->datum stx))))

(define local-expand-memo
  (let ((memo (make-hash)))
    (lambda (stx context-v stop-ids #:reset-memo (reset-memo #f))
      (when reset-memo
        (hash-clear! memo))
      (define key (syntax->hash-key stx))
      ;; #RRstx
      ;; #RR(hash-has-key? memo key)
      ;; #RRkey
      ;; #RRmemo
      ;; FIXME uncomment `if` expression and delete this
      (with-syntax ((value (local-expand stx context-v stop-ids)))
        (hash-set! memo key #'value)
        #'value)
      #;(if (hash-has-key? memo key)
        (hash-ref memo key)
        (with-syntax ((value (local-expand stx context-v stop-ids)))
          (hash-set! memo key #'value)
          #'value)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NOTE: Structs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NOTE: types and interfaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; NOTE: aka box. Box becomes immutable after passing them as synax-parameters, 
; which is not what we want. Hash tables do well 
(define (state/c a) (hash/c fixnum? a))

(define/contract (make-state init-value)  
                 (-> any/c (state/c any/c))
                 (make-hash (list (cons 0 init-value))))

(define/contract (read-state st)
                 (-> (state/c any/c) any/c)
                 (hash-ref st 0))

(define/contract (write-state! st value)
                 (-> (state/c any/c) any/c void?)
                 (hash-set! st 0 value))

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
                 (hash-ref dtbl df-name (thunk (error (format "bug: der-table has no ~a key" df-name)))))


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
                 (integer? (rvt-get-var-size rvt var)))


(define dx-names-set/c (hash/c string? boolean?))
(define grouped-keys-table/c (hash/c var-symbol/c (listof integer?)))


(define dx-mapping-sizes-ht/c (hash/c string? mapping-sizes/c))
(define/contract (dms-member? dms dx-name)
                 (-> dx-mapping-sizes-ht/c string?
                     boolean?)
                 (hash-has-key? dms dx-name))

(define/contract (dms-get-mapping-sizes dms dx-name)
                 (-> dx-mapping-sizes-ht/c string?
                     mapping-sizes/c)
                 (hash-ref dms dx-name (thunk (error (format "bug: dx-mapping-sizes-ht has no ~a key" dx-name)))))


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
                 (unless (var-symbol/c df-name)
                   (error "181"))
                 (cond
                   [(ndt-member? ndt df-name)
                    (let ((dx-name->dx-sizes (ndt-get-dx-mapping-sizes-ht ndt df-name)))
                      (cond
                        ((dms-member? dx-name->dx-sizes dx-name)
                         (dms-get-mapping-sizes dx-name->dx-sizes dx-name))
                        (else #f)))]
                   [else #f]))


(define base-type/c (one-of/c 'real 'int 'dual-l 'dual-r 'int-index))
(define landau-type/c (list/c base-type/c (or/c (list/c integer?) (list/c))))
(define/contract
  (get-type-base type)
  (-> landau-type/c base-type/c)
  (car type))

(define/contract
  (get-type-range type)
  (-> landau-type/c (or/c (list/c integer?) (list/c)))
  (cadr type)) 


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
;; NOTE: Reference to a array or a func-return-value
; nonarray variables are mapped to singleton arrays
(define backrun-ref/c (list/c (or/c 'array-ref 'func-ref) var-symbol/c integer?))
(define functions-symbols/c (hash/c string? symbol?))

(struct mappings-info (.have-minus-ones?) #:prefab)
(define mappings-info/c (struct/c mappings-info boolean?))
(define mappings-table/c (hash/c symbol? (or/c mappings-info/c #f)))
(define/contract (new-mappings-table) (-> mappings-table/c) (make-hash))
(define/contract (get-mappings mappings-table var)
                 (-> mappings-table/c symbol? (or/c mappings-info/c #f))
  (if (hash-has-key? mappings-table var)
    (hash-ref mappings-table var)
    #f))
(define/contract (set-mappings! mappings-table var mappings-info)
                 (-> mappings-table/c symbol? mappings-info/c void?)
  (hash-set! mappings-table var mappings-info))
(define mappings-table-copy hash-copy)
(define (check-if-all-derivatives-are-used ctx dx-size dx-mapping-size mappings-symbol)
  (match mappings-symbol
    (#f #f)
    (_ (let* ((maybe-mappings-info
           (get-mappings (func-context-.mappings-table ctx)
                         mappings-symbol))
         (have-no-minus-ones (if maybe-mappings-info
                                (not (mappings-info-.have-minus-ones? maybe-mappings-info))
                                #f)))
    (and have-no-minus-ones (equal? dx-size dx-mapping-size))
    ))))


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
         .current-arguments
         .self-sufficent-function?
         .mappings-table)
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
            current-arguments/c
            boolean?
            mappings-table/c))

(define/contract 
  (func-context-copy ctx)
  (-> func-context/c func-context/c)
  (func-context
    (deepcopy (func-context-.current-variables ctx))
    (deepcopy (func-context-.function-name ctx))
    (deepcopy (func-context-.dx-names-hash-set ctx))
    (deepcopy (func-context-.der-table ctx))
    (deepcopy (func-context-.real-vars-table ctx))
    (set-copy (func-context-.need-only-value-set ctx))
    (deepcopy (func-context-.need-derivatives-table ctx))
    (deepcopy (func-context-.function-return-value ctx))
    (deepcopy (func-context-.function-return-type ctx))
    (deepcopy (func-context-.current-arguments ctx))
    (deepcopy (func-context-.self-sufficent-function? ctx))
    (mappings-table-copy (func-context-.mappings-table ctx))))

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

(struct func-call-info-pair
  (.func-ret-symbol
   .func-call-info)
  #:prefab)

(define func-call-info-pair/c
  (struct/c func-call-info-pair
            symbol?
            func-call-info/c))

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

(struct function-inline-template
  (.actions-list
   .parameters
   .func-vs)
  #:prefab)

(define function-inline-template/c
  (struct/c function-inline-template
            (listof any/c)
            (listof (list/c var-symbol/c landau-type/c))
            var-symbol/c))

(struct function-inline-semantics-template
  (.name
   .body
   .type
   .parameters
   .return-symbol)
  #:prefab)

(define function-inline-semantics-template/c
  (struct/c function-inline-semantics-template
            (syntax/c symbol?)
            syntax?
            landau-type/c
            (listof symbol?)
            symbol?))

(struct getter-info (type ix-info) #:prefab)
(define getter-info/c
  (struct/c
    getter-info
    (or/c 'cell 'slice 'var 'array)
    (or/c (syntax/c any/c) false/c)))
(define (getter-is-slice? getter-info) (equal? (getter-info-type getter-info) 'slice))
(define (getter-is-cell? getter-info) (equal? (getter-info-type getter-info) 'cell))
(define (getter-is-var? getter-info) (equal? (getter-info-type getter-info) 'var))
(define (getter-is-array? getter-info) (equal? (getter-info-type getter-info) 'array))


(define/contract
  (make-getter-info index slice declated-type)
  (-> (syntax/c any/c) (syntax/c any/c) landau-type/c 
      getter-info/c)
  (when (and (syntax->datum index) (syntax->datum slice))
    (error (format "bug: getter-type: index and slice has non'#f values: ~a ~a" 
                   (syntax->datum index) 
                   (syntax->datum slice))))
  (cond ;; TODO proper getter checks here
    ;; TODO fix all usage
    ((syntax->datum index) 
     (getter-info 'cell index))

    ((syntax->datum slice) 
     (getter-info 'slice #f))

    ((not (empty? (get-type-range declated-type)))
     (getter-info 'array #f))

    (else 
      (getter-info 'var #f)))) 


(struct func-arg-binding
  (.arg-type
   .ref)
  #:prefab)

(define func-arg-binding/c
  (struct/c func-arg-binding
            (symbols 'var 'array 'cell)
            backrun-ref/c))

(define (atom-number/c x)
  (match x
    ((? number?) #t)
    ((? extflonum?) #t)
    (_ #f)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NOTE: functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (new-variables-nesting vars)
  (variables (make-hash) vars))

(define HASH-HIT (make-hash))

(define (search-variable name vars)
  (cond
    ((not vars) #f)
    ((hash-has-key? (variables-current-level vars) name)
     (begin
       (hash-update! HASH-HIT 'hit (lambda (x) (add1 x)) 0)
       (hash-ref (variables-current-level vars) name)))
    ((variables-next-level vars)
     (begin
       (hash-update! HASH-HIT 'miss (lambda (x) (add1 x)) 0)
       (search-variable name (variables-next-level vars))))
    (else #f)))

(define (search-argument name args)
  (hash-ref args name (lambda () #f)))

(define/contract
  (add-variable! vars name type)
  (-> current-variables/c (syntax/c symbol?) type/c
      symbol?)
  (let* ((name_ (syntax->datum name))
         (sym (gensym_ name_))
         (_src-pos (syntax-position name))
         (src-pos (if _src-pos _src-pos 0)))
        ; (println (format "add-variable! ~a ~a" name src-pos))
        (hash-set!
         (variables-current-level vars) name_ (variable sym type src-pos))
        sym))

(define (add-argument! args name type)
  (let ((sym (gensym_ name)))
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
         #`(list #,(string->symbol basic) (list #,size))
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
      (raise-syntax-error #f (format "expect 'int, given: ~a. ~a" type annotation) stx))))

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

(define/contract (to-landau-type stx type-like-anything)
    (-> syntax? any/c landau-type/c)
    (define cant-cast-msg "bug: can't cast ~a to landau-type")
    (define cant-cast-first-item 
      "bug: can't cast ~a to landau-type. Can't cast the first item of a pair to the basic-type/c")
    (define cant-cast-second-item 
      "bug: can't cast ~a to landau-type. Expanded second item is not atom-number")
    (match type-like-anything
      ((? landau-type/c type-like-anything) 
        type-like-anything)
      

      ((? syntax? type-like-anything)
        (to-landau-type stx (syntax->datum type-like-anything)))


      ((list 'list base-type (list 'list range-tree))
       (to-landau-type stx (list base-type (list range-tree))))
      

      ((list 'list base-type (list 'list))
       (to-landau-type stx (list base-type (list))))


      ((list 'list base-type (list 'quote '()))
       (to-landau-type stx (list base-type (list))))
      
      ((? base-type/c type-like-anything)
       (list type-like-anything (list)))

      ((list base-type (list))
       (let* ((base-type (if (syntax? base-type)
                          (syntax->datum base-type)
                          base-type)))
         (unless (base-type/c base-type)
           (raise-error stx (format cant-cast-first-item 
                                    type-like-anything)))
         (list base-type (list))))


      ((list base-type (list type-range-tree))
       (let* ((type-range (if (syntax? type-range-tree) 
                            type-range-tree
                            (datum->syntax stx type-range-tree)))
              (base-type (if (syntax? base-type)
                           (syntax->datum base-type)
                           base-type)))
         (unless (base-type/c base-type)
           (raise-error stx (format cant-cast-first-item
                                    type-like-anything)))
         (if (atom-number type-range)
           (list base-type (list (atom-number type-range)))
           (let ((expanded-range-stx (local-expand 
                                    type-range 'expression '())))
             (if (atom-number expanded-range-stx)
               (list base-type (list (atom-number expanded-range-stx)))
               (raise-error stx (format cant-cast-second-item
                                        type-like-anything)))))))


      (`(,(? base-type/c base-type) ,(? fixnum? type-range-tree))
       (list base-type (list type-range-tree)))
      

      ((list base-type type-range-tree)
       (let* ((type-range (if (syntax? type-range-tree) 
                            (syntax->datum type-range-tree)
                            type-range-tree))
              (base-type (if (syntax? base-type)
                           (syntax->datum base-type)
                           base-type)))
         (unless (base-type/c base-type)
           (raise-error stx (format cant-cast-first-item
                                    type-like-anything)))
         (to-landau-type stx (list base-type type-range))))


      (_ (raise-syntax-error #f (format cant-cast-msg type-like-anything) stx))
      ))

(define/contract
  (add-real-var! name type stx real-vars-table)
  (-> (syntax/c symbol?) type/c (syntax/c any/c) real-vars-table/c
      void?)
  (define src-pos_ (syntax-position name))
  (define src-pos (if src-pos_ src-pos_ 0))
  (define var-landau-type (to-landau-type stx type))
  (define basic-type (get-type-base var-landau-type))
  (define name-vs (var-symbol
                    (symbol->string (syntax->datum name))
                    src-pos))
  (define range (get-type-range var-landau-type))
  (define expanded-range (if (equal? range '())
                           #f
                           (car range)))
  (when (equal? basic-type 'real)
    (hash-set! real-vars-table
               name-vs
               expanded-range)))


(define/contract
  (ref-to-key ref)
  #| (-> backrun-ref/c df-name/c) |#
  (-> any/c df-name/c)
  (match ref
    ((list 'array-ref var-symb idx) (list var-symb idx))
    (else (error (format "bug: unknown ref type: ~a" ref)))))


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

;; NOTE: function argument temporaty variable name
(define/contract
  (make-norm-arg-name func-name arg-number)
  (-> string? integer? 
    string?)
  (format "norm_arg_~a~a" func-name arg-number))


(define/contract (update-current-variables ctx new-curr-var)
                 (-> func-context/c current-variables/c
                     func-context/c)
 (func-context new-curr-var
               (func-context-.function-name ctx)
               (func-context-.dx-names-hash-set ctx)
               (func-context-.der-table ctx)
               (func-context-.real-vars-table ctx)
               (func-context-.need-only-value-set ctx)
               (func-context-.need-derivatives-table ctx)
               (func-context-.function-return-value ctx)
               (func-context-.function-return-type ctx)
               (func-context-.current-arguments ctx)
               (func-context-.self-sufficent-function? ctx)
               (func-context-.mappings-table ctx)))

(define INLINE-VARIABLE-FORMAT "_inl_~a_~a")

(define(make-inlined-variable-name inlined-function-name local-variable-name)
  (format INLINE-VARIABLE-FORMAT inlined-function-name local-variable-name))


