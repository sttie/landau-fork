#lang racket/base

(require 
  (for-syntax racket/base
              syntax/parse
              file/md5
              racket/syntax
              racket/string
              racket/hash
              racket/flonum
              racket/extflonum
              racket/fixnum
              racket/function
              racket/match
              racket/pretty
              racket/set
              racket/serialize
              racket/list
              racket/contract
              #| macro-debugger/stepper |#
              profile
              "environment.rkt"
              "backrun.rkt"
              "target-config.rkt"
              "metalang.rkt"
              (only-in "combinators.rkt" to-c-func-param c-func-arg-decl)
              (except-in "combinators.rkt" to-c-func-param))
  (except-in "combinators.rkt" to-c-func-param)
  "mappings.rkt"
  "type-utils.rkt"
  "target-config.rkt"
  "runtime-defs.rkt"
  "common-for-syntax.rkt"
  "metalang.rkt"
  racket/stxparam
  racket/flonum
  racket/list
  racket/function
  racket/extflonum
  racket/fixnum)

(provide #%module-begin
          #%datum program expr term factor primary element func constant get-value
         var-decl decl-block assignation number if-expr bool-expr bool-term
         bool-factor single-term expr-body func-body for-expr der-apply der-annot
         al_index_name_symbol discard parameter func-call slice_idx print get-derivative
         ;; NOTE: metalang syntax
         _if _if-stm _begin _for _forever _let _let-int _define-var _define-var-with-func-call _decl-real-func _pure-func-call
         _not _or _and _equal?
         _int+ _int- _int* _int/ _int-neg _int= _int> _int>= _int<= _int<
         _rl+ _rl- _rl* _rl/ _rl-neg
         _exact->inexact
         _rl-vector _vector-ref _int-vector-ref _var-ref _vector-set! _set! _func-call 
         _sin _cos _expt _sqr _sqrt
         _0.0 _1.0 _2.0 _0.5
         _break _nothing _empty-statement _local
         )

(define-for-syntax padding 0)

(begin-for-syntax
  (define-syntax-rule
    (print-decorator pair-name expression)
     (begin
       (displayln (string->symbol (format "~a~a-begin" (build-string padding (lambda (x) #\space))
                                          pair-name)))
       (set! padding (add1 padding))
       (let ((e expression))
         (displayln (string->symbol (format "~a~a-end" (build-string padding (lambda (x) #\space))
                                          pair-name)))
       (set! padding (- padding 1))
         e))))

(define-syntax-parameter func-context #f)
(define-syntax-parameter func-call-box 'func-call-box-not-set)
;; NOTE: dx name in currnt assertion loop
(define-syntax-parameter dx-name-in-current-al #f)
;; NOTE: Expansion modes
(define-syntax-parameter typecheck-mode #f)
(define-syntax-parameter expand-value-only #f)
(define-syntax-parameter get-name-mode #f)
;; NOTE: #t if function is called inside an argument of another function. Needs for the func-call ordering
(define-syntax-parameter func-is-called-inside-argument #f)

(define-for-syntax slice-idx-name-GLOBAL 'slice_idx)
(define-for-syntax df-slice-idx-name-GLOBAL 'df_slice_idx)
(define-for-syntax dx-slice-idx-name-GLOBAL 'dx_slice_idx)
(define-for-syntax func-slice-idx-name-GLOBAL 'func_slice_idx)
(define-for-syntax MODULE-FUNCS-TABLE (make-hash))
(define-for-syntax functions-data-ht-GLOBAL #f)
(define-for-syntax funcs-info-GLOBAL #f)

(define-for-syntax debug-compile #f)
(define-for-syntax (log-debug msg [force-log #f])
  (unless (and (not debug-compile) force-log)
    (displayln msg)))

(define/contract-for-syntax
  (get-dx-parameter-size stx parameters dx-name)
  (-> syntax? parameters/c string? (or/c integer? false/c))
  (define range (get-type-range
                  (to-landau-type stx
                                  ;; FIXME dx-name can be not a parameter but an ordinary variable
                                  ;; need to lookup dx set, which is constructed in backrun.rkt/process-actions-list.rkt
                                  (list 'real (hash-ref parameters (string->symbol dx-name) 1)))))
  (match range
    ((list n) n)
    ((list) #f)))

(define/contract-for-syntax
  (get-der-variable-symbol name dx-name-str stx [throw-error? #t])
  (->* (var-symbol/c
         string?
         (syntax/c any/c))
       (boolean?)
       (or/c (syntax/c symbol?) #f))
  (get-symbol-helper 
    stx 
    name
    dx-name-str 
    make-var-der-name 
    (func-context-.current-variables (syntax-parameter-value #'func-context))
    "derivative" 
    throw-error?))

(define/contract-for-syntax
  (get-mappings-variable-symbol name dx-name-str stx [throw-error? #t])
                 (->* (var-symbol/c
                        string?
                        (syntax/c any/c))
                      (boolean?)
                      (or/c (syntax/c symbol?) #f))
                 (get-symbol-helper 
                   stx 
                   name
                   dx-name-str 
                   make-var-mappings-name
                   (func-context-.current-variables (syntax-parameter-value #'func-context))
                   "mapping" 
                   throw-error?))

(define/contract-for-syntax
  (get-dx-idxs-mappings-variable-symbol name dx-name-str stx [throw-error? #t])
                 (->* (var-symbol/c
                        string?
                        (syntax/c any/c))
                      (boolean?)
                      (or/c (syntax/c symbol?) #f))
                 (get-symbol-helper 
                   stx 
                   name
                   dx-name-str
                   make-dx-idxs-mappings-name
                   (func-context-.current-variables (syntax-parameter-value #'func-context))
                   "inv mapping"
                   throw-error?))

;; TODO: make expansion here and return expanded syntax

(define/contract-for-syntax
  (get-der-variable-dx-range var-name dx-name-str stx)
                 (-> var-symbol/c (or/c string? false/c) (syntax/c any/c)
                     (or/c integer? false/c))
                 (get-der-variable-dx-range-helper 
                   var-name 
                   dx-name-str 
                   (func-context-.need-derivatives-table (syntax-parameter-value #'func-context)) 
                   stx))

(define/contract-for-syntax ;; FIXME
  (get-inv-mapping-period df-name dx-name-str stx)
  (-> var-symbol/c string? (syntax/c any/c)
      (or/c integer? false/c))
  (unless (var-symbol/c df-name)
    (error "143"))
  (let ((ndt (func-context-.need-derivatives-table (syntax-parameter-value #'func-context))))
    (cond
      ((not dx-name-str) #f)
      [(ndt-member? ndt df-name)
       (let ((dx-name->dx-sizes (ndt-get-dx-mapping-sizes-ht ndt df-name)))
         (cond
           ((dms-member? dx-name->dx-sizes dx-name-str)
            (mapping-sizes-.inv-mapping-period (dms-get-mapping-sizes dx-name->dx-sizes dx-name-str)))
           (else #f)
           ))]
      [else #f])))


(define-for-syntax (get-real-arg-or-var-type d-name)
  (let
      ((d-name-arg (search-argument
                    (syntax->datum d-name) (func-context-.current-arguments (syntax-parameter-value #'func-context)))))
    (cond
      (d-name-arg
       (let ((base-type (car (argument-type d-name-arg))))
         (if (or (equal? 'real base-type) (equal? 'dual-l base-type))
             (argument-type d-name-arg)
             (raise-syntax-error #f "derivative of non-real in prohibited" d-name))))
        
      (else
       (let 
           ((d-name-var (search-variable
                         (syntax->datum d-name) (func-context-.current-variables (syntax-parameter-value #'func-context)))))
         (cond
           (d-name-var
            (let ((base-type (car (variable-type d-name-var))))
              (if (or (equal? 'real base-type) (equal? 'dual-l base-type))
                  (variable-type d-name-var)
                  (raise-syntax-error #f "derivative of non-real in prohibited" d-name)))) 
           (else (raise-syntax-error #f "name not found; only real function parameters or variables can be used in annotation" d-name))))))))

(define-for-syntax (get-real-arg-or-parameter-type name)
  (let*
      ((name_ (syntax->datum name))
       (type_
        (let ((arg (search-argument
                    name_ (func-context-.current-arguments (syntax-parameter-value #'func-context)))))
          (cond
            (arg
             (argument-type arg))
            ((equal? (func-context-.function-name (syntax-parameter-value #'func-context)) (syntax->datum name))
             (func-context-.function-return-type (syntax-parameter-value #'func-context)))
            ((hash-has-key? parameters name_) (list 'real (hash-ref parameters name_)))
            (else
             (raise-syntax-error #f "name not found" name))))))
    type_))

(define-for-syntax (get-arg-type d-name)
  (get-arg-type-helper d-name (func-context-.current-arguments (syntax-parameter-value #'func-context))))


#| (begin-for-syntax |#
#|  (define/contract (mappend body) |#
#|                   (-> (syntax/c any/c) |#
#|                       (syntax/c any/c)) |#
#|   (match (target-lang TARGET) |#
#|     ('racket #`(begin #,@body)) |#
#|     ;; FIXME: 1st element of list is <void> |#
#|     ('ansi-c #`(apply string-append (filter string? (list #,@body)))))) |#

#|  ;; NOTE: concatinate actions for Racket and C backend |#
#|  (define mbegin |#
#|    (match (target-lang TARGET) |#
#|      ('racket 'begin) |#
#|      ('ansi-c 'string-append)))) |#

(define-for-syntax (write-source stx src-synt)
  (match (target-lang TARGET)
    ('ansi-c
     (quasisyntax/loc stx
       (let*
           ((c-src-out-path (target-c-out-path TARGET))
            (out (open-output-file c-src-out-path #:exists 'replace)))
         (display (prepand-headers #,src-synt) out)
         (close-output-port out)
         (displayln (format "Success.\nThe output is written to ~a" c-src-out-path)))))
    ('racket
     (begin
      (displayln "Success.\n.dau program compiled")
      src-synt))))

(define-for-syntax (make-derivative-getter-func/array stx)
  (with-syntax*
        ((full-idx (datum->syntax stx (gensym_ "full_idx")))
         (dx-mapped-size (datum->syntax stx (gensym_ "dx_mapped_size")))
         (al-index (datum->syntax stx (gensym_ "al_idx")))
         (inv-mapping-period (datum->syntax stx (gensym_ "inv_mapping_period")))
         (dx-idxs-mappings (datum->syntax stx (gensym_ "dx_idx_mappings")))
         (der-vec (datum->syntax stx (gensym_ "der_vec")))
         (argnames (list
                    #'full-idx
                    #'dx-mapped-size
                    #'al-index
                    #'inv-mapping-period
                    #'dx-idxs-mappings
                    #'der-vec))
         (argtypes (list
                    (make-landau-type 'int #f)
                    (make-landau-type 'int #f)
                    (make-landau-type 'int #f)
                    (make-landau-type 'int #f)
                    (make-landau-type 'int -1)
                    (make-landau-type 'real -1)))
         (dx-idx (datum->syntax
                  stx
                  `(_int-vector-ref ,#'dx-idxs-mappings
                                     (_int+ (_int* (_var-ref ,#'inv-mapping-period) (_var-ref ,#'full-idx)) (_var-ref ,#'al-index)))))
         (get_dfdx_cell-name (datum->syntax stx 'get_dfdx_cell))
         (get_dfdx_cell-dx-name (datum->syntax stx 'get_dfdx_cell_dx))
         (get-dfdx-cell
          (datum->syntax
           stx
           `(_decl-real-func ,#'get_dfdx_cell-name ,#'argnames ,#'argtypes
                              (_if (_int< (_var-ref ,#'al-index) (_var-ref ,#'inv-mapping-period))
                                   (_if (_int>= ,#'dx-idx 0)
                                        (_vector-ref ,#'der-vec
                                                      (_int+ (_int* (_var-ref ,#'dx-mapped-size) (_var-ref ,#'full-idx)) ,#'dx-idx))
                                        _0.0)
                                   _0.0))))
         (get-dfdx-cell-dx ;; NOTE: if array is jacobian denominator itself, dx[i]/dx[j] = delta_{ij}
          (datum->syntax
           stx
           `(_decl-real-func ,#'get_dfdx_cell-dx-name ,#'argnames ,#'argtypes
                              (_if (_int< (_var-ref ,#'al-index) (_var-ref ,#'inv-mapping-period))
                                   (_if (_int= (_var-ref ,#'al-index) (_var-ref ,#'full-idx))
                                        _1.0
                                        (_if (_int>= ,#'dx-idx 0)
                                             (_vector-ref ,#'der-vec
                                                           (_int+ (_int* (_var-ref ,#'dx-mapped-size) (_var-ref ,#'full-idx)) ,#'dx-idx))
                                             _0.0))
                                   _0.0)))))
        (values #'get-dfdx-cell #'get-dfdx-cell-dx)))

(define-for-syntax (make-derivative-getter-func stx)
  (with-syntax*
        ((dx-mapped-size (datum->syntax stx (gensym_ "dx_mapped_size")))
        (al-index (datum->syntax stx (gensym_ "al_idx")))
        (inv-mapping-period (datum->syntax stx (gensym_ "inv_mapping_period")))
        (dx-idxs-mappings (datum->syntax stx (gensym_ "dx_idx_mappings")))
        (der-vec (datum->syntax stx (gensym_ "der_vec")))
        (argnames (list
                                 #'al-index
                                 #'inv-mapping-period
                                 #'dx-idxs-mappings
                                 #'der-vec))
        (argtypes (list
                                 (make-landau-type 'int #f)
                                 (make-landau-type 'int #f)
                                 (make-landau-type 'int -1)
                                 (make-landau-type 'real -1)))
        (get-dfdx-var-name (datum->syntax stx 'get_dfdx_var))
        (dx-idx (datum->syntax
                 stx
                 `(_if (_int< (_var-ref ,#'al-index) (_var-ref ,#'inv-mapping-period))
                       (_int-vector-ref ,#'dx-idxs-mappings (_var-ref ,#'al-index))
                       -1)))
        (get-dfdx-var ;; NOTE: get derivative of a nonarray variable
         (datum->syntax
          stx
          `(_decl-real-func ,#'get-dfdx-var-name ,#'argnames ,#'argtypes
                            (_if (_int>= ,#'dx-idx 0)
                                 (_vector-ref ,#'der-vec ,#'dx-idx)
                                 _0.0)))))
        #'get-dfdx-var))


(define-syntax (program stx)
  (syntax-parse stx
    [({~literal program} body ...)
     (let*-values
      (((tik) (current-inexact-milliseconds))
       ((functions-data-assocs-list)
       (parameterize
        ([current-namespace (module->namespace (collection-file-path "backrun.rkt" "landau"))])
        (begin
         (eval stx)))))
      (define tok (current-inexact-milliseconds))
      (set! functions-data-ht-GLOBAL (make-hash (cadr functions-data-assocs-list)))
      (set! funcs-info-GLOBAL (car functions-data-assocs-list))
      (let-values
        ;; NOTE: Make getter functions to access derivatives values
        ;; Each backend has it's own getter implementation
       (((get-dfdx/array get-dfdx/array-dx) (make-derivative-getter-func/array stx))
        ((get-dfdx) (make-derivative-getter-func stx)))
       (with-syntax*
         ((get-dfdx/array get-dfdx/array)
          (get-dfdx/array-dx get-dfdx/array-dx)
          (get-dfdx get-dfdx))
         (write-source stx
                        (datum->syntax stx `(_begin ,#'get-dfdx/array ,#'get-dfdx/array-dx ,#'get-dfdx
                                                          (_begin ,@#'(body ...))))))))]))

(define-syntax (expr stx)
  (syntax-parse stx
    ((_ expr1-unexpanded op:plus-or-minus expr2-unexpanded)
     (let* ((expr1 (local-expand-memo #'expr1-unexpanded 'expression '()))
            (expr2 (local-expand-memo #'expr2-unexpanded 'expression '()))
            (type1 (syntax-property expr1 'landau-type))
            (type2 (syntax-property expr2 'landau-type))
            (op (syntax->datum #'op)))
    (let*-values (((expr1 expr2 type) (binary-op-cast expr1 expr2 stx))
                  ((expr1-atom) (atom-number expr1))
                  ((expr2-atom) (atom-number expr2))
                  ((al-index-name_) 'al_index_name_symbol))
      (cond
        ((or (equal? type 'real) (is-slice-of-type 'real type))
         (if (and expr1-atom expr2-atom)
           (is-type_ type
                     (quasisyntax/loc stx #,((if (equal? op "+") rl+ rl-) expr1-atom expr2-atom)))
           (is-type_ type
                     (datum->syntax stx 
                                    `(,(if (equal? op "+") #'_rl+ #'_rl-) ,expr1 ,expr2)))))
        ((equal? type 'int)
         (is-type_ type
                   (if (and expr1-atom expr2-atom)
                     (quasisyntax/loc stx #,((if (equal? op "+") fx+ fx-) expr1-atom expr2-atom))
                     (datum->syntax stx
                                    `(,(if (equal? op "+") #'_int+ #'_int-) ,expr1 ,expr2)))))
        ((or (equal? type 'dual-b) (is-slice-of-type 'dual-b type))
         (with-syntax*
             ((al-index-name-synt (datum->syntax stx al-index-name_))
              ;; FIXME: where is in used?
              (al-index-name-str (symbol->string al-index-name_)))
           (cond
             ((and expr1-atom expr2-atom)
              (raise-syntax-error #f "bug: expr1-atom expr2-atom casted to dual-b" stx))
             ;; NOTE: The case when only one of the types is a bundle
             ((and (or (equal? type1 'dual-b) (is-slice-of-type 'dual-b type1))
                   (or (or (equal? type2 'real)
                           (equal? type2 'int)) 
                       (is-slice-of-type 'real type2)))
              (let-values
                  (((dual-expr_ n) (if expr2-atom
                                     (values #'expr1-unexpanded expr2-atom )
                                     (values #'expr1-unexpanded expr2))))
; (println (format "~a | ~a | ~a" expr1 expr2 type))
; (println "")
                (is-type_ type
                          (with-syntax* ((dual-expr dual-expr_)
                                         (n n)
                                         (expanded-dual-expr (local-expand-memo #'dual-expr 'expression '()))
                                         (dual-b-value (get-value-stx #'expanded-dual-expr))
                                         (dual-b-derivative (get-derivative-stx #'expanded-dual-expr))
                                         (oprtr (if (equal? op "+") #'_rl+ #'_rl-)))
                                         
                            (quasisyntax/loc stx 
                                  (list #,(datum->syntax stx 
                                         `(,#'oprtr
                                           ,#'dual-b-value
                                           ,#'n))
                                        #,#'dual-b-derivative))))))

            ((and (or (or (equal? type1 'real)
                          (equal? type1 'int))
                      (is-slice-of-type 'real type1))
                  (or (equal? type2 'dual-b) (is-slice-of-type 'dual-b type2)))
              (let-values
                  (((n dual-expr_)
                    (if expr1-atom
                      (values expr1-atom #'expr2-unexpanded)
                      (values expr1 #'expr2-unexpanded))))
                (is-type_ type
                          (with-syntax* ((dual-expr dual-expr_)
                                         (n n)
                                         (expanded-dual-expr (local-expand-memo #'dual-expr 'expression '()))
                                         (dual-b-value (get-value-stx #'expanded-dual-expr))
                                         (dual-b-derivative (get-derivative-stx #'expanded-dual-expr))
                                         (oprtr (if (equal? op "+") #'_rl+ #'_rl-)))
                                         
                            (quasisyntax/loc stx 
                                  (list #,(datum->syntax stx
                                         `(,#'oprtr
                                           ,#'n
                                           ,#'dual-b-value))
                                        #,(datum->syntax 
                                           stx
                                           `(_rl-neg ,#'dual-b-derivative))))))))
            
                                
             ((or (and (equal? type1 'dual-b) (equal? type2 'dual-b))
                  (and (equal? type1 'dual-b) (is-slice-of-type 'dual-b type2))
                  (and (is-slice-of-type 'dual-b type1) (equal? type2 'dual-b))
                  (and (is-slice-of-type 'dual-b type1) (is-slice-of-type 'dual-b type2)))
                  
              (with-syntax* ((expanded-dual-expr-1 (local-expand-memo #'expr1-unexpanded 'expression '()))
                             (dual-b-value-1 (get-value-stx #'expanded-dual-expr-1))
                             (dual-b-derivative-1 (get-derivative-stx #'expanded-dual-expr-1))

                             (expanded-dual-expr-2 (local-expand-memo #'expr2-unexpanded 'expression '()))
                             (dual-b-value-2 (get-value-stx #'expanded-dual-expr-2))
                             (dual-b-derivative-2 (get-derivative-stx #'expanded-dual-expr-2))
                             (oprtr (if (equal? op "+") #'_rl+ #'_rl-))
                             (res 
                              (is-type_ type
                                            (quasisyntax/loc stx
                                              (list
                                               #,(datum->syntax 
                                                  stx
                                                  `(,#'oprtr
                                                    ,#'dual-b-value-1
                                                    ,#'dual-b-value-2))
                                               #,(datum->syntax
                                                  stx 
                                                  `(,#'oprtr
                                                    ,#'dual-b-derivative-1
                                                    ,#'dual-b-derivative-2)))))))
                #'res)))))
        (else
         (raise-syntax-error #f (format "unsupported type: ~v" type) stx))))))

    [(_ "(" expr1 ")")
     (begin
      (syntax-track-origin (syntax/loc stx expr1) #'expr1 #'expr))]

    [(_ expr1)
     (begin
      (syntax-track-origin (syntax/loc stx expr1) #'expr1 #'expr))
     ]))
    

(define-for-syntax (binary-op-cast expr1 expr2 stx) (binary-op-cast-helper expr1 expr2 stx))


;; TODO: throw excetion if zero-denomintor
(define-syntax (term stx)
  (syntax-parse stx
    ((_ expr1-unexpanded op:mul-or-div expr2-unexpanded)

     (begin
       (let* ((expr1_ (local-expand-memo #'expr1-unexpanded 'expression '()))
              (expr2_ (local-expand-memo #'expr2-unexpanded 'expression '()))
              (type1 (syntax-property expr1_ 'landau-type))
              (type2 (syntax-property expr2_ 'landau-type))
              (op (syntax->datum #'op)))

         (let*-values (((expr1 expr2 type) (binary-op-cast expr1_ expr2_ stx))
                      ((expr1-atom) (atom-number expr1))
                      ((expr2-atom) (atom-number expr2))
                      ((al-index-name_) 'al_index_name_symbol))
          (with-syntax ((expr2-stx expr2_))
           (cond
             ((or (equal? type 'real) (is-slice-of-type 'real type))
              (begin
               (if (and expr1-atom expr2-atom)
                  (is-type_ type
                     (quasisyntax/loc stx #,((if (equal? op "*") rl* rl/) expr1-atom expr2-atom)))
                  (is-type_ type
                         (datum->syntax stx 
                            `(,(if (equal? op "*") #'_rl* #'_rl/) ,expr1 ,expr2))))))
             ;; NOTE: Threre are no int slices yet
             ((equal? type 'int)
              (is-int
               (begin
                (if (and expr1-atom expr2-atom)
                (quasisyntax/loc
                 stx
                 #,((if (equal? op "*") fx* fxquotient) expr1-atom expr2-atom))
                (datum->syntax stx
                 `(,(if (equal? op "*") #'_int* #'_int/) ,expr1 ,expr2))))))
             ((or (equal? type 'dual-b) (is-slice-of-type 'dual-b type))
              (begin
               (cond
                 ((and expr1-atom expr2-atom)
                  (raise-syntax-error #f "bug: expr1-atom expr2-atom casted to dual-b" stx))
                 ;; NOTE: The case when only one of the types is a bundle
                 ((xor (or (equal? type1 'dual-b) (is-slice-of-type 'dual-b type1)) 
                       (or (equal? type2 'dual-b) (is-slice-of-type 'dual-b type2)))
                  (let-values
                    (((n_ dual-expr_)
                      (cond
                        ((or (equal? type2 'dual-b) (is-slice-of-type 'dual-b type2))
                         (if expr1-atom
                           (values expr1-atom #'expr2-unexpanded)
                           (values expr1 #'expr2-unexpanded)))
                        (else
                         (if expr2-atom
                           (values expr2-atom #'expr1-unexpanded)
                           (values expr2 #'expr1-unexpanded))))))
                    (is-type_ type
                              (with-syntax* ((dual-expr dual-expr_)
                                             (n n_)
                                             (expanded-dual-expr (local-expand-memo #'dual-expr 'expression '()))
                                             (dual-b-value (get-value-stx #'expanded-dual-expr))
                                             (dual-b-derivative (get-derivative-stx #'expanded-dual-expr)))
                                ;; NOTE: f * n
                                (if (equal? op "*")
                                  
                                  (quasisyntax/loc stx
                                              (list
                                               #,(datum->syntax stx `(_rl* ,#'dual-b-value ,#'n))
                                               #,(datum->syntax stx `(_rl* ,#'dual-b-derivative ,#'n))))
                                  
                                  (if (or (equal? type1 'dual-b) (is-slice-of-type 'dual-b type1))
                                    ;; NOTE: f / n
                                    (begin
                                     
                                     (quasisyntax/loc stx
                                                 (list
                                                  #,(datum->syntax stx `(_rl/ ,#'dual-b-value ,#'n))
                                                  #,(datum->syntax stx `(_rl/ ,#'dual-b-derivative ,#'n)))))
                                    ;; NOTE: n / f
                                    (begin
                                     (quasisyntax/loc stx
                                                 (list
                                                  #,(datum->syntax stx `(_rl/ ,#'n ,#'dual-b-value))
                                                  #,(datum->syntax stx `(_rl* (_rl-neg _1.0) ,#'n (_rl/ _1.0 (_sqr ,#'dual-b-value)) ,#'dual-b-derivative)))))))))))
                 
                 ((or (and (equal? type1 'dual-b) (equal? type2 'dual-b))
                      (and (equal? type1 'dual-b) (is-slice-of-type 'dual-b type2))
                      (and (is-slice-of-type 'dual-b type1) (equal? type2 'dual-b))
                      (and (is-slice-of-type 'dual-b type1) (is-slice-of-type 'dual-b type2)))
                  
                  (with-syntax*
                    ((expanded-dual-expr-1 (local-expand-memo #'expr1-unexpanded 'expression '()))
                     (dual-b-value-1 (get-value-stx #'expanded-dual-expr-1))
                     (dual-b-derivative-1 (get-derivative-stx #'expanded-dual-expr-1))
                     
                     (expanded-dual-expr-2 (local-expand-memo #'expr2-unexpanded 'expression '()))
                     (dual-b-value-2 (get-value-stx #'expanded-dual-expr-2))
                     (dual-b-derivative-2 (get-derivative-stx #'expanded-dual-expr-2))
                     (res 
                      (is-type_ type
                               (has-dual
                                (if (equal? op "*")
                                  (begin
                                   (quasisyntax/loc stx
                                               (list
                                                #,(datum->syntax stx `(_rl*
                                                 ,#'dual-b-value-1
                                                 ,#'dual-b-value-2))
                                                #,(datum->syntax stx `(_rl+ (_rl* ,#'dual-b-value-1 ,#'dual-b-derivative-2)
                                                     (_rl* ,#'dual-b-value-2 ,#'dual-b-derivative-1))))))
                                  
                                  (quasisyntax/loc stx
                                              (list
                                               #,(datum->syntax stx `(_rl/
                                                ,#'dual-b-value-1
                                                ,#'dual-b-value-2))
                                               #,(datum->syntax stx 
                                                 `(_rl/ (_rl- (_rl* ,#'dual-b-derivative-1 ,#'dual-b-value-2)
                                                           (_rl* ,#'dual-b-derivative-2 ,#'dual-b-value-1))
                                                      (_sqr ,#'dual-b-value-2))))))))))
                    #'res)))))))))))
       
    [({~literal term} term1)
     (with-syntax ((expanded (local-expand-memo #'term1 'expression '())))
       (syntax-track-origin (syntax/loc stx expanded) #'term1 #'expr))]))
    
(define-for-syntax (cast-func-args-list func-args-list)
  (for/list ((par (in-list func-args-list)))
    (with-syntax ((par-exp-value-part (extract (local-expand-memo
                                                #`(syntax-parameterize
                                                      ((expand-value-only #t))
                                                    par) 'expression '()))))
      #'par-exp-value-part
      )))

(define/contract-for-syntax
  (map-built-in-function stx func-str args arg-types)
  (-> (syntax/c any/c) string? (listof (syntax/c any/c)) (listof (one-of/c 'int 'real))
      (syntax/c any/c))
  (with-syntax* ((maybe-casted
                   (for/list ((arg (in-list args))
                              (arg-type (in-list arg-types)))
                     (let ((arg-atom (atom-number arg)))
                       (if (equal? arg-type 'int)
                         (if arg-atom 
                           (datum->syntax stx `(->rl ,arg-atom))
                           (datum->syntax stx `(_exact->inexact ,arg)))
                         arg))))
                 (all-atoms (andmap atom-number (syntax-e #'maybe-casted)))
                 (atoms-list (map atom-number (syntax-e #'maybe-casted))))
                (match func-str
                  ("sqr"
                   (is-type_ 'real
                             (if (syntax->datum #'all-atoms)
                               (quasisyntax/loc stx #,(rl-sqr (syntax->datum #'all-atoms)))
                               (datum->syntax stx
                                              `(_sqr ,@#'maybe-casted)))))
                  ("sqrt"
                   (is-type_ 'real
                             (if (syntax->datum #'all-atoms)
                               (quasisyntax/loc stx #,(rl-sqrt (syntax->datum #'all-atoms)))
                               (datum->syntax stx
                                              ;; FIXME: if negative?
                                              `(_sqrt ,@#'maybe-casted)))))
                  ("sin"
                   (is-type_ 'real
                             (if (syntax->datum #'all-atoms)
                               (begin
                                 (quasisyntax/loc stx #,(rl-sin (syntax->datum #'all-atoms))))
                               (datum->syntax stx
                                              `(_sin ,@#'maybe-casted)))))
                  ("cos"
                   (is-type_ 'real
                             (if (syntax->datum #'all-atoms)
                               (quasisyntax/loc stx #,(rl-cos (syntax->datum #'all-atoms)))
                               (datum->syntax stx
                                              `(_cos ,@#'maybe-casted)))))
                  ("pow"
                   (is-type_ 'real
                             (if (syntax->datum #'all-atoms)
                               (quasisyntax/loc stx #,(rl-expt (car (syntax->datum #'atoms-list)) (cadr (syntax->datum #'atoms-list))))
                               (datum->syntax stx
                                              `(_expt ,@#'maybe-casted))))))))


(define/contract-for-syntax
  (bind-parameters-to-arguments inlined-function-name template bindings)
  (-> string? syntax? (hash/c symbol? (or/c symbol? atom-number/c))
      syntax?)
  (define (rebind-name _inlined-function-name _bindings _name)
    (define new-binding 
      (hash-ref _bindings 
                (syntax->datum _name)
                #f))
    (match new-binding
      ;; NOTE: if there is no such name in bindings, then
      ; the name is either a constant or a local variable's name and it should be
      ; prepanded with _ 
      (#f (cond 
            ((hash-has-key? constants (syntax->datum _name)) _name)
            (else (format-id _name
                             INLINE-VARIABLE-FORMAT _inlined-function-name _name 
                             #:source _name #:props _name))))
      ((? atom-number/c new-binding) new-binding)
      (_ (with-syntax ((new-binding-stx new-binding))
           (syntax-track-origin (syntax/loc _name new-binding-stx) _name #'expr)))))

  (define inlined-function 
    (syntax-parse template
      (get-value:get-value-cls
        (let ((new-binding (rebind-name inlined-function-name bindings #'get-value.name)))
          (datum->syntax 
            template
            ;; NOTE: if new-binding is a number then get-value is transformed to
            ; (number new-binding). Otherwise, new-binding is a variable name 
            ; and new get-value is genetated.
            (cond 
              ((atom-number/c new-binding)
               (begin
                 #| (displayln "FIXME: string->number is rounding precision to fit 64bit") |#
                 `(number ,(string->number (to-string new-binding)))))
              (else 
                (with-syntax ((new-name new-binding))
                  (match (map syntax->datum 
                              (list 
                                #'get-value.index 
                                #'get-value.slice-colon 
                                #'get-value.index-start 
                                #'get-value.index-end))
                    ((list #f #f #f #f)
                     `(get-value ,#'new-name))
                    ((list _ #f #f #f)
                     `(get-value ,#'new-name "[" ,#'get-value.index "]"))
                    ((list #f _ #f #f)
                     `(get-value ,#'new-name "[" ,#'get-value.slice-colon "]"))
                    ((list #f _ _ #f)
                     `(get-value ,#'new-name "[" ,#'get-value.index-start ,#'get-value.slice-colon "]"))
                    ((list #f _ #f _)
                     `(get-value ,#'new-name "[" ,#'get-value.slice-colon ,#'get-value.index-end "]"))
                    ((list #f _ _ _)
                     `(get-value ,#'new-name "[" 
                                 ,#'get-value.index-start 
                                 ,#'get-value.slice-colon 
                                 ,#'get-value.index-end "]"))))))))) 

      (((~literal assignation) name:id
                               getter:getter-cls
                               ((~literal mut-assign) op) value:expr)
       (let* ((op_ (syntax->datum #'op))
              (getter-for-splice (syntax-e #'getter))
              (binop (match op_
                       ("+=" "+")
                       ("-=" "-")
                       ("*=" "*")
                       ("/=" "/"))))
         (with-syntax ((new-name (rebind-name inlined-function-name bindings #'name)))
           (match op_
             ((or "+=" "-=")
              (datum->syntax template 
                             `(assignation ,#'new-name ,@getter-for-splice "="
                                           (expr
                                             (expr 
                                               (term 
                                                 (factor 
                                                   (primary 
                                                     (element 
                                                       (get-value ,#'new-name ,@getter-for-splice))))))
                                             ,binop 
                                             ,(bind-parameters-to-arguments inlined-function-name #'value bindings)))))
             ((or "*=" "/=")
              (datum->syntax template 
                             `(assignation ,#'new-name ,@getter-for-splice "="
                                           (expr
                                             (term
                                               (term 
                                                 (factor 
                                                   (primary 
                                                     (element 
                                                       (get-value ,#'new-name ,@getter-for-splice)))))
                                               ,binop
                                               ,(bind-parameters-to-arguments inlined-function-name #'value bindings))))))))))


      (((~literal assignation) name:id 
                               getter:getter-cls
                               "=" value:expr)
       (let ((getter-for-splice (syntax-e #'getter)))
         (with-syntax ((new-name (rebind-name inlined-function-name bindings #'name)))
           (datum->syntax template 
                          ;; NOTE: assignation to parameters is prohibited.
                          ; Left-hand side name is either function name or function local variable.
                          ; It should be renamed.
                          `(assignation ,#'new-name ,@getter-for-splice
                                        "=" ,(bind-parameters-to-arguments inlined-function-name #'value bindings))))))
      (((~literal der-annot) _ "'" _ "=" _)
       (raise-syntax-error #f "Function with derivatives annotation was called from the main function. This feature is not implemented yet." template))

      (((~literal der-apply) _ "=" _ "'" _)
       (raise-syntax-error #f "Function with derivatives application was called from
                           the main function. This feature is not implemented yet." template))

      (((~literal discard) "discard" _ "'" _)
       (raise-syntax-error #f "Function with derivatives discard was called from
                           the main function. This feature is not implemented yet." template))

      (({~literal var-decl}
                  ((~literal type) ((~literal array-type) basic-type "[" num "]")) name:id (~seq "=" value:expr)) 
       (datum->syntax template 
                      `(var-decl 
                         (type 
                           (array-type 
                             ,#'basic-type "[" ,#'num "]"))
                         ;; NOTE: inlined function's local variables are prepanded with _
                         ,(format-id #'name INLINE-VARIABLE-FORMAT inlined-function-name #'name #:source #'name #:props #'name) 
                         "=" ,(bind-parameters-to-arguments inlined-function-name #'value bindings))))

      (({~literal var-decl} type:expr name:id)
       (datum->syntax template `(var-decl ,#'type
                                          ,(format-id #'name INLINE-VARIABLE-FORMAT inlined-function-name #'name #:source #'name #:props #'name))))

      (({~literal var-decl} type:expr name:id "=" value:expr)
       (datum->syntax template `(var-decl ,#'type
                                          ,(format-id #'name INLINE-VARIABLE-FORMAT
                                                      inlined-function-name #'name #:source #'name #:props #'name) 
                                          "=" ,(bind-parameters-to-arguments inlined-function-name #'value bindings))))

      (((~literal par) expr)
       (datum->syntax template `(par ,(bind-parameters-to-arguments inlined-function-name #'expr bindings))))

      (((~literal other-par) "," expr)
       (datum->syntax template `(other-par "," ,(bind-parameters-to-arguments inlined-function-name #'expr bindings))))

      (((~literal func-call) function-name "(" ({~literal parlist} par*:par-spec ...) ")")
       (begin
         (datum->syntax template 
                        `(func-call 
                           ,#'function-name 
                           "(" 
                           (parlist ,@(bind-parameters-to-arguments inlined-function-name #'(par* ...) bindings))
                           ")"))))

      (((~literal func-body) . children-pat)
       ;; NOTE: func-body is inlined and should be in the caller namespace 
       ;; expr-body introduce the new scope layer, so function local variables
       ; do not conflicts with caller variables.
       (with-syntax 
         ((binded-stx 
            (for/list ((child (in-list (syntax-e #'children-pat))))
              (bind-parameters-to-arguments inlined-function-name child bindings)))) 
         #`(expr-body #,@#'binded-stx)))

      ((parent-pat . children-pat)
       (begin
         (with-syntax 
           ((binded-stx 
              (for/list ((child (in-list (syntax-e #'children-pat))))
                (bind-parameters-to-arguments inlined-function-name child bindings)))) 
           #`(#,(bind-parameters-to-arguments inlined-function-name #'parent-pat bindings) #,@#'binded-stx))))

                                                                      (x #'x)))
                                                  inlined-function)


;; TODO:
;; map functions to slices
(define-syntax (func-call stx)
  (syntax-parse stx
    ((_ function-name "(" ({~literal parlist} par*:par-spec ...) ")")
     (let* ((func-str (symbol->string (syntax->datum #'function-name)))
            (fake-src-pos 0)
            (func-vs (var-symbol func-str fake-src-pos))
            (func-pars-list (for/list ((par-value (in-list (syntax-e #'(par*.par-value ...)))))
                                      par-value))
            (par-list-len (length func-pars-list)))
       (if (hash-has-key? BUILT-IN-FUNCTIONS func-str)
           (let* ((expected-arg-types (hash-ref BUILT-IN-FUNCTIONS func-str))
                  (expected-arity (length expected-arg-types)))
             (begin
             (unless (equal? par-list-len expected-arity)
               (raise-syntax-error
                #f
                (format "function ~a expects ~a parameter, but ~a provided" func-str expected-arity par-list-len)
                stx))
             (match expected-arity
              (1 (with-syntax*
                 ((expr (car func-pars-list))
                  (expanded (local-expand-memo #'expr 'expression '() #:reset-memo #t))
                  (expanded-type (syntax-property #'expanded 'landau-type)))
               (let ((expanded-type-datum (syntax->datum #'expanded-type)))
                 (match expanded-type-datum
                   ('dual-b
                    (with-syntax*
                      ((dual-b-value (get-value-stx #'expanded))
                       (dual-b-derivative (get-derivative-stx #'expanded))
                       (applied-to-real (map-built-in-function stx func-str (list #'dual-b-value) (list 'real))))
                      (match func-str
                        ("sqr"
                         (is-type_ 'dual-b
                                   (quasisyntax/loc 
                                     stx
                                     (list
                                       #,#'applied-to-real
                                       #,(datum->syntax
                                           stx
                                           `(_rl* (_rl* _2.0 ,#'dual-b-value) ,#'dual-b-derivative))))))
                        ("sqrt"
                         (is-type_ 'dual-b
                                   (quasisyntax/loc 
                                     stx
                                     (list
                                       ;; FIXME: if negative?
                                       #,#'applied-to-real
                                       #,(datum->syntax
                                           stx
                                           `(_rl* _0.5 (_sqrt (_rl/ _1.0 ,#'dual-b-value)) ,#'dual-b-derivative))))))
                        ("sin"
                         (is-type_ 'dual-b
                                   (quasisyntax/loc 
                                     stx
                                     (list
                                       #,#'applied-to-real
                                       #,(datum->syntax
                                           stx
                                           `(_rl* (_cos ,#'dual-b-value) ,#'dual-b-derivative))))))
                        ("cos"
                         (is-type_ 'dual-b
                                   (quasisyntax/loc 
                                     stx
                                     (list
                                       #,#'applied-to-real
                                       #,(datum->syntax
                                           stx
                                           `(_rl* (_rl-neg _1.0) (_sin ,#'dual-b-value) ,#'dual-b-derivative)))))))))


                   ((or 'real 'int)
                    (map-built-in-function stx func-str (list #'expanded) 
                                           (list expanded-type-datum)))
                   (else
                    (raise-syntax-error #f 
                                        (format "functions ~a do not accept given type: ~a" 
                                                (hash-keys BUILT-IN-FUNCTIONS) (syntax->datum #'expanded-type)) 
                                        stx))))))
           (2 (with-syntax*
                 ((expanded-1 (local-expand-memo (car func-pars-list) 'expression '() #:reset-memo #t))
                  (expanded-2 (local-expand-memo (cadr func-pars-list) 'expression '() #:reset-memo #t)))
               (let ((expanded-type-1 (syntax-property #'expanded-1 'landau-type))
                     (expanded-type-2 (syntax-property #'expanded-2 'landau-type)))
                 (let ()
                  (match (list expanded-type-1 expanded-type-2)
                   ((list 'dual-b 'real)
                    (with-syntax*
                      ((n #'expanded-2)
                       (dual-b-value (get-value-stx #'expanded-1))
                       (dual-b-derivative (get-derivative-stx #'expanded-1))
                       (applied-to-real (map-built-in-function stx func-str 
                                                               (list #'dual-b-value #'n) 
                                                               (list 'real 'real))))
                      (match func-str
                        ("pow"
                         (is-type_ 'dual-b
                                   (quasisyntax/loc 
                                     stx
                                     (list
                                       #,#'applied-to-real
                                       #,(datum->syntax
                                           stx
                                           `(_rl* ,#'n (_expt ,#'dual-b-value (_rl- ,#'n _1.0)) ,#'dual-b-derivative)))))))))

                   ((list 'dual-b 'dual-b)
                    (error "not implemented"))

                   ((list 'real 'dual-b)
                    (error "not implemented"))
 
                   ((list 'real 'real)
                    (map-built-in-function stx func-str (list #'expanded-1 #'expanded-2) (list 'real 'real)))

                   (else
                    (raise-syntax-error #f 
                                        (format "bug: function ~a do not accept given types: ~a" 
                                                (list "pow") 
                                                (list expanded-type-1 expanded-type-2)) 
                                        stx))))))))))
           ;; NOTE: user defined functions
           ;; FIXME: Uncomment later
           (begin
             #| (check-func |#
             #|  stx |# 
             #|  func-vs |#
             #|  func-pars-list |# 
             #|  funcs-info-GLOBAL |#
             #|  (func-context-.current-variables (syntax-parameter-value #'func-context)) |#
             #|  (func-context-.current-arguments (syntax-parameter-value #'func-context))) |#
             (when (equal? (target-lang TARGET) 'racket)
               (raise-syntax-error #f "user defined functions are temporary broken in Racket backend. Please, use C backend instead." stx))
             ;; TODO pass the function-name in inlined-function-list to use its syntax-position
             (with-syntax* 
               ((typecheck-nonsense #'(error "Bug: typecheck-nonsense in runtime"))
                (function-return-variable (format-id 
                                            #'function-name 
                                            "~a~a" 
                                            #'function-name 
                                            (syntax-position #'function-name) 
                                            #:source #'function-name #:props #'function-name)))
               (let* ((typecheck-mode-on (syntax-parameter-value #'typecheck-mode))
                    (subfunc-call-box-value (make-state (list)))
                    (func-args-list-expanded 
                      (for/list ((par (in-list func-pars-list)))
                        (with-syntax* 
                          ((par par)
                           (par-exp-value-part
                             (extract 
                               (local-expand-memo
                                 #`(syntax-parameterize
                                     ((typecheck-mode #,typecheck-mode-on)
                                      ;;  NOTE: whenever a child func-call macro is called,
                                      ;         func-call-box will be populated
                                      (func-call-box #,subfunc-call-box-value)
                                      (func-is-called-inside-argument #t))
                                     par) 'expression '()
                                 #:reset-memo #t))))
                          #'par-exp-value-part)))
                    (func-args-list-casted_ #`(list #,@func-args-list-expanded))
                    (func-slice-range (func-info-output-range (hash-ref funcs-info-GLOBAL func-vs)))
                    (func-base-type (func-info-output-base-type (hash-ref funcs-info-GLOBAL func-vs)))
                    ;; NOTE: if a single real argument is dual-b, function return type should be dual-b
                    ; So it should return dual-b. The easiest way to do that is to genetate get-value.
                    ; func-call macro generates get-value syntax.
                    (func-output-base-type (for/fold 
                                             ((func-type func-base-type))
                                             ((arg (in-list func-args-list-expanded)))
                                             (match (list func-type (syntax-property arg 'landau-type))
                                               ((list _ (list 'dual-b _)) 'dual-b)
                                               ((list _ 'dual-b) 'dual-b)
                                               ((list 'dual-b _) 'dual-b)
                                               ((list t _) t))))

                    (func-return-symbol (gensym_ func-str))
                    ;; NOTE: Read func-call-box, possibly populated by child func-call macro
                    (func-call-box-value (syntax-parameter-value #'func-call-box))
                    (func-argument-names (for/list ((arg (in-list func-args-list-expanded))
                                                    (unexpanded-arg (in-list func-pars-list)))
                                           (define arg-name (syntax-property arg 'get-value-name))
                                           ;; FIXME pass atom number as a property
                                           (define arg-atom-number (atom-number arg))
                                           (match (list arg-name arg-atom-number)
                                             ((list #f #f)
                                              (raise-syntax-error 
                                                #f 
                                                "bug: function's argument is neither a variable nor a numerical constant" stx))
                                             ((list #f _) arg-atom-number)
                                             ((list _ #f) arg-name)
                                             ((list _ _)
                                              (raise-syntax-error 
                                                #f 
                                                "bug: function's argument has a name and resolved 
                                                to a numerical constant" stx))))))

               (define func-template (hash-ref MODULE-FUNCS-TABLE (string->symbol func-str)))
               #| (displayln |# 
               #|   "FIXME: Normalize arguments so their expressions are resolved to a single variable or constant.") |#
               (define func-parameters-symbols func-argument-names)
               (define bindings (make-hash (zip cons 
                                                (function-inline-semantics-template-.parameters func-template)
                                                func-parameters-symbols)))
               (hash-set! bindings 
                          (syntax->datum #'function-name)
                          (syntax->datum #'function-return-variable))

               #| (displayln bindings) |#

               (define inlined-function (bind-parameters-to-arguments 
                                          func-str
                                          (function-inline-semantics-template-.body func-template)
                                          bindings))

               #| (displayln "inlined-function") |#
               #| (pretty-print (syntax->datum inlined-function)) |#

               ;; NOTE: Update pairs-list. It is processed in assignation macro
               (unless (equal? func-call-box-value 'func-call-box-not-set)
                   (begin
                     (define/contract template-value 
                                      function-inline-semantics-template/c
                                      (function-inline-semantics-template
                                        #'function-return-variable
                                        inlined-function
                                        (function-inline-semantics-template-.type func-template)
                                        (function-inline-semantics-template-.parameters func-template)
                                        (function-inline-semantics-template-.return-symbol func-template)))
                     (define updated-call-pairs-list
                       (append (read-state subfunc-call-box-value) ;; func-call list from children 
                               (read-state func-call-box-value) ;; func-call list from neighbors 
                               ;; on the same call depth level 
                               (list template-value)))
                     (write-state! func-call-box-value updated-call-pairs-list)))

               ;; NOTE: check if the function call is inside another function call
               ;; If it is so, return-vector literal should be genetated instead of a slice:
               ; arr[:] = f(g(h(1.0))) -> 
               ;                  h(h_ret, 1.0);
               ;                  g(g_ret, h_ret); 
               ;                  f(f_ret, g_ret);
               ;                  for (int slice_idx = 0; slice_idx < SLICE_LEN; slice_idx++)
               ;                     arr[slice-idx + SLICE_START] = f_ret[slice-idx + SLICE_START];
               (match (list "func return array:" func-slice-range 
                            "func is called inside argument:" (syntax-parameter-value #'func-is-called-inside-argument))
                 ((list "func return array:" #f
                        "func is called inside argument:" #f)
                  ;; FIXME: in C case ,#'func-return-symbol-stx should be transformed to (format "&~a" func-return-symbol)
                  (is-type_ func-output-base-type 
                            (if typecheck-mode-on 
                              (with-syntax-property 'get-value-name 
                                                    (syntax->datum #'function-return-variable)
                                                    #'typecheck-nonsense)
                              
                              (datum->syntax
                                   stx
                                   #| `(_var-ref ,#'func-return-symbol-stx) |#
                                   `(get-value ,#'function-return-variable)))))

                 ((list "func return array:" _ 
                        "func is called inside argument:" #f)
                  (is-type_ (landau-type func-output-base-type func-slice-range)
                            (if typecheck-mode-on
                              (with-syntax-property 'get-value-name 
                                                    (syntax->datum #'function-return-variable)
                                                    #'typecheck-nonsense)
                              
                              (with-syntax*
                                ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                 (func-return-symbol-stx (datum->syntax stx func-return-symbol)))
                                (datum->syntax
                                     stx
                                     #| `(_vector-ref ,#'func-return-symbol-stx (_var-ref ,#'slice-idx)) |#
                                     `(get-value ,#'function-return-variable "[" 0 ":" ,#'func-slice-range "]")
                                     )))))

                 ((list "func return array:" #f 
                        "func is called inside argument:" _)
                  ;; FIXME: in C case ,#'func-return-symbol-stx should be transformed to (format "&~a" func-return-symbol)
                  (is-type_ func-output-base-type
                            (if typecheck-mode-on
                              (with-syntax-property 'get-value-name 
                                                    (syntax->datum #'function-return-variable) 
                                                    #'typecheck-nonsense)
                              (datum->syntax
                                   stx
                                   `(get-value ,#'function-return-variable)))))

                 ((list "func return array:" _ 
                        "func is called inside argument:" _)
                  (is-type_ (landau-type func-output-base-type func-slice-range)
                            (if typecheck-mode-on
                              ;; NOTE: func calls list is populated in typecheck-mode. Function need to
                              ;  know it's arguments names, or in this case, function's return array name
                              (with-syntax-property 'get-value-name 
                                                    (syntax->datum #'function-return-variable)
                                                    #'typecheck-nonsense)
                              (with-syntax*
                                ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                 (func-return-symbol-stx (datum->syntax stx func-return-symbol)))
                                (datum->syntax
                                     stx
                                     `(get-value ,#'function-return-variable)))))))))))))))
           

(define-syntax (factor stx)
  (syntax-parse 
   stx
   [({~literal factor} primary "^" factor)
                 (error "exponention is not implemented")]
   [({~literal factor} number)
    (with-syntax ((expanded (local-expand-memo #'number 'expression '())))
      (syntax-track-origin (syntax/loc stx expanded) #'number #'expr))]))


(define-syntax (primary stx)
  (syntax-parse 
   stx
   [(_ (unop "-") expr-stx)
    (with-syntax*
      ((expanded (local-expand-memo #'expr-stx 'expression '()))
       (expanded-atom (atom-number #'expanded)))
      (let ((expanded-base-type (syntax-property #'expanded 'landau-type)))
        (cond
          [(or (equal? expanded-base-type 'dual-b)
               (is-slice-of-type 'dual-b expanded-base-type))
           (with-syntax*
             ((expanded-value (get-value-stx #'expanded))
              (expanded-derivative (get-derivative-stx #'expanded))
              (result-value (datum->syntax stx `(_rl-neg ,#'expanded-value)))
              (result-derivative (datum->syntax stx `(_rl-neg ,#'expanded-derivative))))
             (is-type_ expanded-base-type
                       (quasisyntax/loc stx
                                   (list #,#'result-value
                                         #,#'result-derivative))))]
          
          [(or (equal? expanded-base-type 'real)
               (is-slice-of-type 'real expanded-base-type))
           (is-type_ expanded-base-type
                     (if (syntax->datum #'expanded-atom)
                       #`#,(rl-neg (syntax->datum #'expanded-atom))
                       (datum->syntax stx
                                      `(_rl-neg ,#'expanded))))]
           ;; NOTE: there are no 'int slices
          [(equal? expanded-base-type 'int)
           (is-type_ expanded-base-type
                     (if (syntax->datum #'expanded-atom) 
                       #`#,(fx- 0 (syntax->datum #'expanded-atom))
                       (datum->syntax stx
                                      `(_int-neg ,#'expanded))))])))]

   [(_ (unop "+") expr-stx)
    (syntax/loc stx expr-stx)]

   [(_ number)
    (with-syntax ((expanded (local-expand-memo #'number 'expression '())))
      (syntax-track-origin (syntax/loc stx expanded) #'number #'expr))]
   ))

(define-syntax (element stx)
  (syntax-parse 
   stx
   [(_ number)
    (with-syntax ((expanded (local-expand-memo #'number 'expression '())))
      (syntax-track-origin (syntax/loc stx expanded) #'number #'expr))] 
   
   [(_ "(" expr-stx ")")
    (begin
     (syntax-track-origin (syntax/loc stx expr-stx) #'expr-stx #'expr))]))


(define-syntax (number stx)
  (syntax-parse stx
    (((~literal number) number-stx)
     (let* ((number-datum (syntax->datum #'number-stx))
            (num-type (if (exact? number-datum)
                        'int
                        'real)))
    (with-syntax
       ((finalized-num
          (match num-type
            ('int #'number-stx)
            ('real (datum->syntax stx (inexact->rl number-datum))))))
    (syntax-property
      #'finalized-num
      'landau-type
      num-type
      #t))))))

(begin-for-syntax
  (define-syntax-class type-spec
    #:attributes (t)
    (pattern
     ((~literal type) t:expr)))

  (define-syntax-class arg-spec
    #:attributes (name type)
    (pattern
     ({~literal argument} type:type-spec name:id))
     
    (pattern
     ({~literal other-argument} "," type:type-spec name:id)))
)


(begin-for-syntax
  (define (process-args! argnames argtypes func-name args add-arg!)
    (for/list ((argname (in-list (syntax-e argnames)))
               (argtype (in-list (syntax->datum argtypes))))
      (check-duplicate-argument-name args func-name
                                     (syntax->datum argname) argname)
      (add-arg! argname argtype args))))


;; NOTE: mutate current-arguments
(define/contract-for-syntax
    (define-func! stx ctx argnames argtypes body)
    (-> (syntax/c any/c)
        func-context/c
        (syntax/c any/c)
        (syntax/c any/c)
        (syntax/c any/c)

        (syntax/c any/c))
    (let* ((func-name (func-context-.function-name ctx))
           (func-return-value (func-context-.function-return-value ctx))
           (func-return-type (func-context-.function-return-type ctx))
           (args (func-context-.current-arguments ctx))
           (self-sufficient? (func-context-.self-sufficent-function? ctx)))
      (match (target-lang TARGET)
        ('racket 
         (let* ((add-arg!
                 (lambda (argname argtype args)
                   (datum->syntax argname (add-argument! args (syntax->datum argname)
                                                         (parse-type argtype)))))
                (args-list (process-args! argnames argtypes func-name args add-arg!)))
           (with-syntax ((ret (datum->syntax stx func-return-value))
                         (name (datum->syntax stx func-name))
                         (instantiated-func-var (instantiate (datum->syntax stx func-return-type))))
             (quasisyntax/loc stx
               (begin
                 (define (name #,@args-list)
                   (syntax-parameterize
                       ((func-context '#,ctx))
                     ;; Note: if func need derivatives then der array is allocated in func-body macro
                     (let ((ret #,#'instantiated-func-var))
                       #,body
                       ret)))
                 (provide name))))))
        ('ansi-c
         (let* ((name-str (symbol->string func-name))
                (func-return-value-str (symbol->string func-return-value))
                (add-arg! (lambda (argname argtype args)
                            (let* ((parsed-type (parse-type argtype))
                                   (arg-symbol (symbol->string (add-argument! args (syntax->datum argname)
                                                                              parsed-type))))
                                  (to-c-func-param parsed-type arg-symbol TARGET))))
                (args-list (process-args! argnames argtypes func-name args add-arg!))
                (arg-decl (c-func-arg-decl argnames args-list))
                (ctx-clone (func-context-copy ctx)))
           (with-syntax* 
             ((ret-str (datum->syntax stx func-return-value-str))
              (name-str-stx name-str)
              (func-ret (to-c-func-param func-return-type (syntax->datum #'ret-str) TARGET))
              (arg-decl arg-decl))

             (quasisyntax/loc
               stx
               (syntax-parameterize ((func-context '#,ctx))
                                    (c-func-decl "int" name-str-stx func-ret arg-decl (thunk #,body))))
             ))))))


;; TODO: this is a redundant function. It should be removed.
(define-for-syntax 
  (is-self-sufficent-function? name)
  (equal? name 'main))


(define-syntax (func stx)
  (syntax-parse stx
    (({~literal func} type:expr name:id "(" ({~literal arglist} arg*:arg-spec ...) ")"
                      "{" body "}")
     (let* ((name_ (syntax->datum #'name))
            (name-str (symbol->string name_))
            (func-return-value (gensym_ name-str))
            (func-return-type (parse-type (syntax->datum #'type)))
            (fake-src-pos 0)
            (func-data-list (hash-ref
                              functions-data-ht-GLOBAL
                              (var-symbol name-str fake-src-pos)))

            (dx-names-hash-set_ (cadr func-data-list))
            (der-table_ (car func-data-list))
            (real-vars-table_local (caddr func-data-list))
            (real-vars-table_with_global_func_names (rvt-new))
            (need-only-value-set_ (mutable-set))
            (need-derivatives-table_ (ndt-new))
            (funcs-return-vars-ht (for/hash ((k (hash-keys funcs-info-GLOBAL)))
                                    (values k (func-info-output-range (hash-ref funcs-info-GLOBAL k)))))
            (args (make-hash)))

       ;; NOTE: Union function variables keys and all functions return keys
       ;; FIXME: keys collision should be resolved (prohibited). Collision allowed 
       ;; only for the function name inside and outside the function
       (hash-union!
         real-vars-table_with_global_func_names
         real-vars-table_local
         funcs-return-vars-ht
         #:combine/key (lambda (k v1 v2) (if (equal? v1 v2)
                                           v1
                                           (raise-syntax-error #f (format "variable names collision: ~a" k) stx))))
       ;; NOTE: self-sufficent functions have der-annot and der-apply inside the body. 
       ;; They can not be called with dual-b args.
       ;; Not self-sufficent functions have not this terms, and may be called with dual-b args.
       ;; Caller should pass derivatives args like derivatives storage and mappings to such function. 
       ;; So a function declaration with derivatives arguments should be genetated. 

       (define self-sufficent-function? (is-self-sufficent-function? name_))
       (define current-variables (new-variables-nesting #f))
       (define/contract ctx func-context/c
                        (func-context
                          current-variables
                          name_
                          dx-names-hash-set_
                          der-table_
                          real-vars-table_with_global_func_names
                          need-only-value-set_
                          need-derivatives-table_
                          func-return-value
                          func-return-type
                          args
                          self-sufficent-function?
                          (new-mappings-table)))
       (define arg-names (for/list ((arg (in-list (syntax-e #'(arg*.name ...)))))
                           (syntax->datum arg)))
       (hash-set! MODULE-FUNCS-TABLE name_ (function-inline-semantics-template #'name 
                                                                               #'body 
                                                                               (to-landau-type stx func-return-type) 
                                                                               arg-names
                                                                               func-return-value))
       (begin
         (define-func! stx ctx #'(arg*.name ...) #'(arg*.type ...) #'body))))))


;; NOTE: Traverse the keys of all needed variables names and allocate mappings vectors if variable need derivatives
;;       Mappings and dx-idx-mappings are allocated for function return value and function arguments
;;       if thee belong to need-derivatives-set set         
;;       (e.g. df-table is empty) df-table :: Hash (String-key, (idx-type (Either (Set Int) Bool))
;; WARNING: Mutation of need-derivatives-set need-only-value-set mappings-table 
(define-for-syntax (make-mappings-decl-list! dx-name-str grouped-keys-table stx)
  (let* ((ctx (syntax-parameter-value #'func-context))
         (lst (make-mappings-decl-list-helper!
              dx-name-str
              grouped-keys-table
              (func-context-.real-vars-table ctx)
              (func-context-.der-table ctx)
              (func-context-.need-only-value-set ctx)
              (func-context-.need-derivatives-table ctx)
              (func-context-.current-variables ctx)
              (func-context-.mappings-table ctx)
              (target-lang TARGET)
              stx)))
    lst))

(define-for-syntax (make-der-decl-list! args stx)
  (let ((ctx (syntax-parameter-value #'func-context)))
    (make-der-decl-list-helper!
      args
      (func-context-.dx-names-hash-set ctx)
      (func-context-.need-derivatives-table ctx)
      (func-context-.current-variables ctx)
      (func-context-.current-arguments ctx)
      (target-lang TARGET)
      stx)))

(begin-for-syntax
  (define (make-function-self-derivatives-declarations stx ctx)
    (let* ((fake-src-pos 0)
           (func-return-type (func-context-.function-return-type ctx))
           (func-name (func-context-.function-name ctx))
           (func-name-vs (var-symbol func-name fake-src-pos))
           (maybe-func-return-size (hash-ref (func-context-.real-vars-table ctx)
                                             func-name-vs))
           (func-return-size (if maybe-func-return-size maybe-func-return-size 1))
           (dx-name->dx-size-table (hash-ref (func-context-.need-derivatives-table ctx)
                                             func-name-vs)))
      (for/list ((dx-name-str (hash-keys dx-name->dx-size-table)))
        (let* ((dx-mapped-size (mapping-sizes-.mapping-period (hash-ref dx-name->dx-size-table dx-name-str)))
               (dfdx-type (make-landau-type 'dual-b (fx* (mapping-sizes-.mapping-period dx-mapped-size) func-return-type)))
               (dual-r-var (add-variable!
                            (func-context-.current-variables ctx)
                            (datum->syntax #f (make-var-der-name func-name-vs 
                                                                 dx-name-str))
                            dfdx-type)))
          (datum->syntax stx
                 `(_define-var ,dual-r-var ,dfdx-type)))))))

(define-syntax (func-body stx)
  (syntax-parse stx
    [(_ body ...)
     (let* ((ctx (syntax-parameter-value #'func-context))
            (args (func-context-.current-arguments ctx))
            (self-sufficent-function? (func-context-.self-sufficent-function? ctx))
            (func-name (func-context-.function-name ctx))
            (func-name-str (symbol->string func-name))
            (fake-src-pos 0)
            (func-name-vs (var-symbol func-name-str fake-src-pos))
            (func-return-type (func-context-.function-return-type ctx))
            (grouped-keys-table (group-by-names (hash-keys (func-context-.der-table ctx)))))
       (with-syntax*
           ;; WARNING: Mutation of need-derivatives-set and need-only-value-set
           ((mappings-decl-list
              ;; NOTE: loop through all dx-names of program
              ;; NOTE: fill need-derivarives-table
              (filter not-void? (flatten (for/list ((dx-name-str (in-hash-keys (func-context-.dx-names-hash-set ctx))))
                                           (make-mappings-decl-list! dx-name-str grouped-keys-table stx)))))
            ;; NOTE: Generate syntax for derivatives storage
            ;; NOTE: Use only needed dx-names here
            (der-decl-list (make-der-decl-list! args stx))
            (i (unless (var-symbol/c func-name-vs)
                 (error "1475")))
            (func-return-type-final
              (cond
                ((set-member? (func-context-.need-only-value-set ctx) func-name-vs)
                 func-return-type)
                ((ndt-member? (func-context-.need-derivatives-table ctx) 
                              func-name-vs)
                 (make-dual-type func-return-type 'dual-l))
                ;; TODO: func value is not needed; do not need to run body
                (else func-return-type)))
            (func-basic-type-final (car (syntax->datum #'func-return-type-final)))
            (maybe-func-derivative-decl-list
              (cond
                ;; NOTE: If func is dual-l then allocate derivative vector and add der variable to variables
                ;; NOTE: Mappings are allocated in mappings-decl-list
                ((equal? (syntax->datum #'func-basic-type-final) 'dual-l)
                 (filter not-void? (make-function-self-derivatives-declarations stx ctx)))

                (else #'("")))))
         (set-need-only-value-set!
          grouped-keys-table 
          (func-context-.der-table ctx)
          (func-context-.need-only-value-set ctx))
         (set-func-context-.function-return-type! ctx (syntax->datum #'func-return-type-final))
        
         (match (target-lang TARGET)
           ('racket 
            (quasisyntax/loc stx
              (begin
                ;; TODO: make tests
                (begin #,@#'der-decl-list)
                (begin #,@#'mappings-decl-list)
                ;; TODO: make tests
                (begin #,@#'maybe-func-derivative-decl-list)
                (begin body ...))))
           ('ansi-c
            (quasisyntax/loc stx
              (string-append
               (string-append
                ;; TODO: make tests
                (apply string-append (filter not-void? (flatten (list #,@#'der-decl-list))))
                (apply string-append (flatten (list #,@#'mappings-decl-list)))
                ;; TODO: make tests
                (apply string-append (flatten (list #,@#'maybe-func-derivative-decl-list))))
               (string-append #,@#'(body ...))
               (c-return "0"))))
           )))]))


(define-for-syntax (need-derivatives? name-str)
  (hash-has-key? (func-context-.der-table (syntax-parameter-value #'func-context)) name-str))


(define/contract-for-syntax
  (add-dual-r-var! ctx var-name var-type)
  (-> func-context/c (syntax/c symbol?) type/c
      (values
        type/c 
        (hash/c string? variable/c)
        boolean?))

  (let* ((var-name-str (syntax->string var-name))
         (name-vs (var-symbol var-name-str 
                              (if (syntax-position var-name) 
                                (syntax-position var-name)
                                (error (format "bug: add-dual-r-var!: syntax-position ~a retured #f" var-name)))))
         (dx-mapped-sizes-table
           (hash-ref (func-context-.need-derivatives-table ctx) 
                     name-vs))
         (real-var-size (hash-ref (func-context-.real-vars-table ctx)
                                  name-vs))
         (var-size (if real-var-size real-var-size 1))
         (dual-r-vars-table
           (for/hash ((dx-name-str (in-list (hash-keys dx-mapped-sizes-table))))
             (let* ((dx-mapped-size (hash-ref dx-mapped-sizes-table dx-name-str))
                    (dfdx-type (make-landau-type 'dual-r 
                                                 (fx* (mapping-sizes-.mapping-period dx-mapped-size) var-size)))
                    (dfdx-var-name (make-var-der-name name-vs dx-name-str))
                    (curr-vars (func-context-.current-variables ctx))
                    (dual-r-sym (add-variable! curr-vars
                                               (datum->syntax #f dfdx-var-name)
                                               dfdx-type))
                    (dfdx-variable (search-variable dfdx-var-name curr-vars)))
               (values dx-name-str dfdx-variable)))))
    (values
      (make-dual-type var-type 'dual-l)
      dual-r-vars-table
      #t)))

(define-syntax (decl-block stx)
  (syntax-parse stx
    [(_ body ...)
     (datum->syntax stx
       `(_begin
         ,@#'(body ...)))]))

(define-syntax (var-decl stx)
  (syntax-parse stx
    (({~literal var-decl};; if expr is array get-value, then emit declaration and assignation syntax objects  
                ((~literal type) ((~literal array-type) basic-type "[" num "]")) name:id (~seq "=" value:expr)) 
     (datum->syntax stx `(decl-block 
                           (var-decl (type (array-type ,#'basic-type "[" ,#'num "]")) ,#'name)
                           (assignation ,#'name "[" ":" "]" "=" ,#'value))))

    (({~literal var-decl}
      type:expr name:id (~optional (~seq "=" value:expr)
                                   #:defaults ((value #'notset))))
     (let* ((name_ (syntax->datum #'name))
            (type (to-landau-type stx (parse-type (syntax->datum #'type))))
            (name-str (symbol->string name_))
            (src-pos (syntax-position #'name))
            (name-vs (var-symbol name-str (if src-pos
                                            src-pos 
                                            (error (format "bug: var-decl: syntax-position ~a retured #f" #'name)))))
            (ctx (syntax-parameter-value #'func-context))
            (base-type (car type))
            (assign-if-set (datum->syntax 
                            stx 
                            (if (equal? (syntax->datum #'value) 'notset)
                              `(_nothing)
                              `(assignation ,#'name "=" ,#'value)))))
       ;; NOTE: checks have been done in backrun 
       #| (check-duplicate-variable-name name_ #'name) |#
       (unless (var-symbol/c name-vs)
         (error (format "bug: var-decl: name-vs is not var-symbol/c: ~a" name-vs)))
       (match base-type
         ('real
          (let*-values
              ;; NOTE: dx-mapped-sizes :: Hash name mapped-size
              ;; NOTE: dual-r-vars :: Hash name dx-gen-symbol
              (((final-type dual-r-vars need-variable)
                (cond
                  ((set-member? (func-context-.need-only-value-set ctx) name-vs)
                   (values type #f #t))
                  ((ndt-member? (func-context-.need-derivatives-table ctx) name-vs)
                   (add-dual-r-var! ctx #'name type))
                  ;; NOTE: Variable value and derivative is not used
                  (else (values type #f #f))))
                       
               ((name-sym) (add-variable!
                       (func-context-.current-variables ctx) #'name final-type))
               ((name-sym-stx) (datum->syntax stx name-sym)))

          #| (displayln "ndt:") |#
          #| (pretty-print (func-context-.need-derivatives-table ctx)) |#
           #| (displayln (format "name-vs final-type dual-r-vars need-variable ~a ~a ~a ~a" |# 
           #|                    name-vs final-type dual-r-vars need-variable)) |# 
            (cond
              ((and need-variable (not dual-r-vars))
               (datum->syntax stx
                              `(_begin
                                (_define-var ,name-sym ,type)
                               ,assign-if-set)))
              
              ((and need-variable (not (hash-empty? dual-r-vars)))
               (with-syntax*
                 ((der-decl-list
                   ;; NOTE: declaration for eash dx-name (which is needed for current name)
                   (for/list ((dx-name-str (in-list (hash-keys dual-r-vars))))
                             (let* ((dual-r-var (hash-ref dual-r-vars dx-name-str))
                                    (dual-r-var-type (variable-type dual-r-var))
                                    (dual-r-var-symb (variable-symbol dual-r-var)))
                                   (datum->syntax stx
                                                  `(_define-var ,dual-r-var-symb ,dual-r-var-type))))))
                ;  (quasisyntax/loc stx
                ;                   (#,mbegin
                ;                    #,(define-var stx name-sym type)
                ;                    (#,mbegin #,@#'der-decl-list)
                ;                    #,assign-if-set))
                 (datum->syntax stx
                                `(_begin
                                  (_define-var ,name-sym ,type)
                                  (_begin ,@#'der-decl-list)
                                  ,assign-if-set))))
              ;; NOTE: skip declaration if value is not used
              (else
               ;; FIXME: keep it For testing
                 (datum->syntax stx
                                `(_begin
                                  (_define-var ,name-sym ,type)
                                 ,assign-if-set))))))

         ('int
          (cond
            ((is-array-type? type) (raise-syntax-error #f "'int arrays are not supported yet" stx))
            (else
             (let* ((name-sym (add-variable!
                               (func-context-.current-variables ctx) #'name type))
                    (name-sym-stx (datum->syntax stx name-sym))) 
               (datum->syntax stx
                              `(_begin
                                (_define-var ,name-sym ,type)
                               ,assign-if-set))))))
         (else (raise-syntax-error #f "unsupported type. Expect: real[INT] STRING | real STRING | int STRING" stx)))))))   


(define-for-syntax (check-duplicate-constant name stx)
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "duplicate constant declaration" stx)))


(define-syntax (constant stx)
  (syntax-parse stx
    (({~literal constant} "const" type name:id "=" (~or value:expr (~seq "{" ({~literal parlist} arr-item*:par-spec ...) "}")))
     (let* ((type (parse-type (syntax->datum #'type)))
            (base-type (car type)))
           
           (check-duplicate-constant (syntax->datum #'name) stx)
           (cond
             ((is-array-type? type)
              (let*
               ((expanded-list
                 (for/list ((arr-item (in-list (syntax-e #'(arr-item*.par-value ...)))))
                           (local-expand-memo arr-item 'expression '())))
                (constant-array-symbol (gensym_ (syntax->datum #'name))))
               (with-syntax ((constant-array-stx (datum->syntax stx constant-array-symbol)))
                 (hash-set! constants (syntax->datum #'name)
                            (constant constant-array-symbol type (list->rl-vector expanded-list)))
                 (quasisyntax/loc stx
                                  (begin
                                   #,(datum->syntax
                                   stx
                                   `(_define-var ,#'constant-array-stx ,type (_rl-vector ,expanded-list) #t)))))))
             (else
              (let* ((expanded (local-expand-memo #'value 'expression '()))
                     (value (atom-number expanded))
                     (value-type (syntax-property expanded 'landau-type))
                     (value
                      (match (list base-type "<-" value-type)
                        ((list 'int "<-" 'int)
                         value)
                        
                        ((list 'int "<-" 'real)
                         (raise-syntax-error #f "real value is not allowed for integer constant" stx))
                        
                        ((list 'real "<-" 'int)
                         (->rl value))
                        
                        ((list 'real "<-" 'real)
                         (inexact->rl value)))))
                    
                    (hash-set! constants (syntax->datum #'name)
                               (constant value type #f))
                    (datum->syntax stx '(_nothing)))))))))


(define-for-syntax (check-duplicate-parameter name stx)
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "duplicate parameter declaration" stx)))


(define-syntax (parameter stx)
  (syntax-parse stx
    (({~literal parameter} "parameter" "[" size "]" name:id)
     
     (check-duplicate-parameter (syntax->datum #'name) stx)
     (let* ((expanded-size (local-expand-memo #'size 'expression '())))
              
       (hash-set! parameters (syntax->datum #'name)
                  (list expanded-size))
       (datum->syntax stx (_empty-statement))))))


(define-syntax (get-value stx)
  (syntax-parse stx
    (get-value:get-value-cls

     (with-syntax ((name #'get-value.name))
       #;(displayln (string->symbol (format "~a~a" (build-string padding (lambda (x) #\space))
                                   'get-value)))
       (define tik (current-inexact-milliseconds))
       (let*-values
           (((name-symb) (syntax->datum #'name))
            ;; NOTE: Check enabled modes
            ;; NOTE: if typecheck-mode-on is #t then syntax is not genetated. Used for type info propagation.
            ((typecheck-mode-on) (syntax-parameter-value #'typecheck-mode)) 
            ((expand-only-value-on) (syntax-parameter-value #'expand-value-only))
            ((get-name-mode-on) (syntax-parameter-value #'get-name-mode))

            ;; NOTE: typechecked in backrun
            ((slice-range) 'range-placeholder)
            ((name-str) (symbol->string name-symb))
            ((dx-name-str-in-current-al) (syntax-parameter-value #'dx-name-in-current-al))
            ((ctx) (syntax-parameter-value #'func-context))
            ((idx) (timeit! TIME-TABLE 'local-expand-memo (thunk (local-expand-memo #'get-value.index 'expression '()))))
            ((value type name-is-dx src-pos)
             (timeit/values! TIME-TABLE 'search-name (thunk (search-name stx ctx name-symb #'name dx-name-str-in-current-al get-name-mode-on idx))))

            ((getter-info) (make-getter-info (if #'get-value.index #'get-value.index #'#f)
                                             (if #'get-value.slice-colon #'get-value.slice-colon #'#f)
                                             (timeit! TIME-TABLE 'to-landau-type (thunk (to-landau-type stx type)))))
            ((is-slice is-cell is-var) (values (getter-is-slice? getter-info)
                                               (getter-is-cell? getter-info)
                                               (getter-is-var? getter-info)))
            ((index-start-expanded) (if is-slice 
                                      (if (syntax->datum #'get-value.index-start)
                                        (timeit! TIME-TABLE 'local-expand-memo-2 (thunk (local-expand-memo #'get-value.index-start 'expression '())))
                                        0)
                                      #f))

            ((name-vs) (var-symbol name-str src-pos))
            ; ((value (resolve-constant-arrays ctx value-unresolved-constant-array #'get-value.index)))
            ((declared-type) (car type))
            ((type) (if (and expand-only-value-on (equal? declared-type 'dual-l))
                      'real
                      declared-type))
            ((base-type) type))
         ; (println (format "debug: get-value: dx-name-str-in-current-al: ~a" dx-name-str-in-current-al))
         (define get-value-head-time (current-inexact-milliseconds))
         (with-syntax* ((value value)
                        (index-start-expanded index-start-expanded)
                        ;; NOTE: nonsense values used to retrieve info in compile time. They should not be in genetated code.
                        ;; If they are, this is a bug.
                        ;; FIXME format -> error
                        (real-nonsense #'(format "bug: real-nonsense in runtime"))
                        (nonsense #'(format "bug: nonsense in runtime"))
                        (type-check-nonsense-val #'(format "bug: type-check-nonsense-val in runtime"))
                        (type-check-nonsense-der #'(format "bug: type-check-nonsense-der in runtime"))
                        (dual-bundle-nonsense
                         (syntax/loc stx (list type-check-nonsense-val type-check-nonsense-der))))
           #| (displayln name-vs) |#
           #| (displayln (format "(list (getter-info-type getter-info) base-type): ~a" |# 
           #|                    (list (getter-info-type getter-info) base-type))) |# 
           #| (displayln (syntax-parameter-value #'func-call-box)) |#
           (with-syntax ((r (with-syntax-property 'get-value-name name-symb
             (with-syntax-property 'getter-info getter-info
               (if get-name-mode-on
                 (datum->syntax stx '(_nothing))
                 (match (list (getter-info-type getter-info) base-type)

                   ((list 'var 'dual-l)
                    (cond
                      (name-is-dx
                        (is-type_ 'dual-b
                                  (if typecheck-mode-on
                                    #'dual-bundle-nonsense
                                    (with-syntax*
                                      ((derivative
                                         (datum->syntax stx '_1.0)))
                                      (quasisyntax/loc stx
                                                       (list
                                                         #,(datum->syntax stx `(_var-ref ,#'value)) 
                                                         #,#'derivative))))))
                      (else          
                        (is-type_
                          'dual-b
                          (if typecheck-mode-on
                            #'dual-bundle-nonsense
                            (with-syntax*
                              ((der-var-synt
                                 (datum->syntax stx (get-der-variable-symbol name-vs dx-name-str-in-current-al stx #f)))
                               (dx-idxs-mappings
                                 (datum->syntax stx (get-dx-idxs-mappings-variable-symbol name-vs dx-name-str-in-current-al stx #f)))
                               (inv-mapping-period (get-inv-mapping-period name-vs dx-name-str-in-current-al stx))
                               (al-index (datum->syntax stx 'al_index_name_symbol))
                               ;; NOTE: genetated function which return derivative value
                               ; using mappings and inverse mapping to get index where derivative is stored
                               (func-name-stx (datum->syntax stx 'get_dfdx_var))
                               ;; FIXME get-der-variable-dx-range is unavailable for dx dx' 
                               (dx-range (get-der-variable-dx-range name-vs dx-name-str-in-current-al stx))
                               #| (dx-range (check-result |# 
                               #|             stx |# 
                               #|             (format "bug2: get-der-variable-dx-range retured #f for ~a ~a" |#
                               #|                     (var-symbol-.name name-vs) dx-name-str-in-current-al) |# 
                               #|             (get-der-variable-dx-range name-vs dx-name-str-in-current-al stx))) |#
                               (dx-parameter-size (get-dx-parameter-size stx parameters dx-name-str-in-current-al))
                               (all-derivatives-are-used? (check-if-all-derivatives-are-used
                                                            ctx
                                                            (syntax->datum #'dx-parameter-size) 
                                                            (syntax->datum #'dx-range)
                                                            (syntax->datum #'dx-idxs-mappings)))
                               (derivative
                                 (match (list 'dont-have-mappings (or (equal? #f (syntax->datum #'dx-idxs-mappings))
                                                                      (equal? #f (syntax->datum #'inv-mapping-period)))
                                              'all-derivatives-are-used (syntax->datum #'all-derivatives-are-used?))
                                   ((list 'dont-have-mappings #true
                                          'all-derivatives-are-used _)
                                    (datum->syntax stx '_0.0))
                                   ((list 'dont-have-mappings #false
                                          'all-derivatives-are-used #false)
                                    (with-syntax
                                      ((args-list (list
                                                    #'al-index
                                                    (datum->syntax stx `(_var-ref ,#'inv-mapping-period))
                                                    (datum->syntax stx `(_var-ref ,#'dx-idxs-mappings))
                                                    (datum->syntax stx `(_var-ref ,#'der-var-synt)))))
                                      (datum->syntax
                                        stx
                                        `(_pure-func-call ,#'func-name-stx ,#'args-list))))
                                   ((list 'dont-have-mappings #false
                                          'all-derivatives-are-used #true)
                                    (datum->syntax
                                        stx
                                        `(_vector-ref ,#'der-var-synt ,#'al-index))))))

                              (quasisyntax/loc stx
                                               (list
                                                 #,(datum->syntax stx `(_var-ref ,#'value)) 
                                                 #,#'derivative))))))))

                   ((list 'array 'dual-l)
                    ;; NOTE: the only case of array varables usage is 
                    ; passing them to a function. Dual array varables passed to
                    ; a function are not really passed, because dual funciton  
                    ; is inlined, so it really does not matter what to emit here.
                    (is-type_ 'dual-b
                                  #'dual-bundle-nonsense))

                   ((list (or 'cell 'slice) 'dual-l)
                    (is-type_ 
                      (if is-slice (landau-type 'dual-b slice-range) 'dual-b)
                      (if typecheck-mode-on
                        #'dual-bundle-nonsense
                        (with-syntax*
                          ((der-vec-synt
                             (datum->syntax stx (get-der-variable-symbol name-vs dx-name-str-in-current-al stx #f)))
                           (dx-idxs-mappings
                             (datum->syntax stx (get-dx-idxs-mappings-variable-symbol name-vs dx-name-str-in-current-al stx #f)))
                           (dx-mapped-size (check-result
                                             stx
                                             (format "bug1: get-der-variable-dx-range returned #f for ~a ~a" #'name dx-name-str-in-current-al)
                                             (get-der-variable-dx-range name-vs 
                                                                        dx-name-str-in-current-al stx)))
                           (inv-mapping-period (get-inv-mapping-period name-vs dx-name-str-in-current-al stx))
                           (name-string (symbol->string name-symb))

                           (dx-parameter-size (get-dx-parameter-size stx parameters dx-name-str-in-current-al))
                           (dx-name-str-in-current-al dx-name-str-in-current-al)
                           (al-index (datum->syntax stx 'al_index_name_symbol))
                           (slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                           (full-idx (if is-slice 
                                       (datum->syntax stx `(_int+ ,#'index-start-expanded (_var-ref ,#'slice-idx)))
                                       #'get-value.index))
                           (dual-b-value (datum->syntax stx `(_vector-ref ,#'value ,#'full-idx)))
                           (all-derivatives-are-used? (check-if-all-derivatives-are-used
                                                        ctx
                                                        (syntax->datum #'dx-parameter-size) 
                                                        (syntax->datum #'dx-mapped-size)
                                                        (syntax->datum #'dx-idxs-mappings)))
                           (dual-b-derivative
                             (cond
                               (name-is-dx
                                 (cond 
                                   ;; NOTE: dx vector is always dual-l, because even if derivatives are not set,
                                   ;; it can be used as r-val when l-val is dual-l and dx[i] ' dx[k] must be I_ik
                                   ((or (equal? #f (syntax->datum #'dx-idxs-mappings))
                                        (equal? #f (syntax->datum #'inv-mapping-period)))
                                    (datum->syntax stx
                                                   `(_if (_int= (_var-ref ,#'al-index) ,#'full-idx) _1.0 _0.0)))
                                   ((syntax->datum #'all-derivatives-are-used?)
                                    (datum->syntax
                                      stx
                                      `(_vector-ref ,#'der-vec-synt (_int+ (_int* ,#'full-idx ,#'dx-mapped-size)
                                                                           (_var-ref ,#'al-index)))))
                                   (else
                                     (with-syntax
                                       ((func-name-stx (datum->syntax stx 'get_dfdx_cell_dx))
                                        (args-list (list
                                                     #'full-idx
                                                     #'dx-mapped-size
                                                     #'al-index
                                                     (datum->syntax stx `(_var-ref ,#'inv-mapping-period))
                                                     (datum->syntax stx `(_var-ref ,#'dx-idxs-mappings))
                                                     (datum->syntax stx `(_var-ref ,#'der-vec-synt)))))
                                       (datum->syntax
                                         stx
                                         `(_pure-func-call ,#'func-name-stx ,#'args-list)))
                                     )))
                               (else
                                 (cond
                                   ((or (equal? #f (syntax->datum #'dx-idxs-mappings))
                                        (equal? #f (syntax->datum #'inv-mapping-period)))
                                    ;; NOTE: Dual-b can has no mappings for current dx, but in has for another
                                    (datum->syntax stx '_0.0))
                                    ;; TODO: move this binding to the assignation term, then use it here from syntax-parameter
                                    ((syntax->datum #'all-derivatives-are-used?)
                                      (datum->syntax
                                         stx
                                         `(_vector-ref ,#'der-vec-synt (_int+ (_int* ,#'full-idx ,#'dx-mapped-size)
                                                                              (_var-ref ,#'al-index)))))
                                    (else
                                      (with-syntax
                                        ((func-name-stx (datum->syntax stx 'get_dfdx_cell))
                                         (args-list (list
                                                      #'full-idx
                                                      #'dx-mapped-size
                                                      #'al-index
                                                      (datum->syntax stx `(_var-ref ,#'inv-mapping-period))
                                                      (datum->syntax stx `(_var-ref ,#'dx-idxs-mappings))
                                                      (datum->syntax stx `(_var-ref ,#'der-vec-synt)))))
                                        (datum->syntax
                                          stx
                                          `(_pure-func-call ,#'func-name-stx ,#'args-list)))))))))
                          (syntax/loc stx
                                      (list 
                                        dual-b-value
                                        dual-b-derivative))))))

                   ((list 'var 'real)
                    (is-type_ type
                              (if typecheck-mode-on
                                #'nonsense
                                (if (atom-number #'value)
                                  #'value
                                  (datum->syntax stx `(_var-ref ,#'value))))))

                   ((list 'array 'real) 
                    ;; FIXME not shoure if I need to dereference in C backend
                    ;; NOTE: the only case of array varables usage is 
                    ; passing them to a function.
                    (is-type_ type
                              (if typecheck-mode-on
                                #'nonsense
                                (if (atom-number #'value)
                                  #'value
                                  (datum->syntax stx `(_var-ref ,#'value))))))

                   ((list (or 'cell 'slice) 'real)
                    ;; FIXME: if resolved constant array cell, just insert it's value
                    (cond
                      ((atom-number #'value)
                       (is-type_ 
                         (if is-slice (landau-type 'real slice-range) 'real) 
                         #'value))
                      (else 
                        (with-syntax* ((value-stx (datum->syntax stx #'value))
                                       (slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                       (full-idx (if is-slice
                                                   (datum->syntax stx `(_int+ ,#'index-start-expanded (_var-ref ,#'slice-idx)))
                                                   #'get-value.index)))
                                      (is-type_
                                        (if is-slice (landau-type 'real slice-range) 'real)
                                        (if typecheck-mode-on
                                          #'real-nonsense
                                          (datum->syntax
                                            stx
                                            `(_vector-ref ,#'value-stx ,#'full-idx))))))))

                   ((list 'var 'int)
                    (is-type_ type
                              ;; NOTE: ignore typecheck mode for indexes
                              (if (atom-number #'value)
                                #'value
                                (datum->syntax stx `(_var-ref ,#'value)))))

                   ((list 'cell 'int)
                    (raise-syntax-error #f "int arrays are not supported, yet" #'name))

                   ((list 'slice 'int)
                    (raise-syntax-error #f "int slices are not supported, yet" #'name)))
                 )))))
             #| (define tok (current-inexact-milliseconds)) |#
             #| (hash-update! TIME-TABLE 'head-time (lambda (old-time) (fl+ old-time (fl- get-value-head-time tik)))) |#
             #| (hash-update! TIME-TABLE 'tail-time (lambda (old-time) (fl+ old-time (fl- tok get-value-head-time)))) |#
             #| (displayln (format "get-value head-time: ~a" (hash-ref TIME-TABLE 'head-time 0.0))) |#
             #| (displayln (format "  get-value search-name time: ~a" (hash-ref TIME-TABLE 'search-name 0.0))) |#
             #| #1| (displayln (format "    get-value search-constant time: ~a" (hash-ref TIME-TABLE 'search-constant 0.0))) |1# |#
             #| (displayln (format "    get-value search-variable time: ~a" (hash-ref TIME-TABLE 'search-variable 0.0))) |#
             #| (displayln (format "    get-value HASH-HIT: ~a" HASH-HIT)) |#
             #| #1| (displayln (format "    get-value search-argument time: ~a" (hash-ref TIME-TABLE 'search-argument 0.0))) |1# |#
             #| (displayln (format "  get-value local-expand-memo time: ~a" (hash-ref TIME-TABLE 'local-expand-memo 0.0))) |#
             #| (displayln (format "  get-value local-expand-memo-2 time: ~a" (hash-ref TIME-TABLE 'local-expand-memo-2 0.0))) |#
             #| (displayln (format "  get-value to-landau-type time: ~a" (hash-ref TIME-TABLE 'to-landau-type 0.0))) |#
             #| (displayln (format "get-value tail-time: ~a" (hash-ref TIME-TABLE 'tail-time 0.0))) |#
             #'r
             )))))))

(define/contract-for-syntax
  (make-set-all-derivatives-to-zero-list stx ctx name-vs al-index-symb)
  (-> (syntax/c any/c) func-context/c var-symbol/c (syntax/c any/c)
      (listof (syntax/c any/c)))
  (let ((dx-names (hash-keys (hash-ref (func-context-.need-derivatives-table ctx)
                                       name-vs))))
    (for/list ((dx-name-str (in-list dx-names)))
      (let ((maybe-mappings (get-mappings-variable-symbol name-vs dx-name-str stx #f)))
        (if maybe-mappings
          (with-syntax*
            ((al-index-symb al-index-symb)
             (der-vec-synt (datum->syntax stx (get-der-variable-symbol name-vs dx-name-str stx)))
             (dfdx-mapping-sizes (hash-ref (func-context-.need-derivatives-table ctx)
                                           name-vs))
             (mapping-vec-len (mapping-sizes-.mapping-period (hash-ref (syntax->datum #'dfdx-mapping-sizes) dx-name-str)))
             (mappings-synt (datum->syntax stx maybe-mappings))
             (dx-name-str_ dx-name-str)
             (int-type (make-landau-type 'int #f))
             (loop-var (datum->syntax stx (gensym_ 'loop_var)))
             (mapped-idx (datum->syntax stx 'mapped_idx)))
            (datum->syntax
              stx
              `(_for ,#'loop-var 0 ,#'mapping-vec-len
                     (_let-int ,#'mapped-idx (_var-ref ,#'loop-var)
                               (_let-int ,#'al-index-symb (_int-vector-ref ,#'mappings-synt (_var-ref ,#'loop-var))
                                         (_vector-set! ,#'der-vec-synt (_var-ref ,#'mapped-idx) _0.0))))))
          ;; TODO: Do not use mapping if dual-l variable has have-mapping flag equal to #f
          (raise-syntax-error #f (format "bug: no mapping found for dual-l variable ~a" (var-symbol-.name name-vs)) stx)))))
  )

(define/contract-for-syntax
  (make-der-assignation-loops-list stx ctx name-vs dx-names assignation-rigth-part al-index-symb)
  (-> (syntax/c any/c) func-context/c var-symbol/c (listof string?) (syntax/c any/c) (syntax/c any/c)
      (listof (syntax/c any/c)))
  (for/list ((dx-name-str (in-list dx-names)))
    (let ((maybe-mappings (get-mappings-variable-symbol name-vs dx-name-str stx #f)))

      (if maybe-mappings
        (with-syntax*
          ((al-index-symb al-index-symb)
           (der-vec-synt (datum->syntax stx (get-der-variable-symbol name-vs dx-name-str stx)))
           (dfdx-mapping-sizes (hash-ref (func-context-.need-derivatives-table ctx)
                                         name-vs))
           (mapping-vec-len (mapping-sizes-.mapping-period (hash-ref (syntax->datum #'dfdx-mapping-sizes) dx-name-str)))
           (mappings-synt (datum->syntax stx maybe-mappings))
           (dx-name-str_ dx-name-str)
           (int-type (make-landau-type 'int #f))
           (loop-var (datum->syntax stx (gensym_ 'loop_var)))
           (mapped-idx (datum->syntax stx 'mapped_idx))
           (dual-b-derivative
             (get-derivative-stx 
               #|   ;; FIXME get-value emited by func-call will fail here |#
               #|   ; add function return variables var-decl before expansion |#
               (extract
                 (local-expand-memo
                   #`(syntax-parameterize
                       ((dx-name-in-current-al '#,dx-name-str))
                       #,assignation-rigth-part) 'expression '() #:reset-memo #t)))
             ))
          (datum->syntax
            stx
            `(_for ,#'loop-var 0 ,#'mapping-vec-len
                   (_begin
                     (_let-int ,#'mapped-idx (_var-ref ,#'loop-var)
                               (_let-int ,#'al-index-symb (_int-vector-ref ,#'mappings-synt (_var-ref ,#'loop-var))
                                         ;; NOTE: value includes references to der-vectors and al-index-symb
                                         (_vector-set! ,#'der-vec-synt (_var-ref ,#'mapped-idx) ,#'dual-b-derivative)))))))
        (raise-syntax-error #f (format "bug: no mapping found for dual-l variable ~a" (var-symbol-.name name-vs)) stx)))))

(define/contract-for-syntax
  (make-set-all-derivatives-to-zero/array
    stx
    ctx
    name-vs
    al-index-symb
    left-hand-getter-info 
    slice-range
    index-start-expanded_
    index-exp)
  (-> (syntax/c any/c)
      func-context/c
      var-symbol/c
      (syntax/c symbol?)
      getter-info/c
      (syntax/c any/c)
      (syntax/c any/c)
      (syntax/c any/c)
      (listof (syntax/c any/c)))
  (let ((dx-names (hash-keys (hash-ref (func-context-.need-derivatives-table ctx)
                                       name-vs))))
    (for/list ((dx-name-str (in-list dx-names)))
      (let ((df-range (hash-ref (func-context-.real-vars-table ctx) name-vs))
            (maybe-mappings (get-mappings-variable-symbol name-vs dx-name-str stx #f)))
        (if maybe-mappings
          (with-syntax*
            ((al-index-symb al-index-symb)
             (slice-range slice-range)
             (index-start-expanded_ index-start-expanded_)
             (der-vec-symb (get-der-variable-symbol name-vs dx-name-str stx))
             (der-vec-synt (datum->syntax stx #'der-vec-symb))
             ;; FIXME get-der-variable-dx-range is unavailable for dx dx' 
             (dx-range (get-der-variable-dx-range name-vs dx-name-str stx))
             #| (dx-range (check-result |# 
             #|             stx |# 
             #|             (format "bug2: get-der-variable-dx-range retured #f for ~a ~a" (var-symbol-.name name-vs) dx-name-str) |# 
             #|             (get-der-variable-dx-range name-vs dx-name-str stx))) |#
             (dx-parameter-size (get-dx-parameter-size stx parameters dx-name-str))
             (all-derivatives-are-used? (check-if-all-derivatives-are-used ctx
                                                                           (syntax->datum #'dx-parameter-size) 
                                                                           (syntax->datum #'dx-range)
                                                                           (syntax->datum maybe-mappings)))
             (mappings-vec-len (fx* (syntax->datum #'dx-range) df-range))
             (mappings-synt (datum->syntax stx maybe-mappings))
             (mappings-full-idx (datum->syntax stx 'mappings_full_idx_symbol))
             (mapped-idx (datum->syntax stx 'mapped_idx))
             (dx-name-str_ dx-name-str)
             (slice-idx (datum->syntax stx slice-idx-name-GLOBAL)))
            (match (list "getter-is-slice?" (getter-is-slice? left-hand-getter-info)
                       "all-derivatives-are-used?" (syntax->datum #'all-derivatives-are-used?))
            ((list "getter-is-slice?" #true
                   "all-derivatives-are-used?" #false)
             (assign-derivative-to-slice stx
                                         #'slice-idx
                                         #'slice-range
                                         #'mapped-idx
                                         #'mappings-full-idx
                                         #'al-index-symb
                                         #'mappings-synt
                                         #'der-vec-synt
                                         #'dx-range
                                         #'index-start-expanded_
                                         #'mappings-vec-len
                                         _0.0))
            ((list "getter-is-slice?" #false
                   "all-derivatives-are-used?" #false)
             (assign-derivative-to-cell stx
                                        #'mapped-idx
                                        index-exp
                                        #'mappings-full-idx
                                        #'al-index-symb
                                        #'mappings-synt
                                        #'der-vec-synt
                                        #'dx-range
                                        #'mappings-vec-len
                                        _0.0))
            ((list "getter-is-slice?" #false
                   "all-derivatives-are-used?" _)
             (assign-derivative-to-cell-dense stx
                                              #'mapped-idx
                                              index-exp
                                              #'mappings-full-idx
                                              #'al-index-symb
                                              #'mappings-synt
                                              #'der-vec-synt
                                              #'dx-range
                                              _0.0))
            ((list "getter-is-slice?" #true
                   "all-derivatives-are-used?" _)
             (assign-derivative-to-slice-dense stx
                                               #'slice-idx
                                               #'slice-range
                                               #'mapped-idx
                                               #'mappings-full-idx
                                               #'al-index-symb
                                               #'mappings-synt
                                               #'der-vec-synt
                                               #'dx-range
                                               #'index-start-expanded_
                                               #'mappings-vec-len
                                               _0.0))))
          (raise-syntax-error #f (format "bug: no mapping found for dual-l variable ~a" (var-symbol-.name name-vs)) stx))))))

(define/contract-for-syntax
  (make-der-assignation-loops-list/array
    stx
    ctx
    name-vs 
    dx-names 
    assignation-rigth-part 
    al-index-symb 
    left-hand-getter-info 
    slice-range
    index-start-expanded_
    index-exp)
  (-> (syntax/c any/c)
      func-context/c
      var-symbol/c
      (listof string?)
      (syntax/c any/c)
      (syntax/c symbol?)
      getter-info/c
      (syntax/c any/c)
      (syntax/c any/c)
      (syntax/c any/c)
      (listof (syntax/c any/c)))
  (for/list ((dx-name-str (in-list dx-names)))
    (let ((df-range (hash-ref (func-context-.real-vars-table ctx) name-vs))
          (maybe-mappings (get-mappings-variable-symbol name-vs dx-name-str stx #f)))
      (if maybe-mappings
        (with-syntax*
          ((al-index-symb al-index-symb)
           (slice-range slice-range)
           (index-start-expanded_ index-start-expanded_)
           (der-vec-symb (get-der-variable-symbol name-vs dx-name-str stx))
           (der-vec-synt (datum->syntax stx #'der-vec-symb))
           (dx-range (check-result 
                       stx 
                       (format "bug: get-der-variable-dx-range retured #f for ~a ~a" 
                               (var-symbol-.name name-vs) 
                               dx-name-str) 
                       (get-der-variable-dx-range name-vs dx-name-str stx)))
           (mappings-vec-len (fx* (syntax->datum #'dx-range) df-range))
           (dx-parameter-size (get-dx-parameter-size stx parameters dx-name-str))
           (all-derivatives-are-used? (check-if-all-derivatives-are-used ctx
                                                                         (syntax->datum #'dx-parameter-size) 
                                                                         (syntax->datum #'dx-range)
                                                                         (syntax->datum maybe-mappings)))
           (mappings-synt (datum->syntax stx maybe-mappings))
           (mappings-full-idx (datum->syntax stx 'mappings_full_idx_symbol))
           (mapped-idx (datum->syntax stx 'mapped_idx))
           (dx-name-str_ dx-name-str)
           (slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
           ;;(dummy1 #RRassignation-rigth-part)
           (dual-b-derivative
             (timeit! TIME-TABLE 'get-derivative-stx
                      (thunk
                        (get-derivative-stx
                          (extract
                            (local-expand
                              #`(syntax-parameterize
                                  ((dx-name-in-current-al '#,dx-name-str))
                                  #,assignation-rigth-part) 'expression '())))
                        ))))
          (match (list "getter-is-slice?" (getter-is-slice? left-hand-getter-info)
                       "all-derivatives-are-used?" (syntax->datum #'all-derivatives-are-used?))
            ((list "getter-is-slice?" #true
                   "all-derivatives-are-used?" #false)
             (assign-derivative-to-slice stx
                                         #'slice-idx
                                         #'slice-range
                                         #'mapped-idx
                                         #'mappings-full-idx
                                         #'al-index-symb
                                         #'mappings-synt
                                         #'der-vec-synt
                                         #'dx-range
                                         #'index-start-expanded_
                                         #'mappings-vec-len
                                         #'dual-b-derivative))
            ((list "getter-is-slice?" #false
                   "all-derivatives-are-used?" #false)
             (assign-derivative-to-cell stx
                                        #'mapped-idx
                                        index-exp
                                        #'mappings-full-idx
                                        #'al-index-symb
                                        #'mappings-synt
                                        #'der-vec-synt
                                        #'dx-range
                                        #'mappings-vec-len
                                        #'dual-b-derivative))
            ((list "getter-is-slice?" #false
                   "all-derivatives-are-used?" _)
             (assign-derivative-to-cell-dense stx
                                              #'mapped-idx
                                              index-exp
                                              #'mappings-full-idx
                                              #'al-index-symb
                                              #'mappings-synt
                                              #'der-vec-synt
                                              #'dx-range
                                              #'dual-b-derivative))
            ((list "getter-is-slice?" #true
                   "all-derivatives-are-used?" _)
             (assign-derivative-to-slice-dense stx
                                               #'slice-idx
                                               #'slice-range
                                               #'mapped-idx
                                               #'mappings-full-idx
                                               #'al-index-symb
                                               #'mappings-synt
                                               #'der-vec-synt
                                               #'dx-range
                                               #'index-start-expanded_
                                               #'mappings-vec-len
                                               #'dual-b-derivative))))
        (raise-syntax-error #f (format "bug: no mapping found for dual-l variable ~a" 
                                       (var-symbol-.name name-vs)) stx)))))


(define-for-syntax
  (assign-derivative-to-slice stx 
                   slice-idx
                   slice-range
                   mapped-idx
                   mappings-full-idx
                   al-index-symb
                   mappings-synt
                   der-vec-synt
                   dx-range
                   index-start-expanded_
                   mappings-vec-len
                   dual-b-derivative)
  #| (-> any/c  any/c  any/c  any/c  any/c |#
  #|     any/c  any/c  any/c  any/c  any/c |#
  #|     any/c |# 
  #|     syntax?) |#
 (datum->syntax
  stx
  `(_for ,slice-idx 0 ,slice-range
         (_forever ,mapped-idx 0
                   (_let-int ,mappings-full-idx
                             (_int+
                               (_var-ref ,mapped-idx)
                               (_int* ,dx-range (_int+ ,index-start-expanded_ (_var-ref ,slice-idx))))
                             (_if-stm
                               (_int>= (_var-ref ,mappings-full-idx) ,mappings-vec-len)
                               (_break)
                               (_let-int ,al-index-symb (_int-vector-ref ,mappings-synt
                                                                           (_var-ref ,mappings-full-idx))
                                         (_if-stm
                                           (_or (_int< (_var-ref ,al-index-symb) 0)
                                                (_int>= (_var-ref ,mapped-idx) ,dx-range))
                                           (_break)
                                           (_vector-set! ,der-vec-synt
                                                         (_int+
                                                           (_int* ,dx-range
                                                                  (_int+ ,index-start-expanded_
                                                                         (_var-ref ,slice-idx)))
                                                           (_var-ref ,mapped-idx))
                                                         ,dual-b-derivative)))))))))

(define-for-syntax (assign-derivative-to-slice-dense stx
                                                     slice-idx
                                                     slice-range
                                                     mapped-idx
                                                     mappings-full-idx
                                                     al-index-symb
                                                     mappings-synt
                                                     der-vec-synt
                                                     dx-range
                                                     index-start-expanded_
                                                     mappings-vec-len
                                                     dual-b-derivative)
(datum->syntax
  stx
  `(_for ,slice-idx 0 ,slice-range
         (_for ,mapped-idx 0 ,dx-range ;; TODO check if mapping block is equal to dx-range
               (_let-int ,mappings-full-idx
                         (_int+
                           (_var-ref ,mapped-idx)
                           (_int* ,dx-range (_int+ ,index-start-expanded_ (_var-ref ,slice-idx))))
                         (_let-int ,al-index-symb (_int-vector-ref ,mappings-synt
                                                                     (_var-ref ,mappings-full-idx))
                                   (_vector-set! ,der-vec-synt
                                                 (_int+
                                                   (_int* ,dx-range
                                                          (_int+ ,index-start-expanded_
                                                                 (_var-ref ,slice-idx)))
                                                   (_var-ref ,mapped-idx))
                                                 ,dual-b-derivative)))
               "vectorize"
               ))))

(define-for-syntax (assign-derivative-to-cell stx
                                              mapped-idx
                                              index-exp
                                              mappings-full-idx
                                              al-index-symb
                                              mappings-synt
                                              der-vec-synt
                                              dx-range
                                              mappings-vec-len
                                              dual-b-derivative)
  (datum->syntax
    stx
    ;; NOTE: it is valid for array not to have mapping at some indexes,
    ;;       because at some indexes derivatives are not needed. But 
    `(_forever ,mapped-idx 0
               (_let-int ,mappings-full-idx (_int+ (_var-ref ,mapped-idx) (_int* ,dx-range ,index-exp))
                         (_if-stm (_int>= (_var-ref ,mappings-full-idx) ,mappings-vec-len)
                                  (_break)
                                  (_let-int ,al-index-symb (_int-vector-ref ,mappings-synt
                                                                              (_var-ref ,mappings-full-idx))
                                            (_if-stm (_or (_int< (_var-ref ,al-index-symb) 0)
                                                          (_int>= (_var-ref ,mapped-idx) ,dx-range))
                                                     (_break)
                                                     (_vector-set!
                                                       ,der-vec-synt
                                                       (_int+ (_int* ,dx-range ,index-exp)
                                                              (_var-ref ,mapped-idx))
                                                       ,dual-b-derivative)))))
               ;; NOTE: dx-idxs-mapping al-index-symb conditions is not checked here
               ;; because in assgnation all indexses are exist
               )))


(define-for-syntax (assign-derivative-to-cell-dense stx
                                                    mapped-idx
                                                    index-exp
                                                    mappings-full-idx
                                                    al-index-symb
                                                    mappings-synt
                                                    der-vec-synt
                                                    dx-range
                                                    dual-b-derivative)
  (datum->syntax
    stx
    ;; NOTE: it is valid for array not to have mapping at some indexes,
    ;;       because at some indexes derivatives are not needed. But 
    `(_for ,mapped-idx 0 ,dx-range ;; TODO check if mapping block is equal to dx-range
           (_let-int ,mappings-full-idx (_int+ (_var-ref ,mapped-idx) (_int* ,dx-range ,index-exp))
                     (_let-int ,al-index-symb (_int-vector-ref ,mappings-synt
                                                                 (_var-ref ,mappings-full-idx))
                               (_vector-set!
                                 ,der-vec-synt
                                 (_int+ (_int* ,dx-range ,index-exp)
                                        (_var-ref ,mapped-idx))
                                 ,dual-b-derivative)))
           "vectorize"
           ;; NOTE: dx-idxs-mapping al-index-symb conditions is not checked here
           ;; because in assgnation all indexses are exist
           )))

(define/contract-for-syntax
  (search-left-hand-side-name stx ctx name)
  (-> (syntax/c any/c) func-context/c (syntax/c any/c)
      (values (syntax/c symbol?) type/c integer?))
  (let ((name_ (syntax->datum name))
        (func-name (func-context-.function-name ctx))
        (fake-src-pos 0))
    #| (displayln "ctx:") |#
    #| (pretty-print ctx) |#
    (cond
      ((hash-has-key? constants name_)
       (raise-syntax-error #f "assignation to a constant is not allowed" name))
      ((search-argument name_ (func-context-.current-arguments ctx))
       (raise-syntax-error #f "assignation to an argument is not allowed" name))
      ((equal? name_ func-name)
       (let ((func-return-value (func-context-.function-return-value ctx)))
         #| (displayln (format "name: ~a fake-src-pos is used" name)) |#
         (values
           (datum->syntax stx func-return-value)
           (func-context-.function-return-type ctx)
           fake-src-pos)))
      (else
        (let ((var (search-variable name_ (func-context-.current-variables ctx))))
          (cond
            (var
              (begin
                #| (displayln (format "name: ~a variable-src-pos = ~a is used." name (variable-src-pos var))) |#
                (values (datum->syntax stx (variable-symbol var))
                        (variable-type var)
                        (variable-src-pos var))))
            (else
              (raise-syntax-error #f "variable not found" name))))))))


(define/contract-for-syntax
  (make-inline-functions-stx stx inlined-function-list)
  (-> (syntax/c any/c) (listof function-inline-semantics-template/c)
      (listof (syntax/c any/c)))
  (for/list ((inl-f (in-list inlined-function-list)))
    (let* ((func-ret-symb (function-inline-semantics-template-.return-symbol inl-f))
           (func-body (function-inline-semantics-template-.body inl-f)) 
           (type (function-inline-semantics-template-.type inl-f))
           (function-name (function-inline-semantics-template-.name inl-f)))
      (with-syntax
        ((function-name function-name)
         (func-ret (datum->syntax stx func-ret-symb))
         (func-body func-body)
         (function-type (datum->syntax
                          stx ;; NOTE: unparse type back
                          (match type
                            ((list basic-type (list num))
                             `(array-type (basic-type ,(symbol->string basic-type)) "[" ,num "]"))
                            ((list basic-type (list))
                             `(basic-type ,(symbol->string basic-type)))))))

        ;; FIXME each funcition variable declaration should be in different namespace. it should be local.
        ; Fix it in expr-body
        ; z = f(x, y)
        ; h = f(x, y)
        ; expands to:
        ; new-variables-nesting
        ;   var-decl f
        ;   assignation f <- x y 
        ;   assignation z <- f 
        ; new-variables-nesting
        ;   var-decl f
        ;   assignation f <- x y 
        ;   assignation h <- f 
        ; but what if same funcition is called multiple times in the single assignation:
        ; z = g(f(x), f(y))
        ; function have distinct source positions. Can we include src-pos in the name?
        ; real-vars-table
        ;; TODO use syntax-track-origin to attach syntax-position to a function-name
        (datum->syntax
          stx
          `(_begin 
             ;; FIXME var-decl will add variable to the current-level, because it is
             ;; local-expand-memoed without new context
             (var-decl (type ,#'function-type) ,#'function-name)
             ,#'func-body))))))


(define-syntax (assignation stx)
  (begin
    #| (displayln (format "assignation time: ~a" (hash-ref TIME-TABLE 'assignation 0.0))) |#
    #| (displayln (format "assignation counts: ~a" (hash-ref TIME-TABLE 'assignation-counts 0))) |#
    #| (displayln (format "  'index-exp: ~a" (hash-ref TIME-TABLE 'index-exp 0))) |#
    #| (displayln (format "  'typecheck-mode: ~a" (hash-ref TIME-TABLE 'typecheck-mode  0))) |#
    #| (displayln (format "  'func-ret-assign: ~a" (hash-ref TIME-TABLE 'func-ret-assign 0))) |#
    #| (displayln (format "  'func-ret-assign_: ~a" (hash-ref TIME-TABLE 'func-ret-assign_ 0))) |#
    #| (displayln (format "  'value-exp: ~a" (hash-ref TIME-TABLE 'value-exp 0))) |#
    #| (displayln (format "  'search-left-hand-side-name ~a" (hash-ref TIME-TABLE 'search-left-hand-side-name 0))) |#
    #| (displayln (format "  'left-hand-getter-info ~a" (hash-ref TIME-TABLE 'left-hand-getter-info 0))) |#
    #| (displayln (format "  'get-slice-start-and-range ~a" (hash-ref TIME-TABLE 'get-slice-start-and-range 0))) |#
    #| (displayln (format "  'assignation-tail ~a" (hash-ref TIME-TABLE 'assignation-tail 0))) |#
    #| (displayln (format "    'make-der-assignation-loops-list ~a" (hash-ref TIME-TABLE 'make-der-assignation-loops-list 0))) |#
    #| (displayln (format "    'make-set-all-derivatives-to-zero-list ~a" (hash-ref TIME-TABLE 'make-set-all-derivatives-to-zero-list 0))) |#
    #| (displayln (format "    'make-der-assignation-loops-list/array ~a" (hash-ref TIME-TABLE 'make-der-assignation-loops-list/array 0))) |#
    #| (displayln (format "      'get-derivative-stx ~a" (hash-ref TIME-TABLE 'get-derivative-stx 0))) |#
    #| (displayln (format "      'get-derivative-stx counts: ~a" (hash-ref TIME-TABLE 'get-derivative-stx-counts 0))) |#
    #| (displayln (format "    'make-set-all-derivatives-to-zero/array ~a" (hash-ref TIME-TABLE 'make-set-all-derivatives-to-zero/array 0))) |#
    
    (timeit! TIME-TABLE 'assignation (thunk (syntax-parse stx

    (((~literal assignation) name:id
                             getter:getter-cls
                             ((~literal mut-assign) op) value:expr)
     (let* ((op_ (syntax->datum #'op))
            (getter-for-splice (syntax-e #'getter))
            (binop (match op_
                     ("+=" "+")
                     ("-=" "-")
                     ("*=" "*")
                     ("/=" "/"))))
           (match op_
             ((or "+=" "-=")
              (datum->syntax stx
                             `(assignation ,#'name ,@getter-for-splice "="
                                           (expr
                                             (expr 
                                               (term 
                                                 (factor 
                                                   (primary 
                                                     (element 
                                                       (get-value ,#'name ,@getter-for-splice)))))) ,binop ,#'value))))
             ((or "*=" "/=")
              (datum->syntax stx
                             `(assignation ,#'name ,@getter-for-splice "="
                                            (expr
                                              (term
                                                (term 
                                                  (factor 
                                                    (primary 
                                                      (element 
                                                        (get-value ,#'name ,@getter-for-splice))))) ,binop ,#'value))))))))

    (((~literal assignation) name:id 
                             getter:getter-cls
                             "=" value:expr)
     (let ((func-call-info-pair-list (make-state (list))))
      (with-syntax* 
        ((index-exp (timeit! TIME-TABLE 'index-exp (thunk (local-expand-memo #'getter.index 'expression '()))))
         ;; NOTE: expand to get the value-type. If value-type is dual then local-expand-memo it for each dx
         ;; Typecheck mode is also needed because before the first expansion some variables are
         ; not declared. For example, inlined function variable. This will cause `name not found`
         ; syntax error if expand without typecheck-mode.
         ;; FIXME uncomment
         (value-exp-typecheck-mode (timeit! TIME-TABLE 'typecheck-mode 
                                            (thunk (extract (local-expand-memo
                                                              #`(syntax-parameterize
                                                                  ((typecheck-mode #t)
                                                                   ;; NOTE: func-call populates this table 
                                                                   (func-call-box #,func-call-info-pair-list))
                                                                  value) 'expression (list)
                                                              #:reset-memo #t)))))
         )

        ;; NOTE: List of inlined functions called in the right-hand side of the assignation. 
        (define inlined-functions (timeit! TIME-TABLE 'func-ret-assign_
                                           (thunk (make-inline-functions-stx stx (read-state func-call-info-pair-list)))))
        ;; NOTE: local-expand-memo it to expand var-decl macro genetated for the function's return
        ; variable and add it to the current variables
        (define funcs-ret-assign 
          (timeit! TIME-TABLE 'func-ret-assign
                   (thunk (for/list ((item (in-list inlined-functions)))
                            (syntax-parse item
                              ((_begin func-var-decl body)
                               (with-syntax
                                 ((func-var-decl-exp (local-expand-memo #'func-var-decl 'expression (list #'define)
                                                                        #:reset-memo #t))
                                  (body-exp (local-expand-memo #'body 'expression (list #'define)
                                                               #:reset-memo #t)))
                                 (datum->syntax stx `(_begin ,#'func-var-decl-exp ,#'body-exp))
                                 )))))))
        (when (syntax->datum #'getter.index)
          (unless (equal? (syntax-property #'index-exp 'landau-type) 'int)
            (raise-syntax-error #f "index must be integer" #'getter.index)))
        (with-syntax*
          (
           ;; NOTE func-call-box is not used because it is already 
           ; populated in typecheck-mode expansion. 
           (value-exp (timeit! TIME-TABLE 'value-exp 
                               (thunk (match (target-lang TARGET)
                                        ('ansi-c
                                         (extract (local-expand-memo
                                                 #`(syntax-parameterize
                                                     ((expand-value-only #t))
                                                     value) 'expression (list)
                                                 #:reset-memo #t)))
                                        ('racket (local-expand-memo
                                                   #`(syntax-parameterize
                                                       ((expand-value-only #t))
                                                       (begin
                                                         ;; FIXME inlined function's arguments are unbond
                                                         #,@funcs-ret-assign
                                                         value)) 'expression (list #'begin))
                                         )
                                        ))))
           (value-exp-value-part #'value-exp))
          (let*-values
            (((ctx) (syntax-parameter-value #'func-context))
             ((name_) (syntax->datum #'name))
             ((name-str) (symbol->string name_))
             ((slice-colon_) (syntax->datum #'getter.slice-colon))
             ((value-type-hash-key) (syntax->hash-key #'value))
             ((value-type) (syntax-property #'value-exp-typecheck-mode 'landau-type))
             ((symbol full-type src-pos) (timeit/values! TIME-TABLE 'search-left-hand-side-name 
                                                         (thunk (search-left-hand-side-name stx ctx #'name))))
             ((left-hand-getter-info) (timeit! TIME-TABLE 'left-hand-getter-info
                                               (thunk (make-getter-info #'getter.index
                                                                        #'getter.slice-colon
                                                                        (to-landau-type stx full-type)))))
             ((name-vs) (var-symbol name-str src-pos))
             ((left-type) (car full-type))
             ((al-index-name_) 'al_index_name_symbol)
             ((maybe-array-range) (datum->syntax stx (cadr full-type)))
             ((index-start-expanded_ slice-range)
              (timeit/values! TIME-TABLE 'get-slice-start-and-range 
                              (thunk (get-slice-start-and-range 
                                       stx 
                                       slice-colon_ 
                                       (if slice-colon_ #'getter.index-start #'#f) (if slice-colon_ #'getter.index-end #'#f) maybe-array-range))))
             )
            #| (hash-set! VALUE-TYPES value-type-hash-key value-type) |#
            #| (pretty-print VALUE-TYPES) |#
            
            (when (and (equal? left-type 'int)
                       (equal? value-type 'real))
              (raise-syntax-error #f "assignation of real to integer is not allowed" stx))
            (timeit! 
              TIME-TABLE 
              'assignation-tail 
              (thunk 
                (with-syntax 
                  ((sym (datum->syntax stx symbol))
                   (type_ value-type)
                   (al-index-symb (datum->syntax stx al-index-name_))
                   (index-start-expanded_ index-start-expanded_)
                   (slice-range slice-range)
                   )
                  #| (displayln "FIXME: make sure that funcs-ret-assign should be always before set-all-derivatives-to-zero") |#
                  #| (displayln "FIXME: make sure that single function return variable is enough for sequential calls to one funcition") |#
                  (cond
                    ((getter-is-var? left-hand-getter-info)
                     ;; NOTE: Non-array case
                     (cond
                       ((equal? left-type 'dual-l)
                        (match value-type
                          ('dual-b
                           (let ((dx-names (check-result
                                             stx 
                                             (format "bug: assignation: need-derivatives-table-GLOBAL: no value for key: ~a" name-str)
                                             (ndt-get-dx-names (func-context-.need-derivatives-table ctx) name-vs))))
                             (with-syntax*
                               ((dual-b-value #'value-exp-value-part)
                                (assertion-loops-list
                                  (timeit! TIME-TABLE 'make-der-assignation-loops-list (thunk (make-der-assignation-loops-list stx ctx name-vs dx-names #'value #'al-index-symb)))
                                  ))
                               (datum->syntax stx
                                              `(expr-body
                                                 (_begin ,@funcs-ret-assign)
                                                 (_begin ,@#'assertion-loops-list)
                                                 (_set! ,#'sym ,#'dual-b-value))))))
                          ('real
                           (with-syntax
                             ((set-all-derivatives-to-zero
                                (timeit! TIME-TABLE 'make-set-all-derivatives-to-zero-list (thunk (make-set-all-derivatives-to-zero-list stx ctx name-vs #'al-index-symb)))
                                ))
                             (datum->syntax stx
                                            `(expr-body
                                               (_begin ,@funcs-ret-assign)
                                               (_begin ,@#'set-all-derivatives-to-zero)
                                               (_set! ,#'sym ,#'value-exp-value-part)))))
                          ('int
                           (with-syntax
                             ((set-all-derivatives-to-zero
                                (timeit! TIME-TABLE 'make-set-all-derivatives-to-zero-list (thunk (make-set-all-derivatives-to-zero-list stx ctx name-vs #'al-index-symb)))
                                ))
                             (datum->syntax stx
                                            `(expr-body
                                               (_begin ,@funcs-ret-assign)
                                               (_begin ,@#'set-all-derivatives-to-zero)
                                               (_set! ,#'sym (_exact->inexact ,#'value-exp-value-part))))))))
                       ;; NOTE: real variable that is not in need-only-value-set-GLOBAL is not used
                       ((equal? left-type 'real)
                        (match value-type
                          ;; FIXME: assignation to unused variable
                          ;  ((and 
                          ;    (not (set-member? (func-context-.need-only-value-set ctx) name-str))
                          ;    (not (equal? name_ func-name)))
                          ;   #'#f)
                          ((or 'dual-b 'real)
                           (datum->syntax stx
                                          `(expr-body
                                             (_begin ,@funcs-ret-assign)
                                             (_set! ,#'sym ,#'value-exp-value-part))))
                          ('int
                           (datum->syntax stx
                                          `(expr-body
                                             (_begin ,@funcs-ret-assign)
                                             (_set! ,#'sym (_exact->inexact ,#'value-exp-value-part)))))))
                       ;; NOTE: 'int <- 'int | 'real | 'dual-b
                       (else
                         (match value-type
                           ('int
                            (datum->syntax stx
                                           `(expr-body
                                              (_begin ,@funcs-ret-assign)
                                              (_set! ,#'sym ,#'value-exp-value-part))))
                           (_
                             (raise-syntax-error
                               #f (format "assignation to 'int is expected to be 'int. Given right-hand side type is ~a" value-type)
                               stx))))))
                    ;; NOTE: Array case
                    ((equal? left-type 'dual-l)
                     (begin
                       (unless (var-symbol/c name-vs)
                         (error "2569"))
                       (match (list (getter-info-type left-hand-getter-info) "<-" value-type)

                         ((or (list 'slice "<-" 'dual-b)
                              (list 'slice "<-" (list 'dual-b _))
                              (list 'cell  "<-" 'dual-b))
                          (let ((dx-names
                                  (cond
                                    ((ndt-member?
                                       (func-context-.need-derivatives-table ctx) name-vs)
                                     (ndt-get-dx-names (func-context-.need-derivatives-table ctx) name-vs))
                                    (else
                                      (raise-syntax-error
                                        #f
                                        (format "bug: assignation: need-derivatives-table: no value for key: ~a" name-vs) stx)))))
                            (with-syntax*
                              ((dual-b-value #'value-exp-value-part)
                               (assertion-loops-list
                                 (timeit! TIME-TABLE 'make-der-assignation-loops-list/array
                                  (thunk
                                   (make-der-assignation-loops-list/array stx
                                    ctx
                                    name-vs
                                    dx-names
                                    #'value
                                    #'al-index-symb
                                    left-hand-getter-info
                                    #'slice-range
                                    #'index-start-expanded_
                                    #'index-exp)))))
                              (if (getter-is-slice? left-hand-getter-info)
                                (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL)))
                                              (datum->syntax stx
                                                             `(expr-body
                                                                (_begin ,@funcs-ret-assign)
                                                                (_begin ,@#'assertion-loops-list)
                                                                (_for ,#'slice-idx 0 ,#'slice-range
                                                                      (_vector-set! ,#'sym (_int+ ,#'index-start-expanded_ (_var-ref ,#'slice-idx)) ,#'dual-b-value)))))
                                (datum->syntax stx
                                               `(expr-body
                                                  (_begin ,@funcs-ret-assign)
                                                  (_begin ,@#'assertion-loops-list)
                                                  (_vector-set! ,#'sym ,#'index-exp ,#'dual-b-value)))))))

                         ((or (list 'slice "<-" 'real)
                              (list 'slice "<-" (list 'real _))
                              (list 'cell  "<-" 'real))
                          (with-syntax
                            ((set-all-derivatives-to-zero 
                               (timeit! TIME-TABLE 'make-set-all-derivatives-to-zero/array (thunk (make-set-all-derivatives-to-zero/array stx
                                                                                                                                          ctx
                                                                                                                                          name-vs
                                                                                                                                          #'al-index-symb
                                                                                                                                          left-hand-getter-info
                                                                                                                                          #'slice-range
                                                                                                                                          #'index-start-expanded_
                                                                                                                                          #'index-exp)))
                               ))
                            (if (getter-is-slice? left-hand-getter-info)
                              (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL)))
                                            (datum->syntax stx
                                                           `(expr-body
                                                              (_begin ,@funcs-ret-assign)
                                                              (_begin ,@#'set-all-derivatives-to-zero)
                                                              (_for ,#'slice-idx 0 ,#'slice-range
                                                                    (_vector-set! ,#'sym (_int+ ,#'index-start-expanded_ (_var-ref ,#'slice-idx)) ,#'value)))))
                              (datum->syntax stx
                                             `(_local
                                                (expr-body
                                                  (_begin ,@funcs-ret-assign)
                                                  (_begin ,@#'set-all-derivatives-to-zero)
                                                  (_vector-set! ,#'sym ,#'index-exp ,#'value)))))))

                         ((or (list 'slice "<-" 'int)
                              (list 'cell  "<-" 'int))
                          (with-syntax
                            ((set-all-derivatives-to-zero
                               (timeit! TIME-TABLE 'make-set-all-derivatives-to-zero/array (thunk (make-set-all-derivatives-to-zero/array stx 
                                                                                                                                          ctx
                                                                                                                                          name-vs
                                                                                                                                          #'al-index-symb
                                                                                                                                          left-hand-getter-info
                                                                                                                                          #'slice-range
                                                                                                                                          #'index-start-expanded_
                                                                                                                                          #'index-exp)))
                               ))
                            (if (getter-is-slice? left-hand-getter-info)
                              (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL)))
                                            (datum->syntax stx
                                                           `(expr-body
                                                              (_begin ,@funcs-ret-assign)
                                                              (_begin ,@#'assertion-loops-list)
                                                              (_for ,#'slice-idx 0 ,#'slice-range
                                                                    (_vector-set! 
                                                                      ,#'sym 
                                                                      (_int+ ,#'index-start-expanded_ (_var-ref ,#'slice-idx)) 
                                                                      (_exact->inexact ,#'value))))))
                              (datum->syntax stx
                                             `(expr-body
                                                (_begin ,@funcs-ret-assign)
                                                (_begin ,@#'set-all-derivatives-to-zero)
                                                (_vector-set! ,#'sym ,#'index-exp (_exact->inexact ,#'value))))))))))

                    ((equal? left-type 'real)
                     (if (getter-is-slice? left-hand-getter-info)
                       ;; NOTE: lvalue is a slice
                       (with-syntax 
                         ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL)))
                         (cond
                           ;; FIXME: assignation to unused variable
                           ; ((and (not (set-member? (func-context-.need-only-value-set ctx) name-str))
                           ;       (not (equal? name_ func-name)))
                           ;  (with-syntax ((name-str name-str))
                           ;     (syntax/loc stx 
                           ;       (begin
                           ;        (flvector-set! sym index-exp value-exp)))))
                           ((or (equal? value-type 'dual-b)
                                (is-slice-of-type 'dual-b value-type)
                                (equal? value-type 'real)
                                (is-slice-of-type 'real value-type))
                            (datum->syntax stx ;; FIXME bug: in-range cause a parse error 
                                           `(_local
                                              (expr-body
                                                (_begin ,@funcs-ret-assign)
                                                (_for ,#'slice-idx 0 ,#'slice-range
                                                      (_vector-set! ,#'sym
                                                                    (_int+ ,#'index-start-expanded_ (_var-ref ,#'slice-idx))
                                                                    ,#'value-exp-value-part))))))
                           ((equal? value-type 'int)
                            (datum->syntax stx
                                           `(_local
                                              (expr-body
                                                (_begin ,@funcs-ret-assign)
                                                (_for ,#'slice-idx 0 ,#'slice-range
                                                      (_vector-set! ,#'sym
                                                                    (_int+ ,#'index-start-expanded_ (_var-ref ,#'slice-idx))
                                                                    (_exact->inexact ,#'value-exp-value-part)))))))))
                       ;; NOTE: lvalue is not a slice
                       (cond
                         ;; FIXME: assignation to unused variable
                         ; ((and (not (set-member? (func-context-.need-only-value-set ctx) name-str))
                         ;       (not (equal? name_ func-name)))
                         ;  (with-syntax ((name-str name-str))
                         ;     (syntax/loc stx 
                         ;       (begin
                         ;        (flvector-set! sym index-exp value-exp)))))

                         ((or (equal? value-type 'dual-b)
                              (equal? value-type 'real))
                          (datum->syntax stx
                                         `(expr-body
                                            (_begin ,@funcs-ret-assign)
                                            (_vector-set! ,#'sym ,#'index-exp ,#'value-exp-value-part))))
                         ((equal? value-type 'int)
                          (datum->syntax stx
                                         `(expr-body
                                            (_local (_begin ,@funcs-ret-assign))
                                            (_vector-set! ,#'sym ,#'index-exp (_exact->inexact ,#'value-exp-value-part)))))
                         ((is-slice? value-type)
                          (raise-syntax-error #f "Bug: assignation: no-a-slice <- slice" #'name)))))
                    ((equal? left-type 'int)
                     (raise-syntax-error #f "int arrays are not supported" #'name))
                    (else (raise-syntax-error 
                            #f 
                            (format "bug: assignation left-types are not matched: ~a <- ~a" left-type value-type)
                            #'name))))))))))))))))

(define-syntax (expr-body stx)
  (syntax-parse stx
    ((_ body ...)
     (let* ((ctx (syntax-parameter-value #'func-context))
            (new-vars-layer (new-variables-nesting
                              (func-context-.current-variables ctx)))
            (new-ctx (update-current-variables ctx new-vars-layer)))
      (quasisyntax/loc stx
                       (syntax-parameterize 
                         ((func-context #,new-ctx))
                         #,(datum->syntax stx
                                          `(_begin
                                             ,@#'(body ...)))))))))


(define-syntax (single-term stx)
  (syntax-parse stx
                [(_ body)
                 (syntax/loc stx body)]))

;; TODO: Maybe need syntax/loc instead of #'
(define-syntax (bool-expr stx)
  (syntax-parse stx
    [(_ b-expr1 "or" b-expr2) (datum->syntax stx `(_or ,#'b-expr1 ,#'b-expr2))]
    [(_ b-term) (syntax/loc stx b-term)]))
    

(define-syntax (bool-term stx)
  (syntax-parse stx
    [(_ b-term1 "and" b-term2) (datum->syntax stx `(_and b-term1 b-term2))]
    [(_ b-factor) (syntax/loc stx b-factor)]))

(define-syntax (bool-factor stx)
  (syntax-parse stx
    [(_ "(" b-expr ")") #'b-expr]
    [(_ "(" n1 ({~literal comp-op} op) n2 ")")
     (with-syntax* (
                    (val1 (local-expand-memo #'n1 'expression '()))
                    (val2 (local-expand-memo #'n2 'expression '())))
       (throw-if-not-int #'val1 "Only 'int allowed inside the `if` condition.")
       (throw-if-not-int #'val2 "Only 'int allowed inside the `if` condition.")
       (quasisyntax/loc stx #,(match (syntax->datum #'op)
                                ["==" (datum->syntax stx `(_equal? ,#'n1 ,#'n2))]
                                ["!=" (datum->syntax stx `(_not (_equal? ,#'n1 ,#'n2)))]
                                [">" (datum->syntax stx `(_int> ,#'n1 ,#'n2))]
                                [">=" (datum->syntax stx `(_int>= ,#'n1 ,#'n2))]
                                ["<=" (datum->syntax stx `(_int<= ,#'n1 ,#'n2))]
                                ["<" (datum->syntax stx `(_int< ,#'n1 ,#'n2))])))]
                    
    [(_ "not" b-factor) (datum->syntax stx `(_not ,#'b-factor))]))
    

; TODO: Do I need to check if it is not using outer variables?
(define-syntax (if-expr stx)
  (syntax-parse stx
    [({~literal if-expr} "if" b-expr "{" body-true ... "}" "else" "{" body-false ... "}")
     (datum->syntax
      stx
      `(_if-stm ,#'b-expr
                 (_begin ,@#'(body-true ...))
                 (_begin ,@#'(body-false ...))))]
         
    [({~literal if-expr} "if" b-expr "{" body ... "}")
     (datum->syntax
      stx
      `(_if-stm ,#'b-expr
                (_begin ,@#'(body ...))))]))
         


(define-syntax (for-expr stx)
  (syntax-parse stx
  [ (_ "for" ({~literal loop-var-decl} name "=" "[" start ":" stop "]") pat:expr-body-cls)
      (let ((name_ (syntax->datum #'name))
            (type (list 'int '())))
        (check-duplicate-variable-name name_ #'name)
        (let* ((ctx (syntax-parameter-value #'func-context))
               (new-vars-layer (new-variables-nesting
                                (func-context-.current-variables ctx)))
               (sym (add-variable! new-vars-layer #'name type))
               (new-ctx (update-current-variables ctx new-vars-layer)))
          (with-syntax* ([symm (datum->syntax stx sym)]
                         [start-val (local-expand #'start 'expression '())]
                         [stop-val (local-expand #'stop 'expression '())]
                         [return-stx (datum->syntax
                                      stx
                                      `(_for ,#'symm ,#'start-val ,#'stop-val ,#'pat.body))])
            (throw-if-not-int #'start-val "Only 'int allowed for ranges")
            (throw-if-not-int #'stop-val "Only 'int allowed for ranges")
            (quasisyntax/loc stx
                             (syntax-parameterize
                              [(func-context #,new-ctx)]
                              #,#'return-stx)))))]))
    


;; TODO Use getter-cls syntax-class to reduce match expression 
;; TODO Do not allocate vector of right-part lenght. Some values have no derivatives
(define-syntax (der-annot stx)
  (syntax-parse stx
    ((_ ({~literal get-value} df-name:id
                              (~optional
                               (~or*
                                (~seq "[" (~optional df-index-start #:defaults ((df-index-start #'#f))) 
                                      df-slice-colon:colon
                                      (~optional df-index-end #:defaults ((df-index-end #'#f))) "]")
                                (~seq "[" df-idx "]"))))
        "'"
        ({~literal get-value} dx-name:id
                              (~optional
                               (~or*
                                (~seq "[" (~optional dx-index-start #:defaults ((dx-index-start #'#f))) 
                                      dx-slice-colon:colon
                                      (~optional dx-index-end #:defaults ((dx-index-end #'#f))) "]")
                                (~seq "[" dx-idx "]"))))
        "="
        der-value)
     (begin
       (let* 
           ((df-type (get-real-arg-or-var-type #'df-name))
            (dx-type (get-real-arg-or-parameter-type #'dx-name))
            (df-array-range (datum->syntax stx (cadr df-type)))
            (dx-array-range (datum->syntax stx (cadr dx-type)))
            (df-name-symbol (syntax->datum #'df-name))
            (df-name-str (symbol->string df-name-symbol))
            (dx-name-str (symbol->string (syntax->datum #'dx-name)))
            (df-slice-colon_ (if (attribute df-slice-colon) #t #f))
            (dx-slice-colon_ (if (attribute dx-slice-colon) #t #f))
            (df-idx (if (attribute df-idx) #'df-idx #'#f))
            (dx-idx (if (attribute dx-idx) #'dx-idx #'#f))
            ;; NOTE: it is #f if curent dx is not needed for current df and der-annot is skipped then
            (ctx (syntax-parameter-value #'func-context)))
         (let*-values
           (((_1 _2 _3 df-src-pos) (search-name stx ctx df-name-symbol #'df-name #f #f #'#f))

            ((df-getter) (make-getter-info df-idx
                                           (if (attribute df-slice-colon) #'df-slice-colon #'#f)
                                           (to-landau-type stx df-type)))

            ((dx-getter) (make-getter-info dx-idx
                                           (if (attribute dx-slice-colon) #'dx-slice-colon #'#f)
                                           ;; NOTE: type is already expaded and partially evaluated.
                                           ; Some types are stored unexpanded and/or unevaluated and it
                                           ; is a mess, that should be resoled in future.
                                           ; All types should be stored expaded and evaluated to the landau-type.
                                           (to-landau-type stx dx-type)))

            ((df-name-vs) (var-symbol df-name-str df-src-pos))
            ((current-dx-maped-size) (get-der-variable-dx-range df-name-vs dx-name-str stx)))
          (cond
           ((and (equal? (car df-type) 'dual-l) current-dx-maped-size)
            (let*
                ((der-vec-symb (get-der-variable-symbol df-name-vs dx-name-str stx))
                 (df-name_ (syntax->datum #'df-name))
                 (dx-name_ (syntax->datum #'dx-name))
                 (df-size_ (cadr df-type)))
              (let*-values
                  (((dx-index-start-expanded dx-slice-range) (get-slice-start-and-range
                                                              stx 
                                                              dx-slice-colon_ 
                                                              (if dx-slice-colon_ 
                                                                #'dx-index-start 
                                                                #'#f) 
                                                              (if dx-slice-colon_ 
                                                                #'dx-index-end 
                                                                #'#f) 
                                                              dx-array-range))
                   ((df-index-start-expanded df-slice-range) (get-slice-start-and-range 
                                                              stx 
                                                              df-slice-colon_ 
                                                              (if df-slice-colon_ 
                                                                #'df-index-start 
                                                                #'#f) 
                                                              (if df-slice-colon_ 
                                                                #'df-index-end 
                                                                #'#f) 
                                                              df-array-range)))
                (with-syntax*
                    ((df-idx df-idx)
                     (dx-idx dx-idx)
                     (df-idx-expanded (local-expand #'df-idx 'expression '()))
                     (dx-idx-expanded (local-expand #'dx-idx 'expression '()))
                     (inv-mapping-period (get-inv-mapping-period df-name-vs dx-name-str stx))
                     (dx-idxs-mappings (get-dx-idxs-mappings-variable-symbol df-name-vs dx-name-str stx))
                     (dx-index-start-expanded dx-index-start-expanded)
                     (dx-slice-range dx-slice-range)
                     (df-index-start-expanded df-index-start-expanded)
                     (df-slice-range df-slice-range)
                     (value-exp (extract (local-expand-memo
                                          #`(syntax-parameterize
                                                ((expand-value-only #t))
                                              der-value) 'expression '())))
                     
                     (current-dx-maped-size current-dx-maped-size)
                     (maybe-mapped-idx (datum->syntax stx 'mapped_idx)))
                  (with-syntax
                      ((der-vec-synt (datum->syntax stx der-vec-symb)))
                    (match (list (getter-info-type df-getter)
                                 (getter-info-type dx-getter) 
                                 "<- der-value")
                      ((list 'var 'var "<- der-value")
                       (datum->syntax
                        stx
                        `(_local 
                           (_let-int 
                             ,#'maybe-mapped-idx (_int-vector-ref ,#'dx-idxs-mappings 0)
                            (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                     (_vector-set! ,#'der-vec-synt (_var-ref ,#'maybe-mapped-idx) 
                                                   ,#'value-exp))))))

                      ((list 'var 'cell "<- der-value")
                       (datum->syntax
                        stx
                        `(_local 
                           (_let-int 
                             ,#'maybe-mapped-idx (_if (_int< ,#'dx-idx-expanded ,#'inv-mapping-period)
                                                      (_int-vector-ref 
                                                        ,#'dx-idxs-mappings 
                                                        ,#'dx-idx-expanded)
                                                      -1)
                             (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                      (_vector-set! ,#'der-vec-synt (_var-ref ,#'maybe-mapped-idx) ,#'value-exp))))))

                      ((list 'var 'slice "<- der-value")
                       (datum->syntax
                        stx
                        `(_for ,#'slice_idx 0 ,#'dx-slice-range
                               (_let-int 
                                 ,#'maybe-mapped-idx 
                                 (_if (_int< (_int+ ,#'dx-index-start-expanded 
                                                    (_var-ref ,#'slice_idx)) 
                                             ,#'inv-mapping-period)
                                      (_int-vector-ref ,#'dx-idxs-mappings 
                                                       (_int+ ,#'dx-index-start-expanded (_var-ref ,#'slice_idx)))
                                      -1)
                                 (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                          (_vector-set! ,#'der-vec-synt (_var-ref ,#'maybe-mapped-idx) ,#'value-exp))))))

                      ((list 'cell 'var "<- der-value")
                       (datum->syntax 
                        stx
                        ;; NOTE: Take value, not der part of bundle
                        `(_local (_let-int ,#'maybe-mapped-idx (_int-vector-ref ,#'dx-idxs-mappings ,#'df-idx-expanded)
                                           (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                    (_vector-set! ,#'der-vec-synt
                                                                  ,#'df-idx-expanded ,#'value-exp))))))

                      ((list 'cell 'cell "<- der-value")
                       (datum->syntax 
                        stx
                        ;; NOTE: Take value, not der part of bundle
                        ;; NOTE: dx-idxs-mapping len is equal to (max_used_dx_index + 1)
                        ;; annotation dx_index can be out of range. Such annotations are not needed because
                        ;; they are not used later
                        `(_local (_let-int ,#'maybe-mapped-idx (_if (_int< ,#'dx-idx-expanded ,#'inv-mapping-period)
                                                                    (_int-vector-ref ,#'dx-idxs-mappings (_int+ (_int* ,#'df-idx-expanded ,#'inv-mapping-period) ,#'dx-idx-expanded))
                                                                    -1)
                                           (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                    (_vector-set! ,#'der-vec-synt
                                                                  (_int+ (_int* ,#'df-idx-expanded ,#'current-dx-maped-size) (_var-ref ,#'maybe-mapped-idx)) ,#'value-exp))))))

                      ((list 'cell 'slice "<- der-value")
                       (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                      (dx-full-idx (datum->syntax stx `(_int+ ,#'dx-index-start-expanded (_var-ref ,#'slice-idx)))))
                         (datum->syntax
                          stx
                          `(_for ,#'slice-idx 0 ,#'dx-slice-range
                                 (_let-int ,#'maybe-mapped-idx (_if (_int< ,#'dx-full-idx ,#'inv-mapping-period)
                                                                    (_int-vector-ref ,#'dx-idxs-mappings (_int+ (_int* ,#'df-idx-epping-period) ,#'dx-full-idx))
                                                                    -1)

                                           (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                    (_vector-set! ,#'der-vec-synt
                                                                  (_int+ (_int* ,#'df-idx-expanded ,#'current-dx-maped-size) (_var-ref ,#'maybe-mapped-idx)) ,#'value-exp)))))))

                      ((list 'slice 'var "<- der-value")
                       (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                      (df-full-idx (datum->syntax stx `(_int+ ,#'df-index-start-expanded (_var-ref ,#'slice-idx)))))
                         (datum->syntax
                          stx
                          `(_for ,#'slice-idx 0 ,#'df-slice-range
                                 (_let-int ,#'maybe-mapped-idx (_if (_int< 0 ,#'inv-mapping-period)
                                                                    (_int-vector-ref ,#'dx-idxs-mappings (_int+ (_int* ,#'df-full-idx ,#'inv-mapping-period) 0))
                                                                    -1)
                                           (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                    (_vector-set! ,#'der-vec-synt
                                                                  (_int+ (_int* ,#'df-full-idx ,#'current-dx-maped-size) (_var-ref ,#'maybe-mapped-idx)) ,#'value-exp)))))))

                      ((list 'slice 'cell "<- der-value")
                       (with-syntax* ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                      (df-full-idx (datum->syntax stx `(_int+ ,#'df-index-start-expanded (_var-ref ,#'slice-idx)))))
                         (datum->syntax
                          stx
                          `(_for ,#'slice-idx 0 ,#'df-slice-range
                                 (_let-int ,#'maybe-mapped-idx (_if (_int< ,#'dx-idx ,#'inv-mapping-period)
                                                                    (_int-vector-ref ,#'dx-idxs-mappings (_int+ (_int* ,#'df-full-idx ,#'inv-mapping-period) ,#'dx-idx))
                                                                    -1)
                                           (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                    (_vector-set! ,#'der-vec-synt
                                                                  (_int+ (_int* ,#'df-full-idx ,#'current-dx-maped-size) (_var-ref ,#'maybe-mapped-idx)) ,#'value-exp)))))))

                      ((list 'slice 'slice "<- der-value")
                       (with-syntax* ((dx-slice-idx (datum->syntax stx dx-slice-idx-name-GLOBAL))
                                      (df-slice-idx (datum->syntax stx df-slice-idx-name-GLOBAL))
                                      (slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                                      (df-full-idx (datum->syntax stx `(_int+ ,#'df-index-start-expanded (_var-ref ,#'slice-idx)))))
                         (datum->syntax
                          stx
                          `(_for ,#'df-slice-idx 0 ,#'df-slice-range
                                 (_for ,#'dx-slice-idx 0 ,#'dx-slice-range
                                       (_let-int ,#'maybe-mapped-idx (_if (_int< (_int+ ,#'dx-index-start-expanded (_var-ref ,#'dx-slice-idx)) ,#'inv-mapping-period)
                                                                          (_int-vector-ref ,#'dx-idxs-mappings (_int+ (_int* (_int+ ,#'df-index-start-expanded (_var-ref ,#'df-slice-idx)) ,#'inv-mapping-period) (_int+ ,#'dx-index-start-expanded (_var-ref ,#'dx-slice-idx))))
                                                                          -1)
                                                 (_if-stm (_int>= (_var-ref ,#'maybe-mapped-idx) 0)
                                                          (_let-int ,#'slice-idx (_int+ (_int* (_var-ref ,#'df-slice-idx) ,#'current-dx-maped-size) (_var-ref ,#'dx-slice-idx))
                                                                    (_vector-set! ,#'der-vec-synt
                                                                                  (_int+ (_int* (_int+ ,#'df-index-start-expanded (_var-ref ,#'df-slice-idx)) ,#'current-dx-maped-size) (_var-ref ,#'maybe-mapped-idx)) ,#'value-exp))))))))))
                    )))))
           (else (datum->syntax stx (_empty-statement))))))))))

;;FIXME: WIP: refactor #'func-name -> func-name-vs 
(define-syntax (der-apply stx)
  (syntax-parse stx
    [(_ ({~literal get-value} func-name:id
                              (~optional
                               (~or*
                                (~seq "[" (~optional func-index-start #:defaults ((func-index-start #'#f))) 
                                      func-slice-colon:colon
                                      (~optional func-index-end #:defaults ((func-index-end #'#f))) "]")
                                (~seq "[" func-ret-idx "]"))))
        "="
        value
        "'"
        ({~literal get-value} dx-name:id
                              (~optional
                               (~or*
                                (~seq "[" (~optional dx-index-start #:defaults ((dx-index-start #'#f)))
                                      dx-slice-colon:colon
                                      (~optional dx-index-end #:defaults ((dx-index-end #'#f))) "]")
                                (~seq "[" dx-idx "]")))))
     
     (let* ((ctx (syntax-parameter-value #'func-context))
            [func-name_ (syntax->datum #'func-name)]
            (func-name-str (symbol->string func-name_))
            [func-name (func-context-.function-name ctx)]
            [func-return-value (func-context-.function-return-value ctx)]

            (func-type (func-context-.function-return-type ctx))
            (dx-type (get-real-arg-or-parameter-type #'dx-name))
            (func-getter (make-getter-info (if (attribute func-ret-idx) #'func-ret-idx #'#f)
                                           (if (attribute func-slice-colon) #'func-slice-colon #'#f)
                                           (to-landau-type stx func-type)))
            (dx-getter (make-getter-info (if (attribute dx-idx) #'dx-idx #'#f)
                                         (if (attribute dx-slice-colon) #'dx-slice-colon #'#f)
                                         (to-landau-type stx dx-type)))

            (func-array-range (datum->syntax stx (cadr func-type)))
            (dx-array-range (datum->syntax stx (cadr dx-type)))
            (al-index-name_ 'al_index_name_symbol)
            (func-slice-colon_ (attribute func-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon))
            (fake-src-pos 0)
            [var-destination
             (if (equal? func-name_ func-name)
                 (datum->syntax stx func-return-value)
                 #f)]
            (throw-impossible-error (thunk
                                     (raise-syntax-error
                                      #f
                                      (format "bug: the right-hand side is slice, but the left-hand one is not") stx))))
       (unless var-destination (raise-syntax-error #f (format "name not found: ~v" #'func-name) stx))
       (with-syntax* ((func-ret-idx (if (attribute func-ret-idx) #'func-ret-idx #'#f))
                      (func-ret-slice-idx (datum->syntax stx 'func_slice_idx))
                      (dx-idx (if (attribute dx-idx) #'dx-idx #'#f))
                      (dx-name-str (symbol->string (syntax->datum #'dx-name)))
                      (expanded-value (extract (local-expand-memo
                                                #`(syntax-parameterize
                                                      ((typecheck-mode #t))
                                                    value) 'expression '())))
                      (df-name (syntax-property #'expanded-value 'get-value-name))
                      [var-destination-sym (datum->syntax stx var-destination)]
                      (debug-msg
                       (if debug
                           (displayln "der-apply: debug is not implemented")
                           (datum->syntax stx '(_nothing))))
                      (al-index-symb (datum->syntax stx al-index-name_)))
         ;(throw-if-not-type 'dual-b #'expanded-dual "This is a bug")
         ;(displayln #'expanded-dual)
         (let* ((df-type (syntax-property #'expanded-value 'landau-type)))
           (let*-values
               (((_1 _2 _3 df-src-pos) (search-name stx ctx (syntax->datum #'df-name) #'df-name #f #f #'#f))
                ((df-value-vs) (var-symbol (syntax->string #'df-name) df-src-pos))
                ((dx-index-start-expanded dx-slice-range) (get-slice-start-and-range
                                                           stx 
                                                           dx-slice-colon_ 
                                                           (if dx-slice-colon_ #'dx-index-start #'#f) 
                                                           (if dx-slice-colon_ #'dx-index-end #'#f) 
                                                           dx-array-range))
                ((func-index-start-expanded func-slice-range) (get-slice-start-and-range 
                                                               stx 
                                                               func-slice-colon_ 
                                                               (if func-slice-colon_ #'func-index-start #'#f)
                                                               (if func-slice-colon_ #'func-index-end #'#f) 
                                                               func-array-range)))
             (cond
               ((or (equal? df-type 'dual-b) (is-slice-of-type 'dual-b df-type))
                (with-syntax*
                    ((dx-index-start-expanded dx-index-start-expanded)
                     (dx-slice-range dx-slice-range)
                     (func-index-start-expanded func-index-start-expanded)
                     (func-slice-range func-slice-range)
                     (dual-b-derivative
                      (let ((dx-name-str (symbol->string (syntax->datum #'dx-name))))
                        (get-derivative-stx
                         (extract
                          (local-expand-memo
                           #`(syntax-parameterize
                                 ((dx-name-in-current-al '#,dx-name-str))
                               value) 'expression '()))))))
                  (match (list (getter-info-type func-getter) "<- df-value /" (getter-info-type dx-getter))
                    ((list 'var "<- df-value /" 'var)
                     (datum->syntax
                      stx
                      `(_local (_let-int ,#'al-index-symb 0
                                         ;; NOTE: al-index-symb is used inside dual-b-derivative
                                         (_set! ,#'var-destination-sym ,#'dual-b-derivative)))))

                    ((list 'var "<- df-value /" 'cell)
                     (datum->syntax
                      stx
                      `(_local (_let-int ,#'al-index-symb ,#'dx-idx
                                         ;; NOTE: al-index-symb is used inside dual-b-derivative
                                         (_set! ,#'var-destination-sym ,#'dual-b-derivative)))))

                    ((list 'var "<- df-value /" 'slice)
                     (throw-impossible-error))

                    ((list 'cell "<- df-value /" 'var)
                     (datum->syntax
                      stx
                      `(_local (_let-int ,#'al-index-symb 0
                                         ;; NOTE: al-index-symb is used inside dual-b-derivative
                                         (_vector-set! ,#'var-destination-sym ,#'func-ret-idx ,#'dual-b-derivative)))))

                    ((list 'cell "<- df-value /" 'cell)
                     (datum->syntax
                      stx
                      `(_local (_let-int ,#'al-index-symb ,#'dx-idx
                                         ;; NOTE: al-index-symb is used inside dual-b-derivative
                                         (_vector-set! ,#'var-destination-sym ,#'func-ret-idx ,#'dual-b-derivative)))))

                    ((list 'cell "<- df-value /" 'slice)
                     (throw-impossible-error))

                    ((list 'slice "<- df-value /" 'var)
                     (cond
                       ((is-slice? df-type)
                        (with-syntax*
                            ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range #'func-slice-range))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_let-int ,#'al-index-symb 0
                                            (_let-int ,#'func-ret-slice-idx (_int+ (_var-ref ,#'slice-idx) ,#'func-index-start-expanded)
                                                      (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative)))))))
                       (else
                        (with-syntax
                            ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range #'func-slice-range))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_let-int ,#'al-index-symb 0
                                            (_let-int ,#'func-ret-slice-idx (_int+ (_var-ref ,#'slice-idx) ,#'func-index-start-expanded)
                                                      (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative)))))))))

                    ((list 'slice "<- df-value /" 'cell)
                     (cond
                       ((is-slice? df-type)
                        (with-syntax*
                            ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range #'func-slice-range))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_let-int ,#'al-index-symb ,#'dx-idx
                                            (_let-int ,#'func-ret-slice-idx (_int+ (_var-ref ,#'slice-idx) ,#'func-index-start-expanded)
                                                      (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative)))))))
                       (else
                        (with-syntax*
                            ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range #'func-slice-range))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_let-int ,#'al-index-symb ,#'dx-idx
                                            (_let-int ,#'func-ret-slice-idx (_int+ (_var-ref ,#'slice-idx) ,#'func-index-start-expanded)
                                                      (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative)))))))))

                    ((list 'slice "<- df-value /" 'slice)
                     (cond
                       ((is-slice? df-type)
                        (with-syntax* 
                            ((dx-slice-idx (datum->syntax stx dx-slice-idx-name-GLOBAL))
                             (slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range (datum->syntax stx `(_int/ ,#'func-slice-range ,#'dx-slice-range)))
                             (dx-period (check-result
                                         stx
                                         (format "bug3: get-der-variable-dx-range retured #f for (~a, ~a)" 
                                                 df-value-vs 
                                                 (symbol->string (syntax->datum #'dx-name)))
                                         (get-der-variable-dx-range
                                          df-value-vs
                                          (symbol->string (syntax->datum #'dx-name))
                                          stx))))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_for ,#'dx-slice-idx 0 ,#'dx-slice-range
                                        (_let-int ,#'al-index-symb (_int+ (_var-ref ,#'dx-slice-idx) ,#'dx-index-start-expanded)
                                                  (_let-int ,#'func-ret-slice-idx (_int+ (_int* ,#'dx-period (_var-ref ,#'slice-idx)) (_var-ref ,#'dx-slice-idx) ,#'func-index-start-expanded)
                                                            (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative))))))))
                       (else
                        (with-syntax*
                            ((slice-idx (datum->syntax stx slice-idx-name-GLOBAL))
                             (value-slice-range #'func-slice-range))
                          (datum->syntax
                           stx
                           `(_for ,#'slice-idx 0 ,#'value-slice-range
                                  (_let-int ,#'al-index-symb (_int+ (_var-ref ,#'slice-idx) ,#'dx-index-start-expanded)
                                            (_let-int ,#'func-ret-slice-idx (_int+ (_var-ref ,#'slice-idx) ,#'func-index-start-expanded)
                                                      (_vector-set! ,#'var-destination-sym (_var-ref ,#'func-ret-slice-idx) ,#'dual-b-derivative))))))))))
                  ))
                    
               ((equal? df-type 'real) (datum->syntax stx `(_nothing)))
               (else (raise-syntax-error #f (format "only real variables or arguments can have derivatives. Given type: ~a" df-type) stx)))))))]))
               

(define-syntax (discard stx)
  (syntax-parse stx
    [(_ "discard" ({~literal get-value} df-name:id
                                        (~optional
                                         (~or*
                                          (~seq "[" (~optional df-index-start #:defaults ((df-index-start #'#f))) 
                                                df-slice-colon:colon
                                                (~optional df-index-end #:defaults ((df-index-end #'#f))) "]")
                                          (~seq "[" df-idx "]"))))
        "'"
        ({~literal get-value} dx-name:id
                              (~optional
                               (~or*
                                (~seq "[" (~optional dx-index-start #:defaults ((dx-index-start #'#f))) 
                                      dx-slice-colon:colon
                                      (~optional dx-index-end #:defaults ((dx-index-end #'#f))) "]")
                                (~seq "[" dx-idx "]")))))
     ;; NOTE: All logic is implemented in backrun
     (datum->syntax stx `(_empty-statement))]))

(define-syntax (print stx)
  (syntax-parse stx
    (({~literal print} "print" (~optional str #:defaults ((str #'debug))) expr)
     (with-syntax* ((line (syntax-line stx))
                    (str (symbol->string (syntax->datum #'str)))
                    (expr-expanded-typecheck (extract (local-expand-memo
                                                       #`(syntax-parameterize
                                                             ((typecheck-mode #t))
                                                           expr) 'expression '())))
                    (type (format "~a" (syntax-property #'expr-expanded-typecheck 'landau-type)))
                    (val (if (equal? #'type #'"dual-b") 
                               (extract (local-expand-memo
                                          #`(syntax-parameterize
                                                ((expand-value-only #t))
                                              expr) 'expression '()))
                               #'expr)))
       (when (equal? (target-lang TARGET) 'ansi-c)
         (raise-syntax-error #f "print is not implemented for 'ansi-c backend" stx))
       (quasisyntax/loc stx (displayln (format "~a: ~a type: ~a, value: ~a" line str type val)))))))

(define-for-syntax (check-duplicate-variable-name name stx-name)
  (check-duplicate-variable-name-helper name stx-name
                                        (func-context-.current-variables (syntax-parameter-value #'func-context))
                                        (func-context-.current-arguments (syntax-parameter-value #'func-context))
                                        (func-context-.function-name (syntax-parameter-value #'func-context))))

(define-syntax (get-derivative stx)
  (syntax-parse stx
    ((_ get-value-1 "'" get-value-2)
     (let* ((df-val (extract (local-expand-memo
                              #`(syntax-parameterize
                                    ((typecheck-mode #t))
                                  get-value-1) 'expression '())))
            ;  (dx-val-typecheck-mode (extract (local-expand-memo
            ;                                        #`(syntax-parameterize
            ;                                           ((typecheck-mode #t))
            ;                                           get-value-2) 'expression '())))
            (dx-val-get-name-mode (extract (local-expand-memo
                                            #`(syntax-parameterize
                                                  ((get-name-mode #t))
                                                get-value-2) 'expression '())))
            (df-type (syntax-property df-val 'landau-type))
            ; (dx-type (syntax-property dx-val-typecheck-mode 'landau-type))
            (df-name (syntax-property df-val 'get-value-name))
            (dx-name-str (fmap symbol->string (syntax-property dx-val-get-name-mode 'get-value-name)))
            (df-getter-info (syntax-property df-val 'getter-info))
            (dx-getter-info (syntax-property dx-val-get-name-mode 'getter-info))
            (df-name-vs (error "bug: get-derivative not implemented"))
            (dx-names (ndt-get-dx-names (func-context-.need-derivatives-table (syntax-parameter-value #'func-context)) df-name-vs))
            (dx-index (getter-info-ix-info dx-getter-info)))

       (when (getter-is-slice? df-getter-info)
         (raise-syntax-error #f "slices are not supported" stx))
       (when (getter-is-slice? dx-getter-info)
         (raise-syntax-error #f "slices are not supported" stx))
       (unless (and df-name dx-name-str)
         (raise-syntax-error #f "both df and dx in df ' dx should not be constants" stx))

       (match df-type
         ('dual-l
           (cond 
           ((member dx-name-str dx-names)
             (is-type_ 'dual-r 
                     (with-syntax
                         ((al-index (datum->syntax stx 'al_index_name_symbol))
                          (dual-b-derivative
                           (get-derivative-stx
                            (extract
                             (local-expand-memo
                              #`(syntax-parameterize
                                    ((dx-name-in-current-al '#,dx-name-str))
                                  get-value-1) 'expression '())))))
                       (datum->syntax
                        stx
                         `(_let-int ,#'al-index ,(if dx-index dx-index 0)
                            ,#'dual-b-derivative)))))
           (else (is-type_ 'message #'"derivative_not_set"))))
(_ (is-type_ 'message #'"has_no_derivatives"))
)))))

