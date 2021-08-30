#lang racket/base
;; TODO move the module to syntax-analyser directory
(require (for-syntax racket/base
                     racket/contract
                     racket/fixnum
                     racket/flonum
                     racket/function
                     racket/list
                     racket/match
                     racket/pretty
                     racket/set
                     racket/syntax
                     syntax/location
                     syntax/parse
                     "environment.rkt"
                     (only-in "process-actions-list.rkt" splice-nested)
                     )
         (only-in "process-actions-list.rkt" process!)
         "environment.rkt"
         "type-utils.rkt"
         racket/contract/region
         racket/function
         racket/match
         racket/stxparam
         racket/flonum
         racket/fixnum
         racket/unsafe/ops
         "common-for-syntax.rkt")


(define-syntax-parameter current-variables #f)
(define-syntax-parameter function-name #f)
(define-syntax-parameter function-return-value #f)
(define-syntax-parameter function-return-type  #f)
(define-syntax-parameter current-arguments #f)
(define-syntax-parameter dx-names-set #f)
(define-syntax-parameter real-vars-table #f)
(define-syntax-parameter module-funcs-table-parameter #f)

;; NOTE: Same logic as in semantics.rkt
(define-syntax-parameter func-call-box 'func-call-box-not-set)
(define-syntax-parameter func-is-called-inside-argument #f)

(define-for-syntax slice-idx-name-GLOBAL 'slice_idx)
(define-for-syntax df-slice-idx-name-GLOBAL 'df_slice_idx)
(define-for-syntax dx-slice-idx-name-GLOBAL 'dx_slice_idx)
(define-for-syntax func-slice-idx-name-GLOBAL 'func_slice_idx)
(define-for-syntax funcs-info-GLOBAL (make-hash))

(define func_slice_idx #f)
(define df_slice_idx #f)
(define dx_slice_idx #f)
(define slice_idx #f)


(begin-for-syntax
  (define-syntax-class plus-or-minus
    (pattern "+")
    (pattern "-"))
  (define-syntax-class mul-or-div
    (pattern "*")
    (pattern "/"))
  (define-syntax-class colon
    (pattern ":")))

(define-for-syntax (copy-syntax-properties src dst)
  (define props (syntax-property-symbol-keys src))
  (foldl 
    (lambda (prop-key result-stx)
      (define prop-value (syntax-property src prop-key))
      (syntax-property result-stx prop-key prop-value))
    dst
    props))

(define-for-syntax (local-expand-opaque stx)
  (local-expand stx 'expression '())
 #| (let-values |#
 #|   (((expanded opaque) (syntax-local-expand-expression stx))) |#
 #|   (define opaque-with-props (copy-syntax-properties expanded opaque)) |#
 #|   (with-syntax ((stx-opaque-with-props opaque-with-props) |#
 #|                 (stx-expanded expanded)) |#
 #|     (displayln (format "~a ~a ~a" stx expanded opaque-with-props)) |#
 #|     (displayln (format "~a ~a" (syntax-property-symbol-keys expanded) |#
 #|                        (syntax-property-symbol-keys opaque-with-props))) |#
 #|     (if (atom-number expanded) |#
 #|       expanded |#
 #|       opaque-with-props))) |#
 )


(define-for-syntax (evaluate expanded-syntax)
        (define ns (make-base-namespace))
        ;; NOTE: Racket generates some inner defined functions during the macro expansion.
        ;; They should be defined in the local namespace.
        (eval '(require (rename-in racket/private/list (reverse alt-reverse))) ns)
        (eval '(require racket/unsafe/ops) ns)
        (eval '(require racket/fixnum) ns)
        (eval '(define (check-range a b step)
                 (unless (real? a) (raise-argument-error 'in-range "real?" a))
                 (unless (real? b) (raise-argument-error 'in-range "real?" b))
                 (unless (real? step) (raise-argument-error 'in-range "real?" step))) ns)
        (eval (syntax->datum expanded-syntax) ns))


(define/contract-for-syntax
  (any->expanded-type v stx)
  (-> any/c syntax?
      landau-type/c)
  (cond
    ((syntax? v) (eval (syntax->datum (expand-type v))))
    (else (eval (syntax->datum (expand-type (datum->syntax stx v)))))))


(define-for-syntax (is-int-index stx)
  (is-type stx 'int-index))

(define-for-syntax (is-func? stx)
  (match (syntax->datum stx)
    ((list 'func _ ...) #t)
    (_ #f)))

(define/contract-for-syntax
  (range-casting-check stx r-value-type lvalue-outer-prod-range)
  (-> syntax? (list/c base-type/c (list/c (or/c syntax? fixnum?))) (or/c syntax? fixnum?)
      syntax?)
  (when (is-slice? r-value-type)
    (with-syntax ((rvalue-slice-range (get-slice-range r-value-type)))
      #`(unless (fx= #,lvalue-outer-prod-range #,#'rvalue-slice-range)
          (raise-syntax-error
            #f 
            (format 
              "can not cast right-hand side range ~v to the left-hand side range ~v" 
              #,#'rvalue-slice-range #,lvalue-outer-prod-range) #'#,stx)))))


(define-syntax (program stx)
  (syntax-parse stx
    [(_ body ...)
     (let* ((program-terms (syntax-e #'(body ...)))
            (top-level-decl (filter (compose not is-func?) program-terms))
            (funcs (filter is-func? program-terms))
            (fake-src-pos 0)
            (tld (if (empty? top-level-decl)
                   #'(void)
                   (local-expand-opaque (datum->syntax stx (cons 'begin top-level-decl))))))
       ;; NOTE: Evaluate top level constants and parameters declarations
       (eval tld (current-namespace))
       (datum->syntax
        stx
        (cons 'list 
              (list funcs-info-GLOBAL 
                    (cons 'list 
                          (let
                            #| module-funcs-table (-> (hash/c symbol? function-inline-template/c)) |#
                            ((module-funcs-table (make-hash)))
                            (for/list ((fnc (in-list funcs)))
                              (let ((fnc-name (symbol->string
                                                (syntax->datum (caddr (syntax-e fnc)))))
                                    (dx-names-set-val (make-hash))
                                    (real-vars-table-val (make-hash)))
                                #`(syntax-parameterize 
                                    ((dx-names-set '#,dx-names-set-val)
                                     (real-vars-table '#,real-vars-table-val)
                                     (module-funcs-table-parameter '#,module-funcs-table))
                                    (cons (var-symbol #,fnc-name #,fake-src-pos)
                                          (call-with-values
                                            (thunk 
                                              (process!
                                                #,fnc
                                                #,dx-names-set-val 
                                                #,real-vars-table-val))
                                            list)))))))))))]))
           
(define (throw-if-bundle-mistype bundl-1 bundl-2)
  (unless (equal? (car bundl-1) (car bundl-2))
    (error "bug: bundles have different type")))
  

(define-for-syntax (expand-range rng-stx rng)
  (if rng
      (begin
        (cadr (syntax->datum
               (local-expand-opaque (datum->syntax #'rng-stx rng)))))
      #f))

(define-for-syntax (check-proper-getter idx-stx slice-colon expanded-range expanded-idx name-str stx)
  (if expanded-range
      (begin
        (unless (or (syntax->datum idx-stx) slice-colon)
          (raise-syntax-error
           #f
           (format "Expect an index or a slice, ~a is an array." name-str) stx))
        (when expanded-idx (throw-if-not-type 'int-index expanded-idx "")))
      (begin
        (when (syntax->datum idx-stx)
          (raise-syntax-error
           #f
           (format "Do not expect index, ~a is not an array." name-str) stx))
        (when slice-colon
          (raise-syntax-error
           #f
           (format "Do not expect slice, ~a is not an array." name-str) stx)))))

(begin-for-syntax
  (define-syntax-class type-spec
    #:attributes (t)
    (pattern
     ((~literal type) t:expr)))
    

  (define-syntax-class arg-spec
    #:attributes (name type)
    (pattern
     (_ type:type-spec name:id))
     
    (pattern
     ({~literal other-argument} "," type:type-spec name:id))))
     

(define-syntax (func stx)
  (syntax-parse stx
    (({~literal func} type:expr name:id "(" ({~literal arglist} arg*:arg-spec ...) ")"
                      "{" body "}")
     (let* ((args (make-hash))
            (func-return-type (parse-type (syntax->datum #'type)))
            (func-name-str (symbol->string (syntax->datum #'name)))
            (func-return-value (gensym_ func-name-str))
            (fake-src-pos 0)
            (func-return-range (atom-number (type-range (expand-type (datum->syntax #'type func-return-type)))))
            (func-args
             (for/list ((argname (in-list (syntax-e #'(arg*.name ...))))
                        (argtype (in-list (syntax->datum #'(arg*.type ...)))))
                       (check-duplicate-argument-name args (syntax->datum #'name)
                                                      (syntax->datum argname) argname)
                       (datum->syntax argname (add-argument! args (syntax->datum argname)
                                                             (parse-type argtype)))
                       (define parsed-type (parse-type-to-syntax argtype))
                       (list (symbol->string (syntax->datum argname)) parsed-type)))
            (func-args-declaration (for/list ((name-type-pair (in-list func-args)))
                                     (match name-type-pair
                                       ((list name-str type)
                                        (begin
                                          (quasisyntax/loc
                                            stx
                                            (list 'var-decl 
                                                  '#,(var-symbol name-str fake-src-pos)
                                                  '#,(to-landau-type stx type))))))))
            (func-args-expanded
             (for/list ((arg (in-list func-args)))
                       (let* ((arg-name (car arg))
                              (arg-type (syntax-e (cadr arg)))
                              (arg-base-type (syntax->datum (cadr arg-type)))
                              (arg-range-expanded
                               (atom-number
                                (local-expand-opaque (datum->syntax stx (syntax->datum (caddr arg-type)))))))
                             (cons arg-base-type arg-range-expanded))))
            (func-parameters-info
             (for/list ((argname (in-list (syntax-e #'(arg*.name ...))))
                        (parsed-arg (in-list func-args)))
               (match parsed-arg
                 ((list _ argtype)
                  (list 
                    (var-symbol 
                      (symbol->string (syntax->datum argname))
                      ;; NOTE: use fake-src-pos because parameters always have unique names inside a single function
                      fake-src-pos)
                    (to-landau-type stx argtype))))))

            (func-return-vs (var-symbol func-name-str fake-src-pos)))

       (when (hash-has-key? BUILT-IN-FUNCTIONS func-name-str)
         (raise-syntax-error #f "function name shadows built-in function" #'name))

       (hash-set! ;; FIXME why use var-symbol with 0? Why dont we use just function name?
        funcs-info-GLOBAL
        func-return-vs 
        (func-info func-return-value (length func-args) func-args-expanded (car func-return-type) func-return-range))
      (with-syntax ((func-return-var-declaration-stx
                      (quasisyntax/loc 
                        stx
                        (list 'var-decl '#,func-return-vs '#,(to-landau-type stx func-return-type))))
                    (func-args-declaration-stx (quasisyntax/loc stx (list #,@func-args-declaration)))
                    (expanded-function-body 
                      (extract 
                        (local-expand
                          #`(syntax-parameterize
                              ((current-variables
                                 (new-variables-nesting
                                   (syntax-parameter-value #'current-variables)))
                               (current-arguments '#,args)
                               (function-name (syntax->datum #'name))
                               (function-return-value '#,func-return-value)
                               (function-return-type '#,func-return-type))
                              body) 'expression (list)))))
        ;; NOTE: Evaluate expanded function body to the nested lists of actions
        ;; TODO do not evaluate. Use semantics approach. write unexpanded stx instead of actions list
        (define func-body-evaluated (evaluate #'expanded-function-body))
        (hash-set! 
          (syntax-parameter-value #'module-funcs-table-parameter) 
          (syntax->datum #'name)
          (function-inline-template func-body-evaluated
                                    func-parameters-info
                                    func-return-vs))

        (quasisyntax/loc 
          stx 
          (list 
            #,#'func-args-declaration-stx 
            #,#'func-return-var-declaration-stx
            #,#'expanded-function-body)))))))


(define/contract-for-syntax
  (search-name-backrun stx name)
  (-> (syntax/c any/c) (syntax/c symbol?)
      ;; NOTE: real constants are not propagated in backrun
      (values (syntax/c (or/c symbol? integer? false/c)) type/c kind/c integer?))
  (let
    ((fake-src-pos 0)
     (is-const #t)
     (is-not-const #f)
     (name_ (syntax->datum name))
     (func-name (syntax-parameter-value #'function-name)))
    (let-values
      (((value type kind src-pos)
        (cond
          ((hash-has-key? constants name_)
           (let ((c (hash-ref constants name_)))
             (values
               (datum->syntax stx
                              (constant-value c))
               (constant-type c)
               'constant
               fake-src-pos)))
          (else
            (begin
              (unless (syntax-parameter-value #'current-variables)
                (raise-syntax-error #f "name not found" name))
              (let ((var (search-variable
                           name_ (syntax-parameter-value #'current-variables))))
                (cond
                  (var
                    (values (datum->syntax stx (variable-symbol var))
                            (variable-type var)
                            'variable
                            (variable-src-pos var)))
                  ((equal? name_ func-name)
                   (let ((func-return-value (syntax-parameter-value #'function-return-value)))
                     (values
                       (datum->syntax stx func-return-value)
                       (syntax-parameter-value #'function-return-type)
                       'function
                       fake-src-pos)))
                  (else
                    (let ((arg (search-argument
                                 name_ (syntax-parameter-value #'current-arguments))))
                      (cond
                        (arg
                          (values
                            (datum->syntax stx (argument-symbol arg))
                            (argument-type arg)
                            'argument
                            fake-src-pos))
                        (else
                          (cond
                            ((hash-has-key? parameters name_) 
                             (let* ((par-size (hash-ref parameters name_))
                                    (par-type (match par-size
                                                ((list size)
                                                 (landau-type 'real size))
                                                ((list)
                                                 (landau-type 'real)))))
                               (values
                                    (datum->syntax stx name_)
                                    par-type
                                    'argument
                                    fake-src-pos))
                             ; (raise-syntax-error #f "parameters can not be used in expression" name)
                             )
                            (else (raise-syntax-error #f "name not found" name))))))))))))))
      ; (println (format "~a -> ~a" name value))
      (values value type kind src-pos))

    ))


;; TODO: Check types. They are available after zero'th compiler run
;; NOTE: Now assume, that all get-value's are real (dual potentially) variables
(define-syntax (get-value stx)
  (syntax-parse stx
    (({~literal get-value} name:id
                           (~optional
                            (~or*
                             (~seq "[" (~optional index-start #:defaults ((index-start #'#f))) 
                                   slice-colon:colon
                                   (~optional index-end #:defaults ((index-end #'#f))) "]")
                             (~seq "[" index "]"))))
     (let*-values
         (((name_) (syntax->datum #'name))
          ((name-str) (symbol->string name_))
          ((index_) (if (attribute index) (syntax->datum #'index) #f))
          ((slice-colon_) (attribute slice-colon))
          ;; NOTE: stx-pos is provided for variables only
          ((value full-type kind src-pos) (search-name-backrun stx #'name))
          ((const?) (equal? kind 'constant))
          ((getter) (let ()
                      (define _type (expand-type-to-datum (datum->syntax stx full-type)))
                      (make-getter-info (if (attribute index) #'index #'#f)
                                      (if (attribute slice-colon) #'slice-colon #'#f)
                                      _type)))

          ((maybe-array-range) (cadr full-type))
          ((base-type) (car full-type))

          ((array-range) (atom-number #`#,(expand-range 
                                            (if (attribute index) #'index #'#f)
                                            (if (equal? maybe-array-range '()) #f maybe-array-range))))
          ((index-expanded) (if (attribute index) 
                              (local-expand-opaque #'index)
                              #f)))
         (when (equal? kind 'parameter)
           (raise-syntax-error #f "parameters can not be used in expression" name_))
         ;;FIXME pass getter-info as syntax property?
         #| (displayln "FIXME: check getter in get-value") |# 
         #| (displayln (format "~a get-value: name: ~a ~a" (syntax-line #'name) name-str (syntax->datum stx))) |#
      ;; FIXME Do not fail array get-value without index of slice. Array can be passed to a function: f(arr)
       #| (check-proper-getter |# 
       #|  (if (attribute index) #'index #'#f) (attribute slice-colon) array-range index-expanded name-str stx) |#
       #| (displayln "full-type:") |#
       #| (displayln (expand-type (datum->syntax stx full-type))) |#
       (with-syntax-property 
         'full-type (to-landau-type stx full-type)
         (with-syntax-property 
           'get-value-name name_
           (with-syntax-property 
             'getter-info getter 
             (with-syntax ((value value)
                           (name-symb (var-symbol (symbol->string name_) src-pos)))
               (let*-values 
                 (((index-start-stx slice-range index-range-checks)
                   (if slice-colon_
                     (slice-helper stx #'index-start #'index-end array-range)
                     (values #f #f #f)))
                  ((result)
                   (match (list (getter-info-type getter) base-type)
                     ((list 'var 'real)
                      (if const?
                        (is-real (syntax/loc stx (list)))
                        (is-real (syntax/loc stx (list (list 'array-ref name-symb 0))))))

                     ((list 'array 'real)
                      (with-syntax-property 
                        'landau-type (to-landau-type stx full-type)
                        (if const?
                          (syntax/loc stx (list))
                          ;; NOTE: idx do not make sense and should not be used. Use -1 to indicate this. 
                          (syntax/loc stx (list (list 'array-ref name-symb -1))))))

                     ((list 'cell 'real)
                      (with-syntax* ((index-expanded-stx index-expanded)
                                     (expanded-range array-range)
                                     (range-check-stx (range-check-stx stx #'index-expanded-stx #'expanded-range)))

                                    ;; NOTE: (cadr (cadr ... because expanded looks like that: '(#%app '4)
                                    ;; NOTE: can't do check in transmormer phase because loop variable is not resolved by that time.
                                    ;;       Need to insert check-bounds code in generated code
                                    (throw-if-not-type 'int-index #'index-expanded-stx 
                                                       "Indexes expected to be expressions with 'int constants and/or loop variables")
                                    (begin
                                      (is-real
                                        (if const?
                                          (quasisyntax/loc
                                            stx
                                            (begin
                                              range-check-stx
                                              (list)))

                                          (quasisyntax/loc
                                            stx
                                            (begin
                                              range-check-stx ;; TODO: add src-loc
                                              (list (list 'array-ref #,#'name-symb #,#'index-expanded-stx))
                                              (list (list 'array-ref #,#'name-symb #,#'index-expanded-stx)))))))))

                     ((list 'slice 'real)
                      (with-syntax*
                        ((slice-range-stx slice-range)
                         (index-range-checks-stx index-range-checks)
                         (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                         (full-idx #`(+ #,index-start-stx slice-idx-synt)))
                        (is-type_ (landau-type 'real #'slice-range-stx)
                                  (if const?
                                    (quasisyntax/loc
                                      stx
                                      (begin
                                        #,#'index-range-checks-stx
                                        (list)))
                                    (quasisyntax/loc
                                      stx
                                      (begin
                                        #,#'index-range-checks-stx
                                        (list (list 'array-ref name-symb full-idx))))))))

                     ((list 'var 'int)
                      (is-int (syntax/loc stx '())))

                     ((list (or 'cell 'slice) (or 'int 'int-index))
                      (raise-syntax-error stx "int arrays are not supported yet" #'name-symb))

                     ((list 'var 'int-index)
                      (is-int-index (syntax/loc stx value))))
                   ))
                 result)))))))))

(define-for-syntax (range-check stx idx range)
  (begin
   (println idx)
   (println range)
   (when (fx>= idx range)
     (raise-syntax-error
      #f
       (format "index ~a out of range: expect [0, ~a)" idx range) stx))))

(define-for-syntax (range-check-stx stx idx range)
  #`(when (fx>= #,idx #,range)
    (raise-syntax-error
     #f
      (format "index ~a out of range: expect [0, ~a)" #,idx #,range) #'#,stx)))

;; TODO: refactor to reduce bindings duplication
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
     (let* ((name-1-dat (symbol->string (syntax->datum #'func-name)))
            (name-2-dat (symbol->string (syntax->datum #'dx-name)))
            (arr-1-range (check-real-func-existence #'func-name 'get-range))
            (arr-2-range (check-real-arg-or-parameter-existence #'dx-name 'get-range))
            (func-type (search-function #'func-name))
            (dx-type (search-arg-or-parameter #'dx-name))
            (func-getter (make-getter-info (if (attribute func-ret-idx) #'func-ret-idx #'#f)
                                           (if (attribute func-slice-colon) #'func-slice-colon #'#f)
                                           (to-landau-type stx func-type)))
            (dx-getter (make-getter-info (if (attribute dx-idx) #'dx-idx #'#f) 
                                         (if (attribute dx-slice-colon) #'dx-slice-colon #'#f) 
                                         (to-landau-type stx dx-type)))
            (func-ret-idx-expanded (if (attribute func-ret-idx) 
                                     (local-expand-opaque #'func-ret-idx)
                                     #f))
            (dx-idx-expanded (if (attribute dx-idx)
                               (local-expand-opaque #'dx-idx)
                               #f))
            (expanded-range 
              (atom-number 
                #`#,(expand-range (if (attribute func-ret-idx) #'func-ret-idx #'#f) 
                                  (if (equal? arr-1-range '()) #f arr-1-range))))
            (expanded-range-dx 
              (atom-number 
                #`#,(expand-range (if (attribute dx-idx) #'dx-idx #'#f) 
                                  (if (equal? arr-2-range '()) #f arr-2-range))))
            (value-type (syntax-property 
                          (local-expand-opaque #'value)
                          'landau-type))
            (func-slice-colon_ (attribute func-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon))
            (throw-slice-cast-error (thunk
                                      (raise-syntax-error
                                        #f
                                        (format "the right-hand side is slice, but the left-hand one is not") stx))))
       (check-proper-getter (if (attribute func-ret-idx) #'func-ret-idx #'#f) 
                            (attribute func-slice-colon) 
                            expanded-range 
                            func-ret-idx-expanded 
                            name-1-dat 
                            stx)
       (check-proper-getter (if (attribute dx-idx) #'dx-idx #'#f) 
                            (attribute dx-slice-colon) 
                            expanded-range-dx 
                            dx-idx-expanded
                            name-2-dat
                            stx)
       (let*-values
         (((func-index-start-stx func-slice-range func-index-range-checks)
           (if func-slice-colon_
             (slice-helper stx #'func-index-start #'func-index-end expanded-range)
             (values #f #f #f)))
          ((dx-index-start-stx dx-slice-range dx-index-range-checks)
           (if dx-slice-colon_
             (slice-helper stx #'dx-index-start #'dx-index-end expanded-range-dx)
             (values #f #f #f))))
         (with-syntax* ((dx-name-str name-2-dat) ;; NOTE: Do not need src-pos for dx
                        (func-index-range-checks func-index-range-checks)
                        (dx-index-range-checks dx-index-range-checks))
                       (match (list (getter-info-type func-getter) "<- value /" (getter-info-type dx-getter))
                         ((list 'var "<- value /" 'var)
                          (begin
                            (throw-if-r-value-is-slice stx #'value)
                            (quasisyntax/loc stx
                                             (list 'der-apply
                                                   (car #,#'value) ;; get-value according to the grammar
                                                   (list 'array-ref #,#'dx-name-str 0)))))

                         ((list 'var "<- value /" 'cell)
                          (begin
                            (throw-if-r-value-is-slice stx #'value)
                            (quasisyntax/loc stx
                                             (begin
                                               #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                                               (list 'der-apply
                                                     (car #,#'value) ;; get-value according to the grammar
                                                     (list 'array-ref #,#'dx-name-str #,dx-idx-expanded))))))

                         ((list 'var "<- value /" 'slice)
                          (throw-slice-cast-error))

                         ((list 'cell "<- value /" 'var)
                          (begin
                            (throw-if-r-value-is-slice stx #'value)
                            (quasisyntax/loc stx
                                             (begin
                                               #,(range-check-stx stx func-ret-idx-expanded expanded-range)
                                               (list 'der-apply
                                                     (car #,#'value) ;; get-value according to the grammar
                                                     (list 'array-ref #,#'dx-name-str 0))))))

                         ((list 'cell "<- value /" 'cell)
                          (begin
                            (throw-if-r-value-is-slice stx #'value)
                            (quasisyntax/loc stx
                                             (begin
                                               #,(range-check-stx stx func-ret-idx-expanded expanded-range)
                                               #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                                               (list 'der-apply
                                                     (car #,#'value) ;; get-value according to the grammar
                                                     (list 'array-ref #,#'dx-name-str #,dx-idx-expanded))))))

                         ((list 'cell "<- value /" 'slice)
                          (throw-slice-cast-error))

                         ((list 'slice "<- value /" 'var)
                          (with-syntax*
                            ((value-slice-range (coerce-to-fixnum (get-slice-range value-type)))
                             (func-slice-range func-slice-range)
                             (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                             (rvalue-outer-prod-range #'value-slice-range)
                             (r-value-type value-type))
                            #| (displayln (format "slice <- value / var: ~a" ) ) |#
                            (quasisyntax/loc
                              stx
                              (begin
                                #,#'func-index-range-checks
                                (unless (fx= rvalue-outer-prod-range 1)
                                  (unless (fx= rvalue-outer-prod-range func-slice-range)
                                    (raise-syntax-error 
                                      #f 
                                      (format "can not cast right-hand side range ~v to the left-hand side range ~v" 
                                              rvalue-outer-prod-range 
                                              func-slice-range) 
                                      #'#,stx)))
                                (for/list ((#,#'slice-idx-synt (in-range 0 #,#'func-slice-range)))
                                  (list 'der-apply
                                        (car #,#'value)
                                        (list 'array-ref #,#'dx-name-str 0)))))))

                         ((list 'slice "<- value /" 'cell)
                          (with-syntax*
                            ((dx-idx-expanded dx-idx-expanded)
                             (func-slice-range func-slice-range)
                             (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                             (full-idx #`(+ #,func-index-start-stx slice-idx-synt))
                             (value-slice-range (coerce-to-fixnum (get-slice-range value-type))))
                            (quasisyntax/loc
                              stx
                              (begin
                                #,#'func-index-range-checks
                                (unless value-slice-range
                                  (unless (fx= value-slice-range func-slice-range) ;; Bug is here
                                    (raise-syntax-error 
                                      #f 
                                      (format "can not cast right-hand side range ~v to the left-hand side range ~v" 
                                              value-slice-range 
                                              func-slice-range) 
                                      #'#,stx)))
                                (for/list ((#,#'slice-idx-synt (in-range 0 #,#'func-slice-range)))
                                  (list 'der-apply
                                        (car #,#'value)
                                        (list 'array-ref #,#'dx-name-str #,#'dx-idx-expanded)))))))

                         ((list 'slice "<- value /" 'slice)
                          (with-syntax*
                            ((func-slice-range func-slice-range)
                             (dx-slice-range dx-slice-range)
                             (value-slice-range (coerce-to-fixnum (get-slice-range value-type)))
                             (rvalue-outer-prod-range #'(fx* value-slice-range dx-slice-range))
                             (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                             (dx-slice-idx-synt (datum->syntax stx dx-slice-idx-name-GLOBAL))
                             (dx-full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt)))
                            (quasisyntax/loc
                              stx
                              (begin
                                #,#'func-index-range-checks
                                #,#'dx-index-range-checks
                                (unless (fx= rvalue-outer-prod-range func-slice-range)
                                  (raise-syntax-error 
                                    #f 
                                    (format "can not cast right-hand side range ~v to the left-hand side range ~v" 
                                            rvalue-outer-prod-range 
                                            func-slice-range) 
                                    #'#,stx))
                                (for/list ((#,#'slice-idx-synt (in-range 0 #,#'value-slice-range)))
                                  (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                    (list 'der-apply
                                          (car #,#'value)
                                          (list 'array-ref #,#'dx-name-str #,#'dx-full-idx)))))))))
                       )))]))

(define-for-syntax (coerce-to-fixnum x)
  (if (equal? #f x)
   1
   x))

;; TODO: refactor to reduce code duplication
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
     (let* ((name-1-dat (symbol->string (syntax->datum #'df-name)))
            (name-2-dat (symbol->string (syntax->datum #'dx-name)))
            (df-range (check-real-var-existence #'df-name 'get-range))
            (dx-range (check-real-arg-or-parameter-existence #'dx-name 'get-range))
            (df-idx-expanded (if (attribute df-idx) 
                               (local-expand-opaque #'df-idx)
                               #f))
            (dx-idx-expanded (if (attribute dx-idx)
                               (local-expand-opaque #'dx-idx)
                               #f))
            (expanded-range-df (atom-number 
                                 #`#,(expand-range 
                                       (if (attribute df-idx) #'df-idx #'#f) 
                                       (if (equal? df-range '()) #f df-range))))
            (expanded-range-dx (atom-number 
                                 #`#,(expand-range 
                                       (if (attribute dx-idx) #'dx-idx #'#f) 
                                       (if (equal? dx-range '()) #f dx-range))))
            (df-slice-colon_ (attribute df-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon)))

       (check-proper-getter 
         (if (attribute df-idx) #'df-idx #'#f) 
         (attribute df-slice-colon) 
         expanded-range-df 
         df-idx-expanded 
         name-1-dat 
         stx)
       (check-proper-getter 
         (if (attribute dx-idx) #'dx-idx #'#f) 
         (attribute dx-slice-colon) 
         expanded-range-dx 
         dx-idx-expanded 
         name-2-dat 
         stx)
       
       (let*-values (((v0 v1 v2 df-name-src-pos) (search-name-backrun stx #'df-name))
                     ((v3 v4 v5 dx-name-src-pos) (search-name-backrun stx #'dx-name))
                     ((dx-index-start-stx dx-slice-range dx-index-range-checks)
                      (if dx-slice-colon_ 
                          (slice-helper stx #'dx-index-start #'dx-index-end expanded-range-dx)
                          (values #f #f #f)))
                     ((df-index-start-stx df-slice-range df-index-range-checks)
                      (if df-slice-colon_ 
                          (slice-helper stx #'df-index-start #'df-index-end expanded-range-df) 
                          (values #f #f #f))))
         (with-syntax* ((df-name-symb (var-symbol name-1-dat df-name-src-pos))
                        (dx-name-symb (var-symbol name-2-dat dx-name-src-pos))
                        (df-slice-range df-slice-range)
                        (dx-slice-range dx-slice-range)
                        (dx-index-range-checks dx-index-range-checks)
                        (df-index-range-checks df-index-range-checks))
           (hash-set! (syntax-parameter-value #'dx-names-set) (symbol->string (syntax->datum #'dx-name)) #t)
           (cond
             ((and (attribute df-idx) (attribute dx-idx))
              (quasisyntax/loc stx
                (begin
                  #,(range-check-stx stx df-idx-expanded expanded-range-df)
                  #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                  (when (and (equal? df-name-symb dx-name-symb) (equal? #,df-idx-expanded #,dx-idx-expanded))
                    (raise-syntax-error
                     #f
                     (format "Can not override annotation: ~a[~a] ' ~a[~a] == 1.0" 
                             df-name-symb #,df-idx-expanded dx-name-symb #,dx-idx-expanded) #'#,stx))
                  (list 'discard
                        (list 'array-ref #,#'df-name-symb #,df-idx-expanded)
                        (list 'array-ref #,#'dx-name-symb #,dx-idx-expanded)))))
             ((and (attribute df-idx) (not (attribute dx-idx)))
              (cond
                (dx-slice-colon_
                 (with-syntax*
                     ((df-idx-expanded df-idx-expanded) 
                      (dx-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'dx-slice-range)
                      (r-value-exp 
                        (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                         (list 'discard
                               (list 'array-ref #,#'df-name-symb #,#'df-idx-expanded)
                               (list 'array-ref #,#'dx-name-symb #,#'full-idx)))))))
                (else
                 (begin
                   (throw-if-r-value-is-slice stx #'der-value)
                   (quasisyntax/loc stx
                     (begin
                       #,(range-check-stx stx df-idx-expanded expanded-range-df)
                       (list 'discard
                             (list 'array-ref #,#'df-name-symb #,df-idx-expanded)
                             (list 'array-ref #,#'dx-name-symb 0))))))))
                          
             ((and (not (attribute df-idx)) (attribute dx-idx))
              (cond
                (df-slice-colon_
                 (with-syntax*
                     ((dx-idx-expanded dx-idx-expanded) 
                      (df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp 
                        (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                         (list 'discard
                               (list 'array-ref #,#'df-name-symb #,#'full-idx)
                               (list 'array-ref #,#'dx-name-symb #,#'dx-idx-expanded)))))))
                (else (quasisyntax/loc stx
                        (begin
                          #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                          (list 'discard
                                (list 'array-ref #,#'df-name-symb 0)
                                (list 'array-ref #,#'dx-name-symb #,dx-idx-expanded)))))))
             ((and (not (attribute df-idx)) (not (attribute dx-idx)))
              (cond
                ((and df-slice-colon_ dx-slice-colon_) 
                 (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx df-slice-idx-name-GLOBAL))
                      (dx-slice-idx-synt (datum->syntax stx dx-slice-idx-name-GLOBAL))
                      (df-full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (dx-full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'(* df-slice-range dx-slice-range))
                      (r-value-exp 
                        (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                         (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                           (list 'discard
                                 (list 'array-ref #,#'df-name-symb #,#'df-full-idx)
                                 (list 'array-ref #,#'dx-name-symb #,#'dx-full-idx))))))))
                (df-slice-colon_
                 (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp
                        (local-expand-opaque #'der-value)))
                     (define r-value-type (syntax-property #'r-value-exp 'landau-type))

                   ;(displayln (format "full-idx ~a\nlvalue-outer-prod-range ~a\nr-value-exp ~a\nr-value-type ~a" #'full-idx #'lvalue-outer-prod-range #'r-value-exp #'r-value-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                         (list 'discard
                               (list 'array-ref #,#'df-name-symb #,#'full-idx)
                               (list 'array-ref #,#'dx-name-symb 0)))))))
                (dx-slice-colon_
                 (with-syntax*
                     ((dx-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'dx-slice-range)
                      (r-value-exp
                        (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                         (list 'discard
                               (list 'array-ref #,#'df-name-symb 0)
                               (list 'array-ref #,#'dx-name-symb #,#'full-idx)))))))
                (else 
                 (begin
                   (throw-if-r-value-is-slice stx #'der-value)
                   (quasisyntax/loc stx
                     (list 'discard
                           (list 'array-ref #,#'df-name-symb 0)
                           (list 'array-ref #,#'dx-name-symb 0)))))))))))]))

(define-for-syntax (get-expanded-slice-index-start index-start)
  (let* ((index-start-expanded (if (syntax->datum index-start)
                                   (local-expand-opaque index-start)
                                   (is-type_ 'int-index #'0)))
         (index-start-expanded_ (if (atom-number index-start-expanded)
                                    (atom-number index-start-expanded)
                                    index-start-expanded)))
    (values index-start-expanded index-start-expanded_)))

(define-for-syntax (get-expanded-slice-index-end stx index-end array-range-stx)
  (let* ((index-end-expanded (if (syntax->datum index-end)
                                 (begin
                                   (local-expand-opaque index-end))
                                 (is-type_ 'int-index (datum->syntax stx array-range-stx))))
         (index-end-expanded_ (if (atom-number index-end-expanded)
                                  (atom-number index-end-expanded)
                                  (begin
                                    ;  (println (format "index-end-expanded ~a" index-end-expanded))
                                    index-end-expanded))))
    (values index-end-expanded index-end-expanded_)))


(define-for-syntax (slice-helper stx dx-index-start dx-index-end expanded-range-dx)
  (let*-values
    (((dx-index-start-stx dx-index-start-number) (get-expanded-slice-index-start dx-index-start))
     ((dx-index-end-stx dx-index-end-number) (get-expanded-slice-index-end stx dx-index-end expanded-range-dx)))

    (throw-if-not-type 'int-index 
                       dx-index-start-stx 
                       "Indexes expected to be expressions with 'int constants and/or loop variables")
    (throw-if-not-type 'int-index 
                       dx-index-end-stx 
                       "Indexes expected to be expressions with 'int constants and/or loop variables")
    #| (define dx-index-end-number-atom (atom-number dx-index-end-number)) |#
    #| (define dx-index-start-number-atom (atom-number dx-index-start-number)) |#
    
    (with-syntax*
      ((slice-range (quasisyntax/loc stx 
                                       (fx- #,dx-index-end-number
                                            #,dx-index-start-number))))

      (when (exact-integer? dx-index-start-number)
        (when (fx< dx-index-start-number 0)
          (raise-syntax-error
            #f
            (format "index must be nonnegative") stx)))
      (when (exact-integer? dx-index-end-number)
        (when (fx> dx-index-end-number expanded-range-dx)
          (raise-syntax-error
            #f
            (format "slice end index ~a is out of range [0, ~a)" dx-index-end-number expanded-range-dx) stx)))
      ; (println "foo")
      (with-syntax ((index-range-checks 
                      (quasisyntax/loc
                        stx
                        (begin
                          ;; NOTE: if index is not expanded to a number in compile time, insert check in runtime
                          #,(when (not (exact-integer? dx-index-start-number))
                              #`(when (fx< #,dx-index-start-number 0)
                                  (raise-syntax-error
                                    #f
                                    (format "index must be nonnegative") #'#,stx)))
                          #,(when (not (exact-integer? dx-index-end-number))
                              #`(when (fx> #,dx-index-end-number #,expanded-range-dx)
                                  (raise-syntax-error
                                    #f
                                    (format "slice end index ~a is out of range [0, ~a)" #,dx-index-end-number #,expanded-range-dx) #'#,stx)))))))
        (values dx-index-start-stx #'slice-range #'index-range-checks)))))

(define-for-syntax (throw-if-r-value-is-slice stx r-value-stx)
  (with-syntax*
      ((r-value-exp
         (local-expand-opaque r-value-stx))
       (r-value-type (syntax-property #'r-value-exp 'landau-type)))
    (when (is-slice? (syntax->datum #'r-value-type))
      (raise-syntax-error #f
                          (format "the right-hand side is slice, but the left-hand one is not") stx))))



;; TODO: refactor to reduce bindings duplication
(define-syntax (der-annot stx)
  (syntax-parse stx
    [(_ ({~literal get-value} df-name:id
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
     (let* ((name-1-dat (symbol->string (syntax->datum #'df-name)))
            (name-2-dat (symbol->string (syntax->datum #'dx-name)))
            (df-range (check-real-var-existence #'df-name 'get-range))
            (dx-range (check-real-arg-or-parameter-existence #'dx-name 'get-range))
            (df-type (search-var #'df-name))
            (dx-type (search-arg-or-parameter #'dx-name))
            (df-getter (make-getter-info (if (attribute df-idx) #'df-idx #'#f)
                                         (if (attribute df-slice-colon) #'df-slice-colon #'#f)
                                         (to-landau-type stx df-type)))
            (dx-getter (make-getter-info (if (attribute dx-idx) #'dx-idx #'#f)
                                         (if (attribute dx-slice-colon) #'dx-slice-colon #'#f)
                                         (to-landau-type stx dx-type)))
            (df-idx-expanded (if (attribute df-idx)
                               (local-expand-opaque #'df-idx)
                               #f))
            (dx-idx-expanded (if (attribute dx-idx)
                               (local-expand-opaque #'dx-idx)
                               #f))
            (expanded-range-df (atom-number #`#,(expand-range 
                                                  (if (attribute df-idx) 
                                                    #'df-idx 
                                                    #'#f) 
                                                  (if (equal? df-range '()) 
                                                    #f 
                                                    df-range))))
            (expanded-range-dx (atom-number #`#,(expand-range 
                                                  (if (attribute dx-idx) 
                                                    #'dx-idx 
                                                    #'#f) 
                                                  (if (equal? dx-range '()) 
                                                    #f 
                                                    dx-range))))
            (df-slice-colon_ (attribute df-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon)))

       (check-proper-getter (if (attribute df-idx) #'df-idx #'#f) 
                            (attribute df-slice-colon) 
                            expanded-range-df 
                            df-idx-expanded 
                            name-1-dat stx)
       (check-proper-getter (if (attribute dx-idx) #'dx-idx #'#f) 
                            (attribute dx-slice-colon) 
                            expanded-range-dx 
                            dx-idx-expanded 
                            name-2-dat stx)
       (let*-values (((v0 v1 v2 df-name-src-pos) (search-name-backrun stx #'df-name))
                    ;  ((v3 v4 v5 dx-name-src-pos) (search-name-backrun stx #'dx-name))
                     ((dx-index-start-stx dx-slice-range dx-index-range-checks)
                      (if dx-slice-colon_ 
                          (slice-helper stx #'dx-index-start #'dx-index-end expanded-range-dx)
                          (values #f #f #f)))
                     ((df-index-start-stx df-slice-range df-index-range-checks)
                      (if df-slice-colon_ 
                          (slice-helper stx #'df-index-start #'df-index-end expanded-range-df)
                          (values #f #f #f))))
         (with-syntax* ((df-name-symb (var-symbol name-1-dat df-name-src-pos))
                        ;; NOTE: Do not need src-pos for dx
                        (dx-name-symb name-2-dat)
                        (df-slice-range df-slice-range)
                        (dx-slice-range dx-slice-range)
                        (dx-index-range-checks dx-index-range-checks)
                        (df-index-range-checks df-index-range-checks))
           (hash-set! (syntax-parameter-value #'dx-names-set) (symbol->string (syntax->datum #'dx-name)) #t)
           (match (list (getter-info-type df-getter)
                        (getter-info-type dx-getter) 
                        "<- der-value")
             ((list 'var 'var "<- der-value")
              (begin
               (throw-if-r-value-is-slice stx #'der-value)
               (quasisyntax/loc stx
                                (begin
                                 der-value
                                 (list 'der-annot
                                       (list 'array-ref #,#'df-name-symb 0)
                                       (list 'array-ref #,#'dx-name-symb 0))))))

             ((list 'var 'cell "<- der-value")
              (begin
               (throw-if-r-value-is-slice stx #'der-value)
               (quasisyntax/loc stx
                        (begin
                          der-value
                          #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                          (list 'der-annot
                                (list 'array-ref #,#'df-name-symb 0)
                                (list 'array-ref #,#'dx-name-symb #,dx-idx-expanded))))))

             ((list 'var 'slice "<- der-value")
              (with-syntax*
                     ((dx-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'dx-slice-range)
                      (r-value-exp (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (let ((#,#'slice-idx-synt 0))
                         der-value
                         (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                           (let ((#,#'slice-idx-synt #,#'dx-slice-idx-synt))
                             (list 'der-annot
                                   (list 'array-ref #,#'df-name-symb 0)
                                   (list 'array-ref #,#'dx-name-symb #,#'full-idx)))))))))

             ((list 'cell 'var "<- der-value")
              (begin
                (throw-if-r-value-is-slice stx #'der-value)
                (quasisyntax/loc stx
                  (begin
                    der-value
                    #,(range-check-stx stx df-idx-expanded expanded-range-df)
                    (list 'der-annot
                          (list 'array-ref #,#'df-name-symb #,df-idx-expanded)
                          (list 'array-ref #,#'dx-name-symb 0))))))

             ((list 'cell 'cell "<- der-value")
              (begin
               (throw-if-r-value-is-slice stx #'der-value)
               (quasisyntax/loc 
                  stx
                (begin
                  der-value
                  #,(range-check-stx stx df-idx-expanded expanded-range-df)
                  #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                  (when (and (equal? df-name-symb dx-name-symb) (equal? #,df-idx-expanded #,dx-idx-expanded))
                    (raise-syntax-error
                     #f
                     (format "Can not override annotation: ~a[~a] ' ~a[~a] == 1.0" df-name-symb #,df-idx-expanded dx-name-symb #,dx-idx-expanded) #'#,stx))
                  (list 'der-annot
                        (list 'array-ref #,#'df-name-symb #,df-idx-expanded)
                        (list 'array-ref #,#'dx-name-symb #,dx-idx-expanded))))))

             ((list 'cell 'slice "<- der-value")
              (with-syntax*
                  ((df-idx-expanded-stx df-idx-expanded) 
                   (dx-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                   (full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                   (lvalue-outer-prod-range #'dx-slice-range)
                   (r-value-exp (local-expand-opaque #'der-value)))
                (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                (quasisyntax/loc
                    stx
                  (begin
                    #,#'dx-index-range-checks
                    #,(range-check-stx stx df-idx-expanded expanded-range-df)
                    #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                    (let ((#,#'slice-idx-synt 0))
                         der-value
                         (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                           (list 'der-annot
                                 (list 'array-ref #,#'df-name-symb #,#'df-idx-expanded-stx)
                                 (list 'array-ref #,#'dx-name-symb #,#'full-idx))))))))

             ((list 'slice 'var "<- der-value")
              (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       (let ((#,#'slice-idx-synt 0))
                            der-value
                            (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                              (list 'der-annot
                                    (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                    (list 'array-ref #,#'dx-name-symb 0))))))))

             ((list 'slice 'cell "<- der-value")
              (with-syntax*
                     ((dx-idx-expanded-stx dx-idx-expanded) 
                      (df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp (local-expand-opaque #'der-value)))
                   (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       ;; assign placeholder value to slice variable which result is ignored. 
                       (let ((#,#'slice-idx-synt 0))
                         der-value
                         (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                           (list 'der-annot
                                 (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                 (list 'array-ref #,#'dx-name-symb #,#'dx-idx-expanded-stx))))))))

             ((list 'slice 'slice "<- der-value")
              (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx df-slice-idx-name-GLOBAL))
                      (dx-slice-idx-synt (datum->syntax stx dx-slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (df-full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (dx-full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'(fx* df-slice-range dx-slice-range))
                      (r-value-exp (local-expand-opaque #'der-value)))
                     (define r-value-type (syntax-property #'r-value-exp 'landau-type))
                     #| (displayln (format "r-value-type ~a" r-value-type)) |#
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx r-value-type #'lvalue-outer-prod-range)
                       ;; assign placeholder value to slice variable which result is ignored. 
                       (let ((#,#'slice-idx-synt 0))
                              der-value ;; to trigger range checks
                              (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                  (list 'der-annot
                                        (list 'array-ref #,#'df-name-symb #,#'df-full-idx)
                                        (list 'array-ref #,#'dx-name-symb #,#'dx-full-idx))))))))))
           )))]))

(define-syntax (single-term stx)
  (syntax-parse stx
                [(_ body)
                 (syntax/loc stx body)]))

;; FIXME Add new context for variables
(define-syntax (expr-body stx)
  (syntax-parse stx
    [(_ body ...)
     (syntax/loc stx
                 (list body ...))]))



(define-syntax (func-body stx)
  (syntax-parse stx
    [(_ body ...)
     (syntax/loc stx
       (list body ...))]))

(define-for-syntax (binary-op-cast op1 op2 stx)
  (let ((type1 (syntax-property op1 'landau-type))
        (type2 (syntax-property op2 'landau-type))
        (op1-atom (atom-number op1))
        (op2-atom (atom-number op2)))
    ; (println (format "~a ~a" type1 type2))
    (match (list type1 type2)
      ((list 'int 'int)
       (values op1 op2 'int))
      ((list 'real 'real)
       (values op1 op2 'real))
      ((list 'int 'real)
       (values (datum->syntax
                op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
               op2 'real))
      ((list 'real 'int)
       (values op1 (datum->syntax
                    op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
               'real))
      ((list 'int-index 'int-index)
       (values op1 op2 'int-index))
      ((list 'int 'int-index)
       (values op1 op2 'int))
      ((list 'int-index 'int)
       (values op1 op2 'int))
      ((list 'real 'int-index)
       (values op1 (datum->syntax
                    op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
               'real))
      ((list 'int-index 'real)
       (values (datum->syntax
                op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
               op2 'real))
      (`(,(? is-slice? type1) ,(? is-slice? type2))
       (let ((range1 (get-slice-range type1))
             (range2 (get-slice-range type2))
             (slice-type1 (get-slice-type type1))
             (slice-type2 (get-slice-type type2)))
         (match (list slice-type1 slice-type2)
           ((list 'real 'real)
            (values
             #`(if (equal? #,range1 #,range2)
                   #,op1 
                   (raise-syntax-error #f (format "cannot cast ranges ~v and ~v" #,range1 #,range2) #'#,stx))
             op2 
             type1))
           (else
            (raise-syntax-error #f (format "cannot1 cast ~v and ~v" slice-type1 slice-type2) stx)))))
      (`(,(? is-slice? type1) ,(? (compose not is-slice?) type2))
       (let ((range1 (get-slice-range type1))
             (slice-type1 (get-slice-type type1)))
         (match (list slice-type1 type2)
           ((list 'real 'int)
            (values 
             op1 
             ;; TODO: check it; compare with 
             (datum->syntax op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
             type1))
           ((list 'real 'int-index)
            (values
             op1
             (datum->syntax op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
             type1))
           ((list 'real 'real) (values op1 op2 type1))
           (else
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type1 type2) stx)))))
      (`(,(? (compose not is-slice?) type1) ,(? is-slice? type2))
       (let ((range2 (get-slice-range type2))
             (slice-type2 (get-slice-type type2)))
         (match (list slice-type2 type1)
           ((list 'real 'int)
            (values 
             (datum->syntax op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
             op2
             type2))
           ((list 'real 'int-index)
            (values
             (datum->syntax op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
             op2
             type2))
           ((list 'real 'real) (values op1 op2 type2))
           (else
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type2 type1) stx)))))
      (else
       (raise-syntax-error #f (format "cannot cast ~v and ~v. op1: ~v, op2: ~v." type1 type2 op1 op2) stx)))))


(define-syntax (expr stx)
  (syntax-parse stx
    ((_ expr1 op:plus-or-minus expr2)
     (let* ((expr1-expanded (local-expand #'expr1 'expression (list)))
            (expr2-expanded (local-expand #'expr2 'expression (list)))
            (expr1-is-int-index (equal? (syntax-property expr1-expanded 'landau-type) 'int-index))
            (expr2-is-int-index (equal? (syntax-property expr2-expanded 'landau-type) 'int-index))
            (op (syntax->datum #'op)))
       (let*-values (((expr1-casted expr2-casted type) (binary-op-cast expr1-expanded expr2-expanded stx))
                     ((expr1-atom) (atom-number expr1-casted))
                     ((expr2-atom) (atom-number expr2-casted)))
         (cond
           ((equal? type 'real)
            (is-real
             (if (and expr1-atom expr2-atom)
                 (quasisyntax/loc stx #,((if (equal? op "+") fl+ fl-) expr1-atom expr2-atom))
                 (with-syntax* ((e1 (if (and (equal? #f expr1-atom) (not expr1-is-int-index))
                                        expr1-casted
                                        #'(list)))
                                (e2 (if (and (equal? #f expr2-atom) (not expr2-is-int-index))
                                        expr2-casted
                                        #'(list))))
                   (syntax/loc stx (append e1 e2))))))
           ((is-slice? type)
            (is-type_ type (with-syntax* ((e1 (if (and (equal? #f expr1-atom) (not expr1-is-int-index))
                                                  expr1-casted
                                                  #'(list)))
                                          (e2 (if (and (equal? #f expr2-atom) (not expr2-is-int-index))
                                                  expr2-casted
                                                  #'(list))))
                             (syntax/loc stx (append e1 e2)))))
           ((equal? type 'int)
            (is-int
             (syntax/loc stx '())))
           ((equal? type 'int-index)
            (is-int-index
               (if (and expr1-atom expr2-atom)
                 (quasisyntax/loc stx #,((if (equal? op "+") fx+ fx-) expr1-atom expr2-atom))
                 (quasisyntax/loc stx (#,(if (equal? op "+") #'fx+ #'fx-) #,expr1-expanded #,expr2-expanded)))))
           (else
            (raise-syntax-error #f (format "unsupported type: ~v" type) stx))))))
    
    [(_ "(" expr1 ")")
     (syntax/loc stx
       expr1)]
    
    [(_ expr1)
     (with-syntax ((expanded (local-expand #'expr1 'expression (list))))
       (syntax-track-origin (syntax/loc stx expanded) #'expr1 #'expr))]))

(define/contract-for-syntax
  (arity-check func-name n-expected n-provided stx)
  (-> string? integer? integer? (syntax/c any/c) void?)
  (unless (equal? n-expected n-provided)
    (raise-syntax-error
      #f
      (format "function ~a expects ~a parameters, but ~a provided" func-name n-expected n-provided)
      stx)))
   
(define/contract-for-syntax
  (ref->var-symbol ref)
  (-> backrun-ref/c var-symbol/c)
  (match ref
    ((list _ vs _) vs)))

(define/contract-for-syntax
  (bind-parameters-to-arguments inlined-function-name actions-list bindings)
  (-> string? (listof any/c) (hash/c var-symbol/c (or/c func-arg-binding/c (symbols 'constant)))
      (listof any/c))
  (define flat-actions-list (reverse (splice-nested actions-list)))
  (define msg-template "bug: not self-sufficent function should not contain ~a term")
  (define (rebind-ref _inlined-function-name ref bindings)
    (define var-name (match ref ((list _ (var-symbol name _) _) name)))
    (define (is-constant? _name)
      (hash-has-key? constants (string->symbol _name)))
    (match ref
      ((list 'array-ref vs idx) 
       ;; NOTE: If there is no such a key, then the key is a 
       ; local variable and should not be changed.
       ;; NOTE: Need to either preserve or override the idx depending on the arg type:
       ;; var -> preserve (really does not matter, because it is 0 in both cases)
       ;; cell -> override
       ;; array -> preserve 
       (match (hash-ref bindings vs (thunk 'not-a-function-parameter))

         ((func-arg-binding 'var (list 'array-ref arg-vs arg-idx))
          (list 'array-ref arg-vs idx))

         ((func-arg-binding 'cell (list 'array-ref arg-vs arg-idx))
          (list 'array-ref arg-vs arg-idx))

         ((func-arg-binding 'array (list 'array-ref arg-vs arg-idx))
          (list 'array-ref arg-vs idx))

         ('not-a-function-parameter
          (cond
            ((is-constant? var-name) (list 'array-ref vs idx))

            ;; NOTE: inlined function's local variables are prepanded with _
            (else (match ref
                    ((list 'array-ref (var-symbol name src-pos) ref-idx)
                     (list 'array-ref (var-symbol (make-inlined-variable-name _inlined-function-name name) src-pos) ref-idx))))))

         ('constant 'constant)))))

  (define (rebind-refs-list _inlined-function-name refs-list bindings)
    (filter (lambda (x) (not (equal? x 'constant))) 
            (for/list ((ref (in-list refs-list)))
              (rebind-ref _inlined-function-name ref bindings))))

  (for/list ((action (in-list flat-actions-list)))
    (match action
      ((list 'func-assign func-ref refs-list)
       (list 'assign 
             (rebind-ref inlined-function-name func-ref bindings)
             (rebind-refs-list inlined-function-name refs-list bindings)))

      ((list 'assign l-val-ref refs-list)
       (list 'assign (rebind-ref inlined-function-name l-val-ref bindings)
             (rebind-refs-list inlined-function-name refs-list bindings)))

      ((list 'var-decl vs type)
       (list 'var-decl 
             ;; NOTE: inlined function's local variables are prepanded with _
             (match vs 
               ((var-symbol name src-pos)
                (var-symbol (make-inlined-variable-name inlined-function-name name) src-pos)))
             type))

      ((list 'der-annot ref-1 ref-2)
       (error (format msg-template 'der-annot)))

      ((list 'der-apply df dx)
       (error (format msg-template 'der-apply)))

      ((list 'discard df dx)
       (error (format msg-template 'discard))))))


(define/contract-for-syntax
  (add-function-return-variable-binding! bindings 
                                         func-name-fake-vs 
                                         func-name-true-vs 
                                         func-getter)
  (-> (hash/c var-symbol/c (or/c func-arg-binding/c (symbols 'constant)))
      var-symbol/c
      var-symbol/c
      getter-info/c
      void?)
  (hash-set! bindings 
             func-name-fake-vs 
             (func-arg-binding (getter-info-type func-getter) (list 'array-ref 
                                                                    func-name-true-vs
                                                                    0))))


(define-syntax (func-call stx)
  (syntax-parse stx
    ((_ function-name "(" ({~literal parlist} par*:par-spec ...) ")")
     (let* ((func-str (symbol->string (syntax->datum #'function-name)))
            (fake-src-pos 0)
            ;; NOTE: funcs-info-GLOBAL use var-symbol/c as keys.
            ; For functions src-pos is always 0. It does not make sense
            ; and should be changed later.
            (func-name-fake-vs (var-symbol func-str fake-src-pos))
            ;; NOTE: True src-pos is used to generate distinct
            ; return variables for each func-call (inlining).
            (true-src-pos (syntax-position #'function-name))
            ;; NOTE: return var name is modified to handle ambiguity.
            ; In semantics.rkt name is modified in the same way
            (function-return-variable-name (format "~a~a" func-str true-src-pos))
            (func-name-true-vs (var-symbol function-return-variable-name true-src-pos))
            (func-pars-list (for/list ((par-value (in-list (syntax-e #'(par*.par-value ...)))))
                              par-value))
            (par-list-len (length func-pars-list))
            (builtin-functions (hash-keys BUILT-IN-FUNCTIONS)))
       ;; FIXME: Check func name, type, casting conditions
       (if (member func-str builtin-functions)
         ;; NOTE: Built-in functions
         (let ((expected-arity (length (hash-ref BUILT-IN-FUNCTIONS func-str))))
           (if (equal? func-str "pow")
             (begin
               (arity-check func-str expected-arity par-list-len stx)
               (with-syntax* ((expr1 (car func-pars-list))
                              (expr2 (cadr func-pars-list))
                              (expanded1 (local-expand #'expr1 'expression (list)))
                              (expanded2 (local-expand #'expr2 'expression (list))))
                             (with-syntax-property 'getter-info (getter-info 'var #f)
                                                   (is-type_ 'real 
                                                             (syntax/loc stx (append expanded1 expanded2))))))
             (begin
               (arity-check func-str expected-arity par-list-len stx)
               (with-syntax* ((expr1 (car func-pars-list))
                              (expanded1 (local-expand #'expr1 'expression (list))))
                             (with-syntax-property 'getter-info (getter-info 'var #f)
                              (is-type_ 'real 
                                        (syntax/loc stx expanded1)))))))
         ;; NOTE: User-defined functions
         ;; FIXME perform checks
         (let*
           ((func-info (if (hash-has-key? funcs-info-GLOBAL func-name-fake-vs)
                         (hash-ref funcs-info-GLOBAL func-name-fake-vs)
                         (raise-syntax-error #f (format "function ~a is not defined." func-str) stx)))
            (type (if (func-info-output-range func-info)
                    (landau-type (func-info-output-base-type func-info)
                                 (func-info-output-range func-info))
                    (landau-type (func-info-output-base-type func-info))))
            (func-getter (if (func-info-output-range func-info)
                           (getter-info 'array #f)
                           (getter-info 'var #f)))
            (subfunc-call-box-value (make-state (list)))
            (func-call-box-value (syntax-parameter-value #'func-call-box)))

           (with-syntax
             ((slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
              (pars-list-expanded 
                (for/list ((par (in-list func-pars-list)))
                  (with-syntax* 
                    ((par par))
                    (extract 
                      (local-expand
                        #`(syntax-parameterize
                            ((func-call-box #,subfunc-call-box-value)
                             (func-is-called-inside-argument #t))
                            par) 'expression (list)))))))

             #| (displayln "FIXME: check if funciton local var are in der-table if they need a derivative after inlining") |#

             ;; NOTE: inline funciton body and bind it's parameters to the passed arguments
             #| (displayln "FIXME: generate normalizied arg variables in semantics.rkt") |#
             (define current-func-arg-normalizations (make-state (list)))
             (define module-funcs-table (syntax-parameter-value #'module-funcs-table-parameter))
             (define args 
               (for/list ((arg-number (in-naturals))
                          (arg (syntax-e #'pars-list-expanded)))
                 (define evaluated-arg (evaluate arg))
                 (match evaluated-arg
                   ((list) 'constant)

                   ((list (list ref vs n)) 
                    (begin
                      (define _getter (syntax-property arg 'getter-info))
                      (unless _getter 
                        (error (format "bug: arg: ~a has no getter property" arg)))
                      (define getter (getter-info-type _getter))
                      (func-arg-binding getter (list ref vs n))))

                   ((list refs ...) 
                    (begin
                      ;; NOTE semantics.rkt can't handle expressions inside an argument
                      ;; TODO handle them
                      (raise-syntax-error 
                        #f 
                        "expressions are not allowed inside function arguments" stx)
                      (define args-norms (read-state current-func-arg-normalizations))
                      (define norm-arg-vs (var-symbol
                                            (make-norm-arg-name func-str arg-number)
                                            (syntax-position arg)))
                      (define norm-arg-ref (list 'array-ref norm-arg-vs 0))
                      (define arg-normalization
                        (list 'assign 
                              norm-arg-ref
                              refs))
                      (write-state! current-func-arg-normalizations
                                    (append args-norms
                                            (list arg-normalization)))
                      (define arg-getter-type 'var)
                      (func-arg-binding arg-getter-type norm-arg-ref))))))

             (define args-types 
               (for/list ((arg (syntax-e #'pars-list-expanded)))
                 (define _arg-getter (syntax-property arg 'getter-info))
                 (define _arg-name (syntax-property arg 'get-value-name))
                 (define _arg-type (to-landau-type stx (syntax-property arg 'landau-type)))
                 (match _arg-getter
                   (#f (list 'constant (to-landau-type stx _arg-type)))
                   (_ (begin
                        (define _arg-full-type
                          (syntax-property arg 'full-type))
                        (list _arg-getter _arg-full-type))))))

             (define func-template (hash-ref module-funcs-table (string->symbol func-str)))
             (define func-template-parameters-vs 
               (for/list ((par (in-list (function-inline-template-.parameters func-template))))
                 (match par
                   ((list par-vs _) par-vs))))

             #| (displayln (format "args-types: ~a" args-types)) |#
             #| (displayln (format "pars-types: ~a" (function-inline-template-.parameters func-template))) |#

             (define (landau-type->string t)
               (match t
                 ((list b (list)) (format "~a" b))
                 ((list b (list n)) (format "~a[~a]" b n))))

             (define (check-arguments _stx _func-str expected-types given-types)
               (define expected-n (length expected-types))
               (define given-n (length given-types))
               (unless (equal? expected-n given-n)
                 (raise-syntax-error 
                   #f 
                   (format "function ~a arity mismatch: expected ~a, given: ~a"
                           _func-str expected-n given-n) _stx))
               (for ((arg-n (in-naturals))
                     (expected-type-info (in-list expected-types))
                     (given-type-info (in-list given-types)))
                 (define expected-type
                   (match expected-type-info
                     ((list _ type) type)))
                 (define-values (given-type given-getter)
                   (match given-type-info
                     ((list getter type) (values type getter))))
                 (when (equal? given-getter 'slice)
                   (raise-syntax-error 
                     #f 
                     (format "function ~a argument type mismatch: argument slices are not implemented"
                             _func-str) _stx))
                 (unless (equal? expected-type given-type)
                   (raise-syntax-error 
                     #f 
                     (format 
                       "function ~a argument type mismatch: argument ~a expected to be ~a, given: ~a"
                       _func-str arg-n (landau-type->string expected-type) (landau-type->string given-type))
                     _stx))))

             (check-arguments stx func-str (function-inline-template-.parameters func-template) args-types)
             (define/contract bindings 
                              (hash/c var-symbol/c (or/c func-arg-binding/c (symbols 'constant)))
                              (make-hash (zip cons 
                                              func-template-parameters-vs
                                              args)))
             (add-function-return-variable-binding! bindings 
                                                    func-name-fake-vs 
                                                    func-name-true-vs 
                                                    func-getter)

             ;; NOTE: function's return variable (unique to an each call) declaration
             (define current-function-return-var-decl
               (list 'var-decl func-name-true-vs (to-landau-type stx type)))

             (define inlined-function 
               (bind-parameters-to-arguments 
                 func-str 
                 (function-inline-template-.actions-list func-template) 
                 bindings)) 

             (unless (equal? func-call-box-value 'func-call-box-not-set)
               (begin
                 (define updated-inlinded-func-actions-list
                   (append 
                     ;; NOTE: func-call list from children   
                     (read-state subfunc-call-box-value) 
                     ;; NOTE: func-call list from neighbors  
                     ;; on the same call depth level 
                     (read-state func-call-box-value) 
                     ;; NOTE: if argument includes multiple variables
                     ; then it is assigned to the new variable and that variable
                     ; substitutes the argument.
                     (read-state current-func-arg-normalizations)
                     (list current-function-return-var-decl)
                     inlined-function))
                 #| (displayln (format "updated-inlinded-func-actions-list: ~a" updated-inlinded-func-actions-list)) |#
                 (write-state! func-call-box-value updated-inlinded-func-actions-list))) 

             #| (displayln (format "~a: func-call-box-value:" func-str)) |#
             #| (pretty-print func-call-box-value) |#

             (with-syntax-property
               'full-type (to-landau-type stx type)
               (with-syntax-property 
                 'get-value-name (string->symbol function-return-variable-name) 
                 (with-syntax-property 
                   'getter-info func-getter
                   (if (syntax-parameter-value #'func-is-called-inside-argument)
                     ;; NOTE: Function called inside the argument of another function.
                     ; Return a link to the function's return-variable.
                     ; Actions of assignation to the variable are appended to the
                     ; list of actions and written to `func-call-box-value`.
                     ; The list of actions is then fetched during the `assignation` macro expansion
                     ; and prepanded to the assignation action. 
                     (match (getter-info-type func-getter)
                       ('array
                        (with-syntax-property 
                          'landau-type (to-landau-type stx type)
                          (quasisyntax/loc 
                            stx
                            (list (list 'array-ref #,func-name-true-vs 0)))))
                       (_
                         (with-syntax-property 
                           'landau-type (func-info-output-base-type func-info)
                           (quasisyntax/loc 
                             stx
                             (list (list 'array-ref #,func-name-true-vs 0))))))
                     ;; NOTE: This function call is top-level and in case of
                     ; function returns an array, return-value is inserted
                     ; into the loop, where the output array is copied to the
                     ; assignation left-hand side variable.
                     (match (getter-info-type func-getter)
                       ('array 
                        (with-syntax-property
                          'landau-type (to-landau-type stx type)
                          (quasisyntax/loc
                            stx
                            (list (list 'array-ref #,func-name-true-vs #,#'slice-idx-synt)))))
                       (_ (with-syntax-property 
                            'landau-type (func-info-output-base-type func-info)
                            (quasisyntax/loc
                              stx
                              (list (list 'array-ref #,func-name-true-vs 0)))))))))))))))))


(define-syntax (term stx)
  (syntax-parse stx
    ((_ expr1 op:mul-or-div expr2)
     (let* ((expr1-expanded (local-expand #'expr1 'expression (list)))
            (expr2-expanded (local-expand #'expr2 'expression (list)))
            (expr1-is-int-index (equal? (syntax-property expr1-expanded 'landau-type) 'int-index))
            (expr2-is-int-index (equal? (syntax-property expr2-expanded 'landau-type) 'int-index))
            (op (syntax->datum #'op)))
       (let*-values (((expr1-casted expr2-casted type) (binary-op-cast expr1-expanded expr2-expanded stx))
                     ((expr1-atom) (atom-number expr1-casted))
                     ((expr2-atom) (atom-number expr2-casted)))
         (cond
           ((equal? type 'real)
            (is-real
             (if (and expr1-atom expr2-atom)
                 (quasisyntax/loc stx #,((if (equal? op "*") fl* fl/) expr1-atom expr2-atom))
                 (with-syntax* ((e1 (if (and (equal? #f expr1-atom) (not expr1-is-int-index))
                                        expr1-casted
                                        #'(list)))
                                (e2 (if (and (equal? #f expr2-atom) #t (not expr2-is-int-index))
                                        expr2-casted
                                        #'(list))))
                   (syntax/loc stx (append e1 e2))))))
           ((is-slice? type)
            (is-type_ type (with-syntax* ((e1 (if (and (equal? #f expr1-atom) (not expr1-is-int-index))
                                                  expr1-casted
                                                  #'(list)))
                                          (e2 (if (and (equal? #f expr2-atom) #t (not expr2-is-int-index))
                                                  expr2-casted
                                                  #'(list))))
                             (syntax/loc stx (append e1 e2)))))
           ((equal? type 'int)
            (is-int
             (syntax/loc stx '())))
           ((equal? type 'int-index)
            (is-int-index
             (if (and expr1-atom expr2-atom)
                 (quasisyntax/loc stx #,((if (equal? op "*") fx* fxquotient) expr1-atom expr2-atom))
                 (quasisyntax/loc stx (#,(if (equal? op "*") #'fx* #'fxquotient) #,expr1-expanded #,expr2-expanded)))))
           (else
            (raise-syntax-error #f (format "unsupported type: ~v" type) stx))))))
       
    [({~literal term} term1)
     (with-syntax ((expanded (local-expand #'term1 'expression (list))))
       ;  (println #'expanded)
       ;  (println (format "term type: ~a" (syntax-property #'expanded 'landau-type)))
       ;  (println (format "term type emitted: ~a" (syntax-property (syntax-track-origin (syntax/loc stx expanded) #'term1 #'expr) 'landau-type)))
       ;  (println "")
       (syntax-track-origin (syntax/loc stx expanded) #'term1 #'expr))]))
    

(define-syntax (element stx)
  (syntax-parse stx
    [(_ number)
     (with-syntax* ((expanded (local-expand #'number 'expression (list)))
                    (atom_ (atom-number #'expanded))
                    (type (syntax-property #'expanded 'landau-type)))
       (let ((type_ (syntax->datum #'type))
             (atom (syntax->datum #'atom_)))
         (match type_
           ('int-index
            (syntax-track-origin (syntax/loc stx expanded) #'number #'expr))
           ('real
            (if atom
                (is-type_ 'real #'(list))
                (is-type_ 'real #'expanded)))
           ('int (is-type_ type_ #'(list)))
           (else
            (cond
              ((is-slice? type_)
               (syntax-track-origin (syntax/loc stx expanded) #'number #'expr))
              (else (error (format "bug: unknown type: ~a" type_))))))))]
    [(_ "(" expr-stx ")")
     (syntax/loc stx
       expr-stx)]
    ))

(define-syntax (factor stx)
  (syntax-parse 
   stx
   [({~literal factor} primary "^" factor)
    (error "^ operator is not implemented")]

   [({~literal factor} primary)
    (with-syntax* ((expanded (local-expand #'primary 'expression (list)))
                   (atom_ (atom-number #'expanded))
                   (type (syntax-property #'expanded 'landau-type)))
      (let ((type_ (syntax->datum #'type))
            (atom (syntax->datum #'atom_)))
        (match type_
          ('int-index ;; FIXME: what is #'expr?
           (syntax-track-origin (syntax/loc stx expanded) #'primary #'expr))
          ('real
           (if atom
             (is-type_ 'real #'(list))
             (is-type_ 'real #'expanded)))
          ('int (is-type_ type_ #'(list)))
          (else
           (cond
             ((is-slice? type_)  ;; FIXME: what is #'expr?
              (syntax-track-origin (syntax/loc stx expanded) #'primary #'expr))
             (else (error (format "bug: unknown type: ~a" type_))))))))]))

(define-syntax (primary stx)
  (syntax-parse stx
                [(_ ((~literal unop) p-or-m) primary)
                 (let ((op (syntax->datum #'p-or-m)))
                  (with-syntax*
                   ((expanded (local-expand #'primary 'expression (list))))
                   (let ((expanded-basic-type (syntax-property #'expanded 'landau-type))
                         (atom (atom-number #'expanded))) 
                    (cond
                      ((equal? expanded-basic-type 'int-index)
                      (is-type_ expanded-basic-type
                         (if atom
                            (quasisyntax/loc stx
                                  #,((if (equal? op "+") fx+ fx-) 0 atom))
                            (datum->syntax stx
                                  `(,(if (equal? op "+") #'fx+ #'fx-) 0 ,#'expanded)))))
                     (else (syntax/loc stx
                                       expanded))))))]

                [(_ primary)
                 (with-syntax* ((expanded (local-expand #'primary 'expression (list)))
                                (atom_ (atom-number #'expanded))
                                (type (syntax-property #'expanded 'landau-type)))
                   (let ((type_ (syntax->datum #'type))
                         (atom (syntax->datum #'atom_)))
                     (match type_
                       ('int-index ;; FIXME: what is #'expr?
                        (syntax-track-origin (syntax/loc stx expanded) #'primary #'expr))
                       ('real
                        (if atom
                          (is-type_ 'real #'(list))
                          (is-type_ 'real #'expanded)))
                       ('int (is-type_ type_ #'(list)))
                       (else
                        (cond
                          ((is-slice? type_)  ;; FIXME: what is #'expr?
                           (syntax-track-origin (syntax/loc stx expanded) #'primary #'expr))
                          (else (error (format "bug: unknown type: ~a" type_))))))))]))
    

(define-syntax (number stx)
  (syntax-parse stx
    (((~literal number) number-stx)
     (syntax-property
      #'number-stx
      'landau-type
      (if (exact? (syntax->datum #'number-stx))
          'int-index
          'real)
      #t))))
      

;; TODO: add checks as in semantics (assignation to a const, handle assignation to int, etc)
;; FIXME: Make array-return funciton assignation logic as in semantics assignation 
(define-syntax (assignation stx)
  (syntax-parse stx

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
              
              (datum->syntax 
                stx
                `(assignation ,#'name ,@getter-for-splice "="
                              (expr
                                (expr 
                                  (term
                                    (factor 
                                      (primary 
                                        (element 
                                          (get-value ,#'name ,@getter-for-splice)))))) ,binop ,#'value))))
             ((or "*=" "/=")
              (datum->syntax 
                stx
                `(assignation ,#'name ,@getter-for-splice "="
                              (expr
                                (term
                                  (term 
                                    (factor 
                                      (primary 
                                        (element 
                                          (get-value ,#'name ,@getter-for-splice))))) ,binop ,#'value))))
              ))))

    (((~literal assignation) name:id 
                             getter:getter-cls
                             "=" value:expr)

     (let*-values (((name_) (syntax->datum #'name))
                   ((func-name) (syntax-parameter-value #'function-name))
                   ((func-basic-type) (car (syntax-parameter-value #'function-return-type)))
                   ((slice-colon_) (syntax->datum #'getter.slice-colon))
                   ((fake-src-pos) (values 0))
                   ((symbol type src-pos)
                    (cond
                      ((hash-has-key? constants name_)
                       (raise-syntax-error #f "assignation to constant is not allowed" #'name))
                      ((search-argument name_ (syntax-parameter-value #'current-arguments))
                       (raise-syntax-error #f "assignation to argument is not allowed" #'name))
                      ((hash-has-key? parameters name_) (raise-syntax-error #f "assignation to parameter is not allowed" #'name))
                      ((equal? name_ func-name)
                       (let ((func-return-value (syntax-parameter-value #'function-return-value)))
                         (values
                          (datum->syntax stx func-return-value)
                          (syntax-parameter-value #'function-return-type)
                          fake-src-pos)))
                      (else
                       (let ((var (search-variable name_ (syntax-parameter-value #'current-variables))))
                         (cond
                           (var
                            (values (datum->syntax stx (variable-symbol var))
                                    (variable-type var)
                                    (variable-src-pos var)))
                           (else
                            (raise-syntax-error #f "variable not found" #'name)))))))
                   ((maybe-array-range) (cadr type))
                   ((func-inline-list) (make-state (list))))
                   
       (with-syntax* ((index-exp (local-expand-opaque #'getter.index))
                      (index-start-expanded (if slice-colon_ 
                                              (if (syntax->datum #'getter.index-start)
                                                (local-expand-opaque #'getter.index-start)
                                                (is-type_ 'int-index #'0))
                                              #f))
                      (index-end-expanded (if slice-colon_
                                            (if (syntax->datum #'getter.index-end)
                                              (local-expand-opaque #'getter.index-end)
                                              (is-type_ 'int-index (datum->syntax stx (car maybe-array-range))))
                                            #f))
                      (value-exp_ (extract 
                                    (local-expand
                                      #`(syntax-parameterize
                                          ;; NOTE: func-call mutate list inside state. Resulted list 
                                          ; constains spliced actions list of inlined functions 
                                          ; called in the right-hand side of the current assignation
                                          ((func-call-box #,func-inline-list)) 
                                          value) 'expression (list))))
                      (value-exp #'value-exp_))

         (define value-type (syntax-property #'value-exp_ 'landau-type))
         (define-values 
           (value-base-type value-range) 
           (match value-type
                              ((list b (list n)) (values b n))
                              (b (values b #f))))

         #| (displayln (format "value-type ~a value-base-type ~a value-range ~a" value-type value-base-type value-range)) |#

         ;; NOTE: Indexes are resolved at zeroth pass
         (let* ((err-msg "index must be an integer expression of constants and loop variables")
                (throw-if-not-int-index
                 (lambda (idx-stx idx-expanded-stx)
                   (unless (equal? (syntax-property idx-expanded-stx 'landau-type) 'int-index)
                     (raise-syntax-error #f err-msg idx-stx)))))
           (when (syntax->datum #'getter.index)
             (throw-if-not-int-index #'getter.index #'index-exp))
           (when (and (not slice-colon_) (is-slice? value-type))
             (raise-syntax-error 
               #f 
               "the right-hand expression is a slice, but the left-hand one is not. A slice must be assigned to a slice"
               stx))
           (when slice-colon_
             (when (syntax->datum #'index-start-expanded)
               (throw-if-not-int-index #'getter.index-start #'index-start-expanded))
             (when (syntax->datum #'index-end-expanded)
               (throw-if-not-int-index #'getter.index-end #'index-end-expanded))))

         #| (displayln (format "~a func-inline-list:" name_)) |#
         #| (pretty-print (read-state func-inline-list)) |#

         ;; FIXME use func-inline-list in generation
         ;; NOTE: code generation part
         (with-syntax (;; NOTE: list of assignations of inlined functions called in the right-hand side
                       (func-inline-list-stx (datum->syntax stx `(#%datum ,(read-state func-inline-list))))
                       (name-symb (var-symbol 
                                   (symbol->string (syntax->datum #'name))
                                   src-pos)))
           #| (displayln "func-inline-list-stx") |#
           #| (pretty-print #'func-inline-list-stx) |#
           (let ((assign-label
                  (cond
                    ((and (equal? func-name name_) (equal? func-basic-type 'real))
                     'func-assign)
                    (else 'assign))))
             (cond
               ((syntax->datum #'getter.index)
                (begin
                  (when (equal? maybe-array-range '())
                    (raise-syntax-error #f (format "Variable ~a is not an array." name_) #'name))
                  (with-syntax* ((index-expanded
                                  (local-expand-opaque #'getter.index))
                                 (expanded-range 
                                   (cadr 
                                     (syntax->datum 
                                       (local-expand
                                         (datum->syntax #'getter.index maybe-array-range) 'expression (list))))))
                    ;; NOTE: (cadr (cadr ... because expanded looks like that: '(#%app '4)
                    ;; NOTE: can't do check in transmormer phase because loop variable is not resolved by that time.
                    ;;       Need to insert check-bounds code in generated code
                    (quasisyntax/loc
                        stx
                      (if (fx>= #,#'index-expanded #,#'expanded-range)
                          (raise-syntax-error
                           #f
                           (format "index ~a out of range: expect [0, ~a)" 
                                   #,#'index-expanded 
                                   #,#'expanded-range) 
                           #'#,stx)
                          (list 
                            #,#'func-inline-list-stx
                            (list '#,assign-label (list 'array-ref name-symb #,#'index-expanded) value-exp)))))))
               (slice-colon_
                (let ((slice-idx-symb 'slice_idx))
                  (when (equal? maybe-array-range '())
                    (raise-syntax-error #f (format "Variable ~a is not an array." name_) #'name))
                  (with-syntax* ((slice-idx-synt (datum->syntax stx slice-idx-symb))
                                 (expanded-range 
                                   (cadr 
                                     (syntax->datum 
                                       (local-expand (datum->syntax #'stx maybe-array-range) 'expression (list)))))
                                 (index-start-expanded_ (if (atom-number #'index-start-expanded)
                                                            (atom-number #'index-start-expanded)
                                                            #'index-start-expanded))
                                 (index-end-expanded_ (if (atom-number #'index-end-expanded)
                                                          (atom-number #'index-end-expanded)
                                                          #'index-end-expanded))
                                 (lvalue-slice-range #'(fx- index-end-expanded_
                                                            index-start-expanded_))
                                 (rvalue-slice-range value-range))
                    
                    (quasisyntax/loc
                        stx
                      (begin
                        (when (fx< #,#'index-start-expanded 0) 
                          (raise-syntax-error 
                            #f 
                            (format "index must be nonnegative, given: ~v" 
                                    #,#'index-start-expanded) 
                            #'#,stx))
                        (when #,#'rvalue-slice-range
                          (unless (fx= #,#'lvalue-slice-range #,#'rvalue-slice-range)
                            (raise-syntax-error 
                              #f 
                              (format "can not cast right-hand side range ~v to the left-hand side range ~v" 
                                      #,#'rvalue-slice-range 
                                      #,#'lvalue-slice-range) 
                              #'#,stx)))
                        (when (fx> #,#'index-end-expanded #,#'expanded-range)
                          (raise-syntax-error
                           #f
                           (format "slice end index ~a is out of range: expect [0, ~a)" 
                                   #,#'index-end-expanded 
                                   #,#'expanded-range) 
                           #'#,stx))                          
                        (list 
                          #,#'func-inline-list-stx
                          (for/list ((slice-idx-synt (in-range 0 #,#'lvalue-slice-range)))
                            (list '#,assign-label 
                                  (list 
                                    'array-ref 
                                    name-symb 
                                    (fx+ #,#'index-start-expanded slice-idx-synt)) 
                                  value-exp))))))))
               (else
                (quasisyntax/loc stx
                  (list 
                    #,#'func-inline-list-stx
                    (list '#,assign-label (list 'array-ref name-symb 0) value-exp))))))))))))


(define-syntax (decl-block stx)
  (syntax-parse stx
    [(_ body ...)
     (syntax/loc stx
                 (list body ...))]))

(define-syntax (var-decl stx)
  (syntax-parse stx
    (({~literal var-decl};; if expr is array get-value, then emit declaration and assigantion syntax objects  
                ((~literal type) ((~literal array-type) basic-type "[" num "]")) name:id (~seq "=" value:expr)) 
     (datum->syntax stx `(decl-block 
                           (var-decl (type (array-type ,#'basic-type "[" ,#'num "]")) ,#'name)
                           (assignation ,#'name "[" ":" "]" "=" ,#'value))))

    (({~literal var-decl}
      type:expr name:id (~optional (~seq "=" value:expr)
                                   #:defaults ((value #'notset))))
     (let* ((name_ (syntax->datum #'name))
            (type-datum (syntax->datum #'type))
            (type (to-landau-type stx (parse-type type-datum)))
            (name-vs (var-symbol (symbol->string name_) (syntax-position #'name))))
       (check-duplicate-variable-name name_ #'name)
       (add-variable! (syntax-parameter-value #'current-variables) #'name type)
       #| (add-real-var! #'name type stx (syntax-parameter-value #'real-vars-table)) |#
       
       (if (equal? (syntax->datum #'value) 'notset)
           (quasisyntax/loc stx 
                            (list 'var-decl #,name-vs '#,(to-landau-type stx type)))
           (quasisyntax/loc stx
             (list
               (list 'var-decl #,name-vs '#,(to-landau-type stx type))
               #,(datum->syntax
                   stx `(assignation ,#'name "=" ,#'value)) )))))))
                         

;; NOTE: Copied from semantics. Maybe need to reduce code duplication
(define-for-syntax (check-duplicate-constant name stx)
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "duplicate constant declaration" stx)))

(define-for-syntax (check-duplicate-parameter name stx)
  (when (hash-has-key? parameters name)
    ; (print parameters)
    (raise-syntax-error #f "duplicate parameter declaration" stx)))

;; TODO: change int to int-index and prohibit int arrays
(define-syntax (constant stx)
  (syntax-parse stx
    (({~literal constant} "const" type name:id "=" (~or value:expr (~seq "{" ({~literal parlist} arr-item*:par-spec ...) "}")))
     (let* ((type-parsed (parse-type (syntax->datum #'type)))
            (base-type (car type-parsed))
            (type (if (and (not (is-array-type? type-parsed))
                           (equal? base-type 'int))
                    (list 'int-index (cadr type-parsed))
                    type-parsed)))
          ;  (displayln #'(list par ...))
           
           (check-duplicate-constant (syntax->datum #'name) stx)
           (cond
             ((is-array-type? type-parsed)
              (begin
               (when (or (equal? base-type 'int)
                         (equal? base-type 'int-index))
                 (raise-syntax-error #f "integer constant arrays are not implemented" stx))
               (let*
                ((expanded-list
                  (for/list ((arr-item (in-list (syntax-e #'(arr-item*.par-value ...)))))
                            (begin
                             (local-expand-opaque arr-item)))
                  ))
                (hash-set! constants (syntax->datum #'name)
                           (constant 'constant-name-placeholder type 'constant-array-placeholder))
                (syntax/loc stx '()))))
             (else
              (let* ((expanded (local-expand-opaque #'value))
                     (value (atom-number expanded))
                     (value-type (syntax-property expanded 'landau-type))
                     (value
                      (cond
                        ((and (equal? base-type 'int-index) (equal? value-type 'real))
                         (raise-syntax-error #f "real value is not allowed for integer constant" stx))
                        ((and (equal? base-type 'real) (equal? value-type 'int-index))
                         (exact->inexact value))
                        (else value))))
                    (hash-set! constants (syntax->datum #'name)
                               (constant value type #f))
                    (syntax/loc stx '())))
             )))))
    

(define-syntax (parameter stx)
  (syntax-parse stx
    (({~literal parameter} "parameter" "[" size "]" name:id)
     
     (check-duplicate-parameter (syntax->datum #'name) stx)
     (let* ((expanded-size (local-expand-opaque #'size)))
              
       (hash-set! parameters (syntax->datum #'name)
                  (list expanded-size))
       (syntax/loc stx '())))))

(define-syntax (print stx)
  (syntax-parse stx
    (({~literal print} "print" str expr)
     (begin
       ;; FIXME: Check if slice
       (syntax/loc stx '())))
    (({~literal print} "print" expr)
     (begin
       ;; FIXME: Check if slice
       (syntax/loc stx '())))))

(define-syntax (for-expr stx)
  (syntax-parse stx
[ (_ "for" ({~literal loop-var-decl} name:id "=" "[" start:expr ":" stop:expr "]") pat:expr-body-cls)
      (let ((name_ (syntax->datum #'name))
            (type (list 'int-index '())))
        (check-duplicate-variable-name name_ #'name)
        (let* ((new-vars-layer (new-variables-nesting
                                (syntax-parameter-value #'current-variables)))
               (sym (add-variable! new-vars-layer #'name type))
              
               (err-msg "Only integer expressions of constant and loop variables are allowed for ranges"))
          (with-syntax ([symm (datum->syntax stx sym)]
                        [start-val (local-expand-opaque #'start)]
                        [stop-val (local-expand-opaque #'stop)])
            (throw-if-not-type 'int-index #'start-val err-msg)
            (throw-if-not-type 'int-index #'stop-val err-msg)
            (quasisyntax/loc stx
              (begin
                ;(writeln (format "stop-val: ~v" (syntax->datum #'stop-val)))
                (syntax-parameterize
                    [(current-variables
                      '#,new-vars-layer)]
                  (for/list ([symm (in-range start-val stop-val)])
                    pat.body)))))))]))
    

(define-syntax (bool-expr stx)
  (syntax-parse stx
    [(_ b-expr1 "or" b-expr2) #'(or b-expr1 b-expr2)]
    [(_ b-term) #'b-term]))
    

(define-syntax (bool-term stx)
  (syntax-parse stx
    [(_ b-term1 "and" b-term2) #'(and b-term1 b-term2)]
    [(_ b-factor) #'b-factor]))

(define-syntax (bool-factor stx)
  (syntax-parse stx
    [(_ "(" b-expr ")") #'b-expr]
    [(_ "(" n1 ({~literal comp-op} op) n2 ")")
     (with-syntax* (
                    (val1 (local-expand-opaque #'n1))
                    (val2 (local-expand-opaque #'n2)))
       (throw-if-not-type 'int-index #'val1 "Only 'int allowed inside the `if` condition.")
       (throw-if-not-type 'int-index #'val2 "Only 'int allowed inside the `if` condition.")
       (quasisyntax/loc stx #,(match (syntax->datum #'op)
                                ["==" (syntax/loc stx (equal? n1 n2))]
                                ["!=" (syntax/loc stx (not(equal? n1 n2)))]
                                [">" (syntax/loc stx (> n1 n2))]
                                [">=" (syntax/loc stx (>= n1 n2))]
                                ["<=" (syntax/loc stx (<= n1 n2))]
                                ["<" (syntax/loc stx (< n1 n2))])))]
                    
    [(_ "not" b-factor) #'(not b-factor)]))
    

; TODO: Do I need to check if it is not using outer variables?
(define-syntax (if-expr stx)
  (syntax-parse stx
    [({~literal if-expr} "if" b-expr "{" body-true ... "}" "else" "{" body-false ... "}")
     (syntax/loc stx
       (match b-expr
         [#f (list body-false ...)]
         [#t (list body-true ...)]
         [else (raise-syntax-error #f "not a boolean in if-expr")]))]
         
    [({~literal if-expr} "if" b-expr "{" body ... "}")
     (syntax/loc stx
       (match b-expr
         [#f '()]
         [#t (list body ...)]
         [else (raise-syntax-error #f "not a boolean in if-expr")]))]))
         
;; FIXME: check variable shadowing other functions names
(define-for-syntax (check-duplicate-variable-name name stx-name)
  #| (displayln "FIXME: commented argument shadowing check") |#
  (when (search-variable name (syntax-parameter-value #'current-variables))
    (raise-syntax-error #f "duplicate variable declaration" stx-name))
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "variable name shadows constant" stx-name))
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "variable name shadows parameter" stx-name))
  #| (when (hash-has-key? (syntax-parameter-value #'current-arguments) name) |#
  #|   (raise-syntax-error #f "variable name shadows argument" stx-name)) |#
  (when (equal? (syntax-parameter-value #'function-name) name)
    (raise-syntax-error #f "variable name shadows function" stx-name)))
  

(define-for-syntax (check-duplicate-argument-name args funcname name stx-name)
  (when (equal? funcname name)
    (raise-syntax-error #f "argument name shadows function" stx-name))
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "variable name shadows parameter" stx-name))
  (when (hash-has-key? args name)
    (raise-syntax-error #f "duplicate argument" stx-name))
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "argument name shadows constant" stx-name)))


(define-for-syntax (search-function name)
  (let*
      ((name_ (syntax->datum name))
       (err-msg "Derivative can be written only to the real function return variable")
       (type_
        (cond
          ((hash-has-key? constants name_)
           (raise-syntax-error
            #f
            err-msg #'name))
          (else
           (let ((var (search-variable
                       name_ (syntax-parameter-value #'current-variables))))
             (cond
               (var
                (raise-syntax-error
                 #f
                 err-msg #'name))
               (else
                (let ((arg (search-argument
                            name_ (syntax-parameter-value #'current-arguments))))
                  (cond
                    (arg
                     (raise-syntax-error
                      #f
                      err-msg #'name))
                    ((equal? (syntax-parameter-value #'function-name) (syntax->datum name))
                     (syntax-parameter-value #'function-return-type))
                    (else
                     (raise-syntax-error #f "name not found" name)))))))))))
    type_))
  


(define-for-syntax (search-arg-or-parameter name)
  (let*
      ((name_ (syntax->datum name))
       (type_
        (cond
          ((hash-has-key? constants name_)
           (raise-syntax-error
            #f
            "constant can not be inside the derivative anotation or application. Expect a real argument." name))
          (else
           (let ((var (search-variable
                       name_ (syntax-parameter-value #'current-variables))))
             (cond
               (var
                (raise-syntax-error
                 #f
                 "Variable can not be inside the derivative anotation or application. Expect a real argument." name))
               (else
                (let ((arg (search-argument
                            name_ (syntax-parameter-value #'current-arguments))))
                  (cond
                    (arg
                     (argument-type arg))
                    ((equal? (syntax-parameter-value #'function-name) (syntax->datum name))
                     (syntax-parameter-value #'function-return-type))
                    ((hash-has-key? parameters name_) (list 'real (hash-ref parameters name_)))
                    (else
                     (raise-syntax-error #f "name not found" name)))))))))))
    type_))
  

(define-for-syntax (search-var name)
  (let*
      ((name_ (syntax->datum name))
       (type_
        (cond
          ((hash-has-key? constants name_)
           (raise-syntax-error
            #f
            "constant can not be inside the derivative annotation or discard. Expect a real argument or a variable."
            name))
          (else
           (let ((var (search-variable
                       name_ (syntax-parameter-value #'current-variables))))
             (cond
               (var
                (variable-type var))
               (else
                (let ((arg (search-argument
                            name_ (syntax-parameter-value #'current-arguments))))
                  (cond
                    (arg
                     (argument-type arg))
                    ((equal? (syntax-parameter-value #'function-name) (syntax->datum name))
                     (syntax-parameter-value #'function-return-type))
                    (else
                     (raise-syntax-error #f "name not found" name)))))))))))
      type_))
 

(define-for-syntax (check-if-real name type err-msg)
  (unless (equal? type 'real)
      (raise-syntax-error #f err-msg name)))


(define-for-syntax (check-real-func-existence name (get-range #f))
  (define err-msg "Derivative can be written only to the real function return variable")
  (define type (any->expanded-type (search-function name) name))
  (define type-base (get-type-base type))
  (define type-range (get-type-range type))
  (check-if-real name type-base err-msg)
  (if get-range
    (if (equal? type-range '()) 
      #f 
      type-range)
    type))


;; TODO same function as check-real-func-existence except err-msg. Refactor.
(define-for-syntax (check-real-arg-or-parameter-existence name [get-range #f])
 (define err-msg "only real function arguments are allowed be inside the derivative anotation or application")
 (define type-like-anything (search-arg-or-parameter name))
 (define type (any->expanded-type type-like-anything name))
 (define type-base (get-type-base type))
 (define type-range (get-type-range type))
 (check-if-real name type-base err-msg)
 (if get-range
   (if (equal? type-range '()) 
     #f 
     type-range)
   type))


;; TODO same function as check-real-func-existence except err-msg. Refactor.
(define-for-syntax (check-real-var-existence name [get-range #f])
 (define err-msg "only real arguments or variables are allowed be inside the derivative discard")
 (define type (any->expanded-type (search-var name) name))
 (define type-base (get-type-base type))
 (define type-range (get-type-range type))
 (check-if-real name type-base err-msg)
 (if get-range
   (if (equal? type-range '()) 
     #f 
     type-range)
   type))


(define-syntax (get-derivative stx)
  (syntax-parse stx
    ((_ get-value-1 "'" get-value-2)
     (syntax/loc stx (list 'get-derivative get-value-1 get-value-1)))))
;; NOTE:
;; non-array dx arguments are casted to single-cell array

