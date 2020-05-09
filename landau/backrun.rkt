#lang racket
(require (for-syntax racket/base
                     syntax/parse
                     racket/syntax
                     "environment.rkt"
                     racket/flonum
                     racket/function
                     racket/fixnum
                     racket/match
                     racket/set
                     racket/list
                     racket/contract)
         (only-in "process-actions-list.rkt" process)
         "environment.rkt"
         "type-utils.rkt"
         racket/contract/region
         racket/stxparam
         racket/flonum
         racket/fixnum
         "common-for-syntax.rkt")

(define-syntax-parameter current-variables #f)
(define-syntax-parameter function-name #f)
(define-syntax-parameter function-return-value #f)
(define-syntax-parameter function-return-type  #f)
(define-syntax-parameter current-arguments #f)
(define-syntax-parameter dx-names-set #f)
(define-syntax-parameter real-vars-table #f)

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

(define-for-syntax (is-int-index stx)
  (is-type stx 'int-index))

(define-for-syntax (is-func? stx)
  (match (syntax->datum stx)
    ((list 'func _ ...) #t)
    (_ #f)))


(define-syntax (program stx)
  (syntax-parse stx
    [(_ body ...)
     (let* ((program-terms (syntax-e #'(body ...)))
            (top-level-decl (filter (compose not is-func?) program-terms))
            (funcs (filter is-func? program-terms))
            (fake-src-pos 0)
            (tld (if (empty? top-level-decl)
                   #'(void)
                   (local-expand (datum->syntax stx (cons 'begin top-level-decl)) 'expression '()))))
       ;; NOTE: Evaluate top level constants and parameters declarations
       (eval tld (current-namespace))
       (datum->syntax
        stx
        (cons 'list 
              (list funcs-info-GLOBAL 
                    (cons 'list 
                          (for/list ((fnc (in-list funcs)))
                            (let ((fnc-name (symbol->string
                                             (syntax->datum (caddr (syntax-e fnc)))))
                                  (dx-names-set-val (make-hash))
                                  (real-vars-table-val (make-hash)))
                              #`(syntax-parameterize ((dx-names-set '#,dx-names-set-val)
                                                      (real-vars-table '#,real-vars-table-val))
                                  (cons (var-symbol #,fnc-name #,fake-src-pos) (call-with-values
                                                    (thunk 
                                                     (process 
                                                      ;; FIXME: Dont know why, but it works
                                                      (list (list) #,fnc) 
                                                      #,dx-names-set-val 
                                                      #,real-vars-table-val))
                                                    list)))))))))
       )]))
           
(define (throw-if-bundle-mistype bundl-1 bundl-2)
  (unless (equal? (car bundl-1) (car bundl-2))
    (error "bug: bundles have different type")))
  

(define-for-syntax (expand-range rng-stx rng)
  (if rng
      (begin
        (cadr (syntax->datum
               (local-expand (datum->syntax #'rng-stx rng) 'expression '()))))
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
     ({~literal argument} type:type-spec name:id))
     
    (pattern
     ({~literal other-argument} "," type:type-spec name:id))))
     

(define-syntax (func stx)
  (syntax-parse stx
    (({~literal func} type:expr name:id "(" ({~literal arglist} arg*:arg-spec ...) ")"
                      "{" body "}")
     (let* ((args (make-hash))
            (func-return-type (parse-type (syntax->datum #'type)))
            (func-name-str (symbol->string (syntax->datum #'name)))
            (func-return-value (gensym func-name-str))
            (fake-src-pos 0)
            (func-return-range (atom-number (type-range (expand-type (datum->syntax #'type func-return-type)))))
            (func-args
             (for/list ((argname (in-list (syntax-e #'(arg*.name ...))))
                        (argtype (in-list (syntax->datum #'(arg*.type ...)))))
                       (check-duplicate-argument-name args (syntax->datum #'name)
                                                      (syntax->datum argname) argname)
                       (datum->syntax argname (add-argument! args (syntax->datum argname)
                                                             (parse-type argtype)))
                       (add-real-var! (datum->syntax #f (syntax->datum argname)) 
                                      (parse-type argtype) 
                                      stx 
                                      (syntax-parameter-value #'real-vars-table))

                       (list (symbol->string (syntax->datum argname)) (parse-type-to-syntax argtype))))
            (func-args-expanded
             (for/list ((arg (in-list func-args)))
                       (let* ((arg-name (car arg))
                              (arg-type (syntax-e (cadr arg)))
                              (arg-base-type (syntax->datum (cadr arg-type)))
                              (arg-range-expanded
                               (atom-number
                                (local-expand
                                 (datum->syntax stx (syntax->datum (caddr arg-type))) 'expression '()))))
                             (cons arg-base-type arg-range-expanded)))))

       (when (hash-has-key? BUILT-IN-FUNCTIONS func-name-str)
         (raise-syntax-error #f "function name shadows built-in function" #'name))

       (hash-set!
        funcs-info-GLOBAL
        (var-symbol func-name-str fake-src-pos)
        (func-info func-return-value (length func-args) func-args-expanded (car func-return-type) func-return-range))
       (unless (equal? (car func-return-type) 'real)
         (raise-syntax-error #f "function can return only reals" #'type))
       (add-real-var! #'name func-return-type stx (syntax-parameter-value #'real-vars-table))
       
       (quasisyntax/loc stx
         (syntax-parameterize
             ((current-variables
               (new-variables-nesting
                (syntax-parameter-value #'current-variables)))
              (current-arguments '#,args)
              (function-name (syntax->datum #'name))
              (function-return-value '#,func-return-value)
              (function-return-type '#,func-return-type))
           body))))))

(begin-for-syntax
 (define/contract 
  (search-name-backrun stx name)
  (-> (syntax/c any/c) (syntax/c symbol?)
      ;; NOTE: real constants are not propagated in backrun
      (values (syntax/c (or/c symbol? integer? false/c)) type/c boolean? integer?))
  (let
   ((fake-src-pos 0)
    (is-const #t)
    (is-not-const #f)
    (name_ (syntax->datum name))
    (func-name (syntax-parameter-value #'function-name)))
     (let-values
      (((value type const? src-pos)
        (cond
          ((hash-has-key? constants name_)
           (let ((c (hash-ref constants name_)))
             (values
              (datum->syntax stx
                             (constant-value c))
              (constant-type c)
              is-const
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
                         (variable-type var) #f (variable-src-pos var)))
                ((equal? name_ func-name)
                 (let ((func-return-value (syntax-parameter-value #'function-return-value)))
                   (values
                    (datum->syntax stx func-return-value)
                    (syntax-parameter-value #'function-return-type)
                    is-not-const
                    fake-src-pos)))
                (else
                 (let ((arg (search-argument
                             name_ (syntax-parameter-value #'current-arguments))))
                   (cond
                     (arg
                      (values
                       (datum->syntax stx (argument-symbol arg))
                       (argument-type arg)
                       is-not-const
                       fake-src-pos))
                     (else
                      (cond
                        ((hash-has-key? parameters name_) (raise-syntax-error #f "parameters can not be used in expression" name))
                        (else (raise-syntax-error #f "name not found" name))))))))))))))
      ; (println (format "~a -> ~a" name value))
      (values value type const? src-pos))

)))

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
         (((getter) (make-getter-info (if (attribute index) #'index #'#f)
                                      (if (attribute slice-colon) #'slice-colon #'#f)))
          ((name_) (syntax->datum #'name))
          ((name-str) (symbol->string name_))
          ((index_) (if (attribute index) (syntax->datum #'index) #f))
          ((slice-colon_) (attribute slice-colon))
          ;; NOTE: stx-pos is provided for variables only
          ((value full-type const? src-pos) (search-name-backrun stx #'name))
          ((maybe-array-range) (cadr full-type))
          ((base-type) (car full-type))
          
          ((array-range) (atom-number #`#,(expand-range 
                                                 (if (attribute index) #'index #'#f)
                                                 (if (equal? maybe-array-range '()) #f maybe-array-range))))
          ((index-expanded) (if (attribute index) 
                              (local-expand 
                                       #'index
                                       'expression '())
                              #f)))
      ;  (displayln (format "~a ~a" #'name (syntax-position #'name)))
       (check-proper-getter 
        (if (attribute index) #'index #'#f) (attribute slice-colon) array-range index-expanded name-str stx)
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
           result))))))

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
     (let* ((func-getter (make-getter-info (if (attribute func-ret-idx) #'func-ret-idx #'#f)
                                           (if (attribute func-slice-colon) #'func-slice-colon #'#f)))
            (dx-getter (make-getter-info (if (attribute dx-idx) #'dx-idx #'#f)
                                     (if (attribute dx-slice-colon) #'dx-slice-colon #'#f)))
            (name-1-dat (symbol->string (syntax->datum #'func-name)))
            (name-2-dat (symbol->string (syntax->datum #'dx-name)))
            (arr-1-range (check-real-func-existence #'func-name 'get-range))
            (arr-2-range (check-real-arg-or-parameter-existence #'dx-name 'get-range))
            (func-ret-idx-expanded (if (attribute func-ret-idx) (local-expand #'func-ret-idx 'expression '()) #f))
            (dx-idx-expanded (if (attribute dx-idx) (local-expand #'dx-idx 'expression '()) #f))
            (expanded-range (atom-number #`#,(expand-range (if (attribute func-ret-idx) #'func-ret-idx #'#f) (if (equal? arr-1-range '()) #f arr-1-range))))
            (expanded-range-dx (atom-number #`#,(expand-range (if (attribute dx-idx) #'dx-idx #'#f) (if (equal? arr-2-range '()) #f arr-2-range))))
            (value-type (syntax-property (local-expand #'value 'expression '()) 'landau-type))
            (func-slice-colon_ (attribute func-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon))
            (throw-slice-cast-error (thunk
                                     (raise-syntax-error
                                      #f
                                       (format "the right-hand side is slice, but the left-hand one is not") stx))))
       ;  (println "der-apply")
       (check-proper-getter (if (attribute func-ret-idx) #'func-ret-idx #'#f) (attribute func-slice-colon) expanded-range func-ret-idx-expanded name-1-dat stx)
       (check-proper-getter (if (attribute dx-idx) #'dx-idx #'#f) (attribute dx-slice-colon) expanded-range-dx dx-idx-expanded name-2-dat stx)
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
                      (quasisyntax/loc
                          stx
                        (begin
                          #,#'func-index-range-checks
                          (unless (fx= rvalue-outer-prod-range 1)
                            (unless (fx= rvalue-outer-prod-range func-slice-range)
                             (raise-syntax-error #f 
                                                (format "can not cast right-hand side range ~v to the left-hand side range ~v" rvalue-outer-prod-range func-slice-range) #'#,stx)))
                          (list 'nested (for/list ((#,#'slice-idx-synt (in-range 0 #,#'func-slice-range)))
                                          (list 'der-apply
                                                (car #,#'value)
                                                (list 'array-ref #,#'dx-name-str 0))))))))

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
                                (raise-syntax-error #f 
                                                (format "can not cast right-hand side range ~v to the left-hand side range ~v" value-slice-range func-slice-range) #'#,stx)))
                          (list 'nested (for/list ((#,#'slice-idx-synt (in-range 0 #,#'func-slice-range)))
                                          (list 'der-apply
                                                (car #,#'value)
                                                (list 'array-ref #,#'dx-name-str #,#'dx-idx-expanded))))))))

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
                            (raise-syntax-error #f 
                                                (format "can not cast right-hand side range ~v to the left-hand side range ~v" rvalue-outer-prod-range func-slice-range) #'#,stx))
                          (list 'nested (for/list ((#,#'slice-idx-synt (in-range 0 #,#'value-slice-range)))
                                          (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                                          (list 'der-apply
                                                                (car #,#'value)
                                                                (list 'array-ref #,#'dx-name-str #,#'dx-full-idx)))))))))))
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
            (df-idx-expanded (if (attribute df-idx) (local-expand #'df-idx 'expression '()) #f))
            (dx-idx-expanded (if (attribute dx-idx) (local-expand #'dx-idx 'expression '()) #f))
            (expanded-range-df (atom-number #`#,(expand-range (if (attribute df-idx) #'df-idx #'#f) (if (equal? df-range '()) #f df-range))))
            (expanded-range-dx (atom-number #`#,(expand-range (if (attribute dx-idx) #'dx-idx #'#f) (if (equal? dx-range '()) #f dx-range))))
            (df-slice-colon_ (attribute df-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon)))

       (check-proper-getter (if (attribute df-idx) #'df-idx #'#f) (attribute df-slice-colon) expanded-range-df df-idx-expanded name-1-dat stx)
       (check-proper-getter (if (attribute dx-idx) #'dx-idx #'#f) (attribute dx-slice-colon) expanded-range-dx dx-idx-expanded name-2-dat stx)
       
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
                     (format "Can not override annotation: ~a[~a] ' ~a[~a] == 1.0" df-name-symb #,df-idx-expanded dx-name-symb #,dx-idx-expanded) #'#,stx))
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
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                       (list 'discard
                                             (list 'array-ref #,#'df-name-symb #,#'df-idx-expanded)
                                             (list 'array-ref #,#'dx-name-symb #,#'full-idx))))))))
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
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'discard
                                             (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                             (list 'array-ref #,#'dx-name-symb #,#'dx-idx-expanded))))))))
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
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                                       (list 'discard
                                                             (list 'array-ref #,#'df-name-symb #,#'df-full-idx)
                                                             (list 'array-ref #,#'dx-name-symb #,#'dx-full-idx))))))))))
                (df-slice-colon_
                 (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   ;(displayln (format "full-idx ~a\nlvalue-outer-prod-range ~a\nr-value-exp ~a\nr-value-type ~a" #'full-idx #'lvalue-outer-prod-range #'r-value-exp #'r-value-type))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'discard
                                             (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                             (list 'array-ref #,#'dx-name-symb 0))))))))
                (dx-slice-colon_
                 (with-syntax*
                     ((dx-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'dx-slice-range)
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                       (list 'discard
                                             (list 'array-ref #,#'df-name-symb 0)
                                             (list 'array-ref #,#'dx-name-symb #,#'full-idx))))))))
                (else 
                 (begin
                   (throw-if-r-value-is-slice stx #'der-value)
                   (quasisyntax/loc stx
                     (list 'discard
                           (list 'array-ref #,#'df-name-symb 0)
                           (list 'array-ref #,#'dx-name-symb 0)))))))))))]))

(define-for-syntax (get-expanded-slice-index-start index-start)
  (let* ((index-start-expanded (if (syntax->datum index-start)
                                   (local-expand index-start 'expression '())
                                   (is-type_ 'int-index #'0)))
         (index-start-expanded_ (if (atom-number index-start-expanded)
                                    (atom-number index-start-expanded)
                                    index-start-expanded)))
    (values index-start-expanded index-start-expanded_)))

(define-for-syntax (get-expanded-slice-index-end stx index-end array-range-stx)
  (let* ((index-end-expanded (if (syntax->datum index-end)
                                 (begin
                                   ; (println index-end)
                                   (local-expand index-end 'expression '()))
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
    (with-syntax*
        ((slice-range (quasisyntax/loc stx 
                        (fx- #,dx-index-end-number
                             #,dx-index-start-number))))
      (throw-if-not-type 'int-index dx-index-start-stx "Indexes expected to be expressions with 'int constants and/or loop variables")
      (throw-if-not-type 'int-index dx-index-end-stx "Indexes expected to be expressions with 'int constants and/or loop variables")
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
      ((r-value-exp (local-expand r-value-stx 'expression '()))
       (r-value-type (syntax-property #'r-value-exp 'landau-type)))
    (when (is-slice? (syntax->datum #'r-value-type))
      (raise-syntax-error #f
                          (format "the right-hand side is slice, but the left-hand one is not") stx))))

(define-for-syntax (range-casting-check stx r-value-type lvalue-outer-prod-range)
  (when (is-slice? (syntax->datum r-value-type))
    (with-syntax ((rvalue-slice-range (get-slice-range (syntax-e r-value-type))))
      #`(unless (fx= #,lvalue-outer-prod-range #,#'rvalue-slice-range)
          (raise-syntax-error #f 
                              (format "can not cast right-hand side range ~v to the left-hand side range ~v" #,#'rvalue-slice-range #,lvalue-outer-prod-range) #'#,stx)))))

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
     (let* ((df-getter (make-getter-info (if (attribute df-idx) #'df-idx #'#f)
                                     (if (attribute df-slice-colon) #'df-slice-colon #'#f)))
            (dx-getter (make-getter-info (if (attribute dx-idx) #'dx-idx #'#f)
                                     (if (attribute dx-slice-colon) #'dx-slice-colon #'#f)))
            (name-1-dat (symbol->string (syntax->datum #'df-name)))
            (name-2-dat (symbol->string (syntax->datum #'dx-name)))
            (df-range (check-real-var-existence #'df-name 'get-range))
            (dx-range (check-real-arg-or-parameter-existence #'dx-name 'get-range))
            (df-idx-expanded (if (attribute df-idx) (local-expand #'df-idx 'expression '()) #f))
            (dx-idx-expanded (if (attribute dx-idx) (local-expand #'dx-idx 'expression '()) #f))
            (expanded-range-df (atom-number #`#,(expand-range (if (attribute df-idx) #'df-idx #'#f) (if (equal? df-range '()) #f df-range))))
            (expanded-range-dx (atom-number #`#,(expand-range (if (attribute dx-idx) #'dx-idx #'#f) (if (equal? dx-range '()) #f dx-range))))
            (df-slice-colon_ (attribute df-slice-colon))
            (dx-slice-colon_ (attribute dx-slice-colon)))

       (check-proper-getter (if (attribute df-idx) #'df-idx #'#f) (attribute df-slice-colon) expanded-range-df df-idx-expanded name-1-dat stx)
       (check-proper-getter (if (attribute dx-idx) #'dx-idx #'#f) (attribute dx-slice-colon) expanded-range-dx dx-idx-expanded name-2-dat stx)
       ;  (println "der-annot")
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
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (let ((#,#'slice-idx-synt 0))
                         der-value
                         (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                       (let ((#,#'slice-idx-synt #,#'dx-slice-idx-synt))
                                          (list 'der-annot
                                             (list 'array-ref #,#'df-name-symb 0)
                                             (list 'array-ref #,#'dx-name-symb #,#'full-idx))))))))))

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
                   (r-value-exp (local-expand #'der-value 'expression '()))
                   (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                (quasisyntax/loc
                    stx
                  (begin
                    #,#'dx-index-range-checks
                    #,(range-check-stx stx df-idx-expanded expanded-range-df)
                    #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                    (let ((#,#'slice-idx-synt 0))
                         der-value
                         (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                    (list 'der-annot
                                          (list 'array-ref #,#'df-name-symb #,#'df-idx-expanded-stx)
                                          (list 'array-ref #,#'dx-name-symb #,#'full-idx)))))))))

             ((list 'slice 'var "<- der-value")
              (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       (let ((#,#'slice-idx-synt 0))
                            der-value
                            (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'der-annot
                                             (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                             (list 'array-ref #,#'dx-name-symb 0)))))))))

             ((list 'slice 'cell "<- der-value")
              (with-syntax*
                     ((dx-idx-expanded-stx dx-idx-expanded) 
                      (df-slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (lvalue-outer-prod-range #'df-slice-range)
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,(range-check-stx stx dx-idx-expanded expanded-range-dx)
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       ;; assign placeholder value to slice variable which result is ignored. 
                       (let ((#,#'slice-idx-synt 0))
                         der-value
                         (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'der-annot
                                             (list 'array-ref #,#'df-name-symb #,#'full-idx)
                                             (list 'array-ref #,#'dx-name-symb #,#'dx-idx-expanded-stx)))))))))

             ((list 'slice 'slice "<- der-value")
              (with-syntax*
                     ((df-slice-idx-synt (datum->syntax stx df-slice-idx-name-GLOBAL))
                      (dx-slice-idx-synt (datum->syntax stx dx-slice-idx-name-GLOBAL))
                      (slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                      (df-full-idx #`(+ #,df-index-start-stx df-slice-idx-synt))
                      (dx-full-idx #`(+ #,dx-index-start-stx dx-slice-idx-synt))
                      (lvalue-outer-prod-range #'(fx* df-slice-range dx-slice-range))
                      (r-value-exp (local-expand #'der-value 'expression '()))
                      (r-value-type (syntax-property #'r-value-exp 'landau-type)))
                   (quasisyntax/loc
                       stx
                     (begin
                       #,#'df-index-range-checks
                       #,#'dx-index-range-checks
                       #,(range-casting-check stx #'r-value-type #'lvalue-outer-prod-range)
                       ;; assign placeholder value to slice variable which result is ignored. 
                       (let ((#,#'slice-idx-synt 0))
                              der-value ;; to trigger range checks
                              (list 'nested (for/list ((#,#'df-slice-idx-synt (in-range 0 #,#'df-slice-range)))
                                       (list 'nested (for/list ((#,#'dx-slice-idx-synt (in-range 0 #,#'dx-slice-range)))
                                                               (list 'der-annot
                                                                     (list 'array-ref #,#'df-name-symb #,#'df-full-idx)
                                                                     (list 'array-ref #,#'dx-name-symb #,#'dx-full-idx))))))))))))
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
       (list 'nested
             (list body ...)))]))

(define-syntax (func-body stx)
  (syntax-parse stx
    [(_ body ...)
     (syntax/loc stx
       (list 'nested
             (list body ...)))]))

(define-for-syntax (binary-op-cast op1 op2 stx)
  (let ((type1 (syntax-property op1 'landau-type))
        (type2 (syntax-property op2 'landau-type))
        (op1-atom (atom-number op1))
        (op2-atom (atom-number op2)))
    ; (println (format "~a ~a" type1 type2))
    (cond
      ((and (equal? type1 'int) (equal? type2 'int))
       (values op1 op2 'int))
      ((and (equal? type1 'real) (equal? type2 'real))
       (values op1 op2 'real))
      ((and (equal? type1 'int) (equal? type2 'real))
       (values (datum->syntax
                op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
               op2 'real))
      ((and (equal? type1 'real) (equal? type2 'int))
       (values op1 (datum->syntax
                    op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
               'real))
      ((and (equal? type1 'int-index) (equal? type2 'int-index))
       (values op1 op2 'int-index))
      ((and (equal? type1 'int) (equal? type2 'int-index))
       (values op1 op2 'int))
      ((and (equal? type1 'int-index) (equal? type2 'int))
       (values op1 op2 'int))
      ((and (equal? type1 'real) (equal? type2 'int-index))
       (values op1 (datum->syntax
                    op2 (if op2-atom (exact->inexact op2-atom) `(exact->inexact ,op2)) op2)
               'real))
      ((and (equal? type1 'int-index) (equal? type2 'real))
       (values (datum->syntax
                op1 (if op1-atom (exact->inexact op1-atom) `(exact->inexact ,op1)) op1)
               op2 'real))
      ((and (is-slice? type1) (is-slice? type2))
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
            (raise-syntax-error #f (format "cannot cast ~v and ~v" slice-type1 slice-type2) stx)))))
      ((and (is-slice? type1) (not (is-slice? type2)))
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
      ((and (is-slice? type2) (not (is-slice? type1)))
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
       (raise-syntax-error #f (format "cannot cast ~v and ~v. op1: ~v, op2: ~v" type1 type2 op1 op2) stx)))))


(define-syntax (expr stx)
  (syntax-parse stx
    ((_ expr1 op:plus-or-minus expr2)
     (let* ((expr1-expanded (local-expand #'expr1 'expression '()))
            (expr2-expanded (local-expand #'expr2 'expression '()))
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
     (with-syntax ((expanded (local-expand #'expr1 'expression '())))
       (syntax-track-origin (syntax/loc stx expanded) #'expr1 #'expr))]))

(begin-for-syntax
  (define/contract 
    (arity-check func-name n-expected n-provided stx)
    (-> string? integer? integer? (syntax/c any/c) void?)
    (unless (equal? n-expected n-provided)
      (raise-syntax-error
       #f
       (format "function ~a expects ~a parameters, but ~a provided" func-name n-expected n-provided)
       stx))))
    

(define-syntax (func-call stx)
  (syntax-parse stx
    ((_ function-name "(" ({~literal parlist} par*:par-spec ...) ")")
     (let* ((func-str (symbol->string (syntax->datum #'function-name)))
            (fake-src-pos 0)
            (func-name-vs (var-symbol func-str fake-src-pos))
            (func-pars-list (for/list ((par-value (in-list (syntax-e #'(par*.par-value ...)))))
                                      par-value))
            (par-list-len (length func-pars-list))
            (builtin-functions (hash-keys BUILT-IN-FUNCTIONS)))
       ;; FIXME: Check func name, type, casting conditions
       
       (if (member func-str builtin-functions)
           (let ((expected-arity (length (hash-ref BUILT-IN-FUNCTIONS func-str))))
             (if (equal? func-str "pow")
               (begin
                (arity-check func-str expected-arity par-list-len stx)
                (with-syntax* ((expr1 (car func-pars-list))
                               (expr2 (cadr func-pars-list))
                               (expanded1 (local-expand #'expr1 'expression '()))
                               (expanded2 (local-expand #'expr2 'expression '())))
                  (is-type_ 'real (syntax/loc stx (append expanded1 expanded2)))))
               (begin
                (arity-check func-str expected-arity par-list-len stx)
                (with-syntax* ((expr1 (car func-pars-list))
                               (expanded1 (local-expand #'expr1 'expression '())))
                  (is-type_ 'real (syntax/loc stx expanded1))))))
           ;; FIXME: perform checks
           (let*
               ((func-info (if (hash-has-key? funcs-info-GLOBAL func-name-vs)
                             (hash-ref funcs-info-GLOBAL func-name-vs)
                             (raise-syntax-error #f (format "function ~a is not defined." func-str) stx)))
                (type (if (func-info-output-range func-info)
                        (landau-type (func-info-output-base-type func-info)
                                     (func-info-output-range func-info))
                        (landau-type (func-info-output-base-type func-info)))))
             (with-syntax
                 ((slice-idx-synt (datum->syntax stx slice-idx-name-GLOBAL))
                  (pars-list-expanded (for/list ((par (in-list func-pars-list)))
                                                (local-expand par 'expression '()))))
               (is-type_ type
                         (quasisyntax/loc
                                    stx
                                    (append #,@#'pars-list-expanded))))))))))


(define-syntax (term stx)
  (syntax-parse stx
    ((_ expr1 op:mul-or-div expr2)
     (let* ((expr1-expanded (local-expand #'expr1 'expression '()))
            (expr2-expanded (local-expand #'expr2 'expression '()))
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
     (with-syntax ((expanded (local-expand #'term1 'expression '())))
       ;  (println #'expanded)
       ;  (println (format "term type: ~a" (syntax-property #'expanded 'landau-type)))
       ;  (println (format "term type emitted: ~a" (syntax-property (syntax-track-origin (syntax/loc stx expanded) #'term1 #'expr) 'landau-type)))
       ;  (println "")
       (syntax-track-origin (syntax/loc stx expanded) #'term1 #'expr))]))
    

(define-syntax (element stx)
  (syntax-parse stx
    [(_ number)
     (with-syntax* ((expanded (local-expand #'number 'expression '()))
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
    ; (let* ((primary-expanded (local-expand #'primary 'expression '()))
    ;        (factor-expanded (local-expand #'factor 'expression '()))
    ;        (primary-is-int-index (equal? (syntax-property expr1-expanded 'landau-type) 'int-index))
    ;        (factor-is-int-index (equal? (syntax-property expr2-expanded 'landau-type) 'int-index)))
    ;    (let*-values (((primary-casted factor-casted type) (binary-op-cast primary-expanded factor-expanded stx))
    ;                  ((primary-atom) (atom-number primary-casted))
    ;                  ((factor-atom) (atom-number factor-casted)))
    ;      (cond
    ;        ((equal? type 'real)
    ;         (is-real
    ;          (if (and primary-atom factor-atom)
    ;              (quasisyntax/loc stx #,(rlexpt primary-atom factor-atom))
    ;              (with-syntax* ((e1 (if (and (equal? #f primary-atom) (not primary-is-int-index))
    ;                                     primary-casted
    ;                                     #'(list)))
    ;                             (e2 (if (and (equal? #f factor-atom) (not factor-is-int-index))
    ;                                     factor-casted
    ;                                     #'(list))))
    ;                (syntax/loc stx (append e1 e2))))))
    ;        ((is-slice? type)
    ;         (is-type_ type (with-syntax* ((e1 (if (and (equal? #f primary-atom) (not primary-is-int-index))
    ;                                               primary-casted
    ;                                               #'(list)))
    ;                                       (e2 (if (and (equal? #f factor-atom) (not factor-is-int-index))
    ;                                               factor-casted
    ;                                               #'(list))))
    ;                          (syntax/loc stx (append e1 e2)))))
    ;        ((equal? type 'int)
    ;         (is-int
    ;          (syntax/loc stx '())))
    ;        ((equal? type 'int-index)
    ;         (is-int-index
    ;            (if (and primary-atom factor-atom)
    ;              (quasisyntax/loc stx #,((if (equal? op "+") fx+ fx-) primary-atom factor-atom))
    ;              (quasisyntax/loc stx (#,(if (equal? op "+") #'fx+ #'fx-) #,primary-expanded #,factor-expanded)))))
    ;        (else
    ;         (raise-syntax-error #f (format "unsupported type: ~v" type) stx)))))
    (error "^ operator is not implemented")]

   [({~literal factor} primary)
    (with-syntax* ((expanded (local-expand #'primary 'expression '()))
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
                   ((expanded (local-expand #'primary 'expression '())))
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
                 (with-syntax* ((expanded (local-expand #'primary 'expression '()))
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
              
              (datum->syntax stx
                             `(assignation ,#'name ,@getter-for-splice
                                            "="
                                            (expr
                                             (expr (term (factor (primary (element (get-value ,#'name ,@getter-for-splice)))))) ,binop ,#'value))))
             ((or "*=" "/=")
              (datum->syntax stx
                             `(assignation ,#'name ,@getter-for-splice "="
                                            (expr
                                             (term
                                              (term (factor (primary (element (get-value ,#'name ,@getter-for-splice))))) ,binop ,#'value))))
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
                   ((maybe-array-range) (cadr type)))
                   
       (with-syntax* ((index-exp (local-expand #'getter.index 'expression '()))
                      (index-start-expanded (if slice-colon_ 
                                                (if (syntax->datum #'getter.index-start)
                                                    (local-expand #'getter.index-start 'expression '())
                                                    (is-type_ 'int-index #'0))
                                                #f))
                      (index-end-expanded (if slice-colon_
                                              (if (syntax->datum #'getter.index-end)
                                                  (local-expand #'getter.index-end 'expression '())
                                                  (is-type_ 'int-index (datum->syntax stx (car maybe-array-range))))
                                              #f))
                      (value-exp_ (local-expand #'value 'expression '()))
                      (value-exp #'value-exp_)
                      (value-type (syntax-property #'value-exp_ 'landau-type)))
         ;; NOTE: Indexes are resolved at zeroth pass
         (let* ((err-msg "index must be an integer expression of constants and loop variables")
                (throw-if-not-int-index
                 (lambda (idx-stx idx-expanded-stx)
                   (unless (equal? (syntax-property idx-expanded-stx 'landau-type) 'int-index)
                     (raise-syntax-error #f err-msg idx-stx)))))
           (when (syntax->datum #'getter.index)
             (throw-if-not-int-index #'getter.index #'index-exp))
           (when (and (not slice-colon_) (is-slice? (syntax->datum #'value-type)))
             (raise-syntax-error #f "the right-hand expression is a slice, but the left-hand one is not. A slice must be assigned to a slice" stx))
           (when slice-colon_
             (when (syntax->datum #'index-start-expanded)
               (throw-if-not-int-index #'getter.index-start #'index-start-expanded))
             (when (syntax->datum #'index-end-expanded)
               (throw-if-not-int-index #'getter.index-end #'index-end-expanded))))
         (with-syntax ((name-symb (var-symbol 
                                   (symbol->string (syntax->datum #'name))
                                   src-pos)))
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
                                  (local-expand #'getter.index 'expression '()))
                                 (expanded-range (cadr (syntax->datum (local-expand (datum->syntax #'getter.index maybe-array-range) 'expression '())))))
                    ;; NOTE: (cadr (cadr ... because expanded looks like that: '(#%app '4)
                    ;; NOTE: can't do check in transmormer phase because loop variable is not resolved by that time.
                    ;;       Need to insert check-bounds code in generated code
                    (quasisyntax/loc
                        stx
                      (if (fx>= #,#'index-expanded #,#'expanded-range)
                          (raise-syntax-error
                           #f
                           (format "index ~a out of range: expect [0, ~a)" #,#'index-expanded #,#'expanded-range) #'#,stx)
                          (list '#,assign-label (list 'array-ref name-symb #,#'index-expanded) value-exp))))))
               (slice-colon_
                (let ((slice-idx-symb 'slice_idx))
                  (when (equal? maybe-array-range '())
                    (raise-syntax-error #f (format "Variable ~a is not an array." name_) #'name))
                  (with-syntax* ((slice-idx-synt (datum->syntax stx slice-idx-symb))
                                 (expanded-range (cadr (syntax->datum (local-expand (datum->syntax #'stx maybe-array-range) 'expression '()))))
                                 (index-start-expanded_ (if (atom-number #'index-start-expanded)
                                                            (atom-number #'index-start-expanded)
                                                            #'index-start-expanded))
                                 (index-end-expanded_ (if (atom-number #'index-end-expanded)
                                                          (atom-number #'index-end-expanded)
                                                          #'index-end-expanded))
                                 (lvalue-slice-range #'(fx- index-end-expanded_
                                                            index-start-expanded_))
                                 (rvalue-slice-range (get-slice-range (syntax-e #'value-type))))
                    
                    (quasisyntax/loc
                        stx
                      (begin
                        (when (fx< #,#'index-start-expanded 0) (raise-syntax-error #f (format "index must be nonnegative, given: ~v" #,#'index-start-expanded) #'#,stx))
                        (when #,#'rvalue-slice-range
                          (unless (fx= #,#'lvalue-slice-range #,#'rvalue-slice-range)
                            (raise-syntax-error #f 
                                                (format "can not cast right-hand side range ~v to the left-hand side range ~v" #,#'rvalue-slice-range #,#'lvalue-slice-range) #'#,stx)))
                        (when (fx> #,#'index-end-expanded #,#'expanded-range)
                          (raise-syntax-error
                           #f
                           (format "slice end index ~a is out of range: expect [0, ~a)" #,#'index-end-expanded #,#'expanded-range) 
                           #'#,stx))                          
                        (list 'nested
                              (for/list ((slice-idx-synt (in-range 0 #,#'lvalue-slice-range)))
                                (list '#,assign-label (list 'array-ref name-symb (fx+ #,#'index-start-expanded slice-idx-synt)) value-exp))))))))
               (else
                (quasisyntax/loc stx
                  (list '#,assign-label (list 'array-ref name-symb 0) value-exp)))))))))))
               

                 

(define-syntax (var-decl stx)
  (syntax-parse stx
    (({~literal var-decl}
      type:expr name:id (~optional (~seq "=" value:expr)
                                   #:defaults ((value #'notset))))
     (let* ((name_ (syntax->datum #'name))
            (type-datum (syntax->datum #'type))
            (type (parse-type type-datum)))
       (check-duplicate-variable-name name_ #'name)
       (add-variable! (syntax-parameter-value #'current-variables) #'name type)
       (add-real-var! #'name type stx (syntax-parameter-value #'real-vars-table))
       
       (if (equal? (syntax->datum #'value) 'notset)
           (syntax/loc stx '())
           (quasisyntax/loc stx
             #,(datum->syntax
                stx `(assignation ,#'name "=" ,#'value))))))))
                         

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
                             (local-expand arr-item 'expression '())))
                  ))
                (hash-set! constants (syntax->datum #'name)
                           (constant 'constant-name-placeholder type 'constant-array-placeholder))
                (syntax/loc stx '()))))
             (else
              (let* ((expanded (local-expand #'value 'expression '()))
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
     (let* ((expanded-size (local-expand #'size 'expression '())))
              
       (hash-set! parameters (syntax->datum #'name)
                  expanded-size)
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
                        [start-val (local-expand #'start 'expression '())]
                        [stop-val (local-expand #'stop 'expression '())])
            (throw-if-not-type 'int-index #'start-val err-msg)
            (throw-if-not-type 'int-index #'stop-val err-msg)
            (quasisyntax/loc stx
              (begin
                ;(writeln (format "stop-val: ~v" (syntax->datum #'stop-val)))
                (syntax-parameterize
                    [(current-variables
                      '#,new-vars-layer)]
                  (list 'nested
                        (for/list ([symm (in-range start-val stop-val)])
                          pat.body))))))))]))
    

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
                    (val1 (local-expand #'n1 'expression '()))
                    (val2 (local-expand #'n2 'expression '())))
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
         [#f (list 'nested (list body-false ...))]
         [#t (list 'nested (list body-true ...))]
         [else (raise-syntax-error #f "not a boolean in if-expr")]))]
         
    [({~literal if-expr} "if" b-expr "{" body ... "}")
     (syntax/loc stx
       (match b-expr
         [#f '()]
         [#t (list 'nested (list body ...))]
         [else (raise-syntax-error #f "not a boolean in if-expr")]))]))
         
;; FIXME: check variable shadowing other functions names
(define-for-syntax (check-duplicate-variable-name name stx-name)
  (when (search-variable name (syntax-parameter-value #'current-variables))
    (raise-syntax-error #f "duplicate variable declaration" stx-name))
  (when (hash-has-key? constants name)
    (raise-syntax-error #f "variable name shadows constant" stx-name))
  (when (hash-has-key? parameters name)
    (raise-syntax-error #f "variable name shadows parameter" stx-name))
  (when (hash-has-key? (syntax-parameter-value #'current-arguments) name)
    (raise-syntax-error #f "variable name shadows argument" stx-name))
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


(define-for-syntax (check-real-func-existence name [get-range #f])
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
                     (raise-syntax-error #f "name not found" name))))))))))
                     
       (type (car type_))
       (maybe-range (cadr type_)))
       
    (unless (equal? type 'real)
      (raise-syntax-error
       #f
       err-msg #'name))
    (if get-range
        (if (equal? maybe-range '()) #f maybe-range)
        type_)))
  

(define-for-syntax (check-real-arg-or-parameter-existence name [get-range #f])
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
                     (raise-syntax-error #f "name not found" name))))))))))
                     
       (type (car type_))
       (maybe-range (cadr type_)))
       
    (unless (equal? type 'real)
      (raise-syntax-error
       #f
       "only real function arguments are allowed be inside the derivative anotation or application" name))
    (if get-range
        (if (equal? maybe-range '()) #f maybe-range)
        type_)))
  

(define-for-syntax (check-real-var-existence name [get-range #f])
  (let*
      ((name_ (syntax->datum name))
       (type_
        (cond
          ((hash-has-key? constants name_)
           (raise-syntax-error
            #f
            "constant can not be inside the derivative annotation or discard. Expect a real argument or a variable." name))
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
                     (raise-syntax-error #f "name not found" name))))))))))
                     
       (maybe-range (cadr type_))
       (type (car type_)))
    (unless (equal? type 'real)
      (raise-syntax-error
       #f
       "only real arguments or variables are allowed be inside the derivative discard" #'name))
    (if get-range
        (if (equal? maybe-range '()) #f maybe-range)
        type_)))
  
(define-syntax (get-derivative stx)
  (syntax-parse stx
    ((_ get-value-1 "'" get-value-2)
     (syntax/loc stx (list 'get-derivative get-value-1 get-value-1)))))
;; NOTE:
;; non-array dx arguments are casted to single-cell array
