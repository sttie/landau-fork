#lang racket
#| INFO: Routines for generating direct and inverse mappings |#
(require 
  (only-in "common-for-syntax.rkt" vec->fxvec fxvector->vector define/contract-for-syntax)
  (only-in "combinators.rkt" c-define-array)
  (for-syntax racket/syntax
              racket/base
              racket/vector
              racket/match
              racket/set
              racket/contract
              "environment.rkt"
              "target-config.rkt"
              "combinators.rkt"
              "metalang.rkt"
              (only-in "common-for-syntax.rkt" fxvec->vec)
              racket/stxparam
              racket/list
              racket/flonum
              racket/extflonum
              racket/fixnum))

(provide (for-syntax make-mappings-decl-list-helper!))

(define-for-syntax debug #f)

(define/contract-for-syntax
  ;; NOTE: generate mappings for basic variable (not array)
  ;; return mappings (straigt and inverse) vectors and their symbols
  ;; WARNING: mutate need-only-value-set
  ;; WARNING: mutate current-variables 
  ;; WARNING: mutate need-derivatives-table
  (generate-mappings-for-basic-variable! derivatives-info
                                         current-variables
                                         var-name-key 
                                         dx-name-str )
  (-> derivatives-info/c current-variables/c var-symbol/c dx-name/c 
      (values
        (or/c symbol? #f)
        (or/c fxvector? #f)
        (or/c symbol? #f)
        (or/c fxvector? #f)))
  ;; NOTE: ditry trick with passing ref to ref-to-key. der-table include only (list/c df-name/c integer?) keys
  (let* ((der-table (derivatives-info-.der-table derivatives-info))
         (need-only-value-set (derivatives-info-.need-only-value-set derivatives-info))
         (need-derivatives-table (derivatives-info-.need-derivatives-table derivatives-info))
         (df-table (dtbl-get-df-table 
                     der-table
                     (ref-to-key (list 'array-ref var-name-key 0)))))
    (cond
      ((hash-empty? df-table)
       (set-add! need-only-value-set var-name-key)
       (values #f #f #f #f))

      (else
        (cond
          ((hash-has-key? df-table dx-name-str)
           (let* ((der-bundle (hash-ref df-table dx-name-str))
                  (dx-indexes-sorted
                    (sort (if (equal? 'dx-idxs (car der-bundle))
                            (set->list (cadr der-bundle))
                            (error "bug: car is not 'dx-idxs")) <))
                  (mapped-dx-size (length dx-indexes-sorted))
                  (mapping (for/fxvector ((i (in-list dx-indexes-sorted))) i))
                  (inv-mapping-period-value (fx+ (apply fxmax dx-indexes-sorted) 1))
                  (dx-idx-mapping (make-fxvector inv-mapping-period-value -1)))

             ;; TODO: refactor
             ;(hash-set! mapped-dx-size-GLOBAL (make-mapped-dx-size-key var-name-key dx-name-str) mapped-dx-size)
             (unless (var-symbol/c var-name-key)
               (error "71")) 
             (if (ndt-member? need-derivatives-table var-name-key)
               (let ((dx-names-sizes (hash-ref need-derivatives-table var-name-key)))
                 (hash-set! dx-names-sizes dx-name-str (mapping-sizes mapped-dx-size inv-mapping-period-value)))
               (hash-set!
                 need-derivatives-table
                 var-name-key
                 (make-hash (list (cons dx-name-str (mapping-sizes mapped-dx-size inv-mapping-period-value))))))

             (for ((dx-idx (in-list dx-indexes-sorted))
                   (mapped-dx-index (in-naturals)))
               (fxvector-set! dx-idx-mapping dx-idx mapped-dx-index))

             (values

               (add-variable!
                 current-variables
                 (datum->syntax #f (make-var-mappings-name var-name-key dx-name-str))
                 'mappings)

               mapping

               (add-variable!
                 current-variables
                 (datum->syntax #f (make-dx-idxs-mappings-name var-name-key dx-name-str))
                 'dx-mappings)

               dx-idx-mapping)))

          (else (values #f #f #f #f)))))))


;; NOTE: WARNING: mutation derivatives-info
(define/contract-for-syntax
  (generate-mappings-for-array-variable! derivatives-info
                                         current-variables 
                                         var-name-key
                                         dx-name-str 
                                         n-mappings
                                         grouped-keys-table)
  (-> derivatives-info/c
      current-variables/c
      var-symbol/c
      dx-name/c
      fixnum?
      grouped-keys-table/c 
      (values 
        (or/c symbol? #f)
        (or/c fxvector? #f)
        (or/c symbol? #f)
        (or/c fxvector? #f)))
  (let* ((der-table (derivatives-info-.der-table derivatives-info))
         (need-only-value-set (derivatives-info-.need-only-value-set derivatives-info))
         (need-derivatives-table (derivatives-info-.need-derivatives-table derivatives-info))
         (df-tables-are-empty
           (for/and ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
             (let ((df-table (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))
               (hash-empty? df-table)))))

    (if df-tables-are-empty
      (begin
        (set-add! need-only-value-set var-name-key)
        (values #f #f #f #f))
      ;; FIXME: at some indexses (hash-ref df-table dx-name-str) can hae no value
      (let ((df-tables-has-no-dx-key
              (for/and ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
                (let ((df-table (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))
                  (not (hash-has-key? df-table dx-name-str))))))


        (if df-tables-has-no-dx-key
          (values #f #f #f #f)
          (let* ((needed-dx-idx-count
                   (for/list ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
                     (let
                       ((df-table (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))
                       (cond
                         ((hash-empty? df-table) 0)
                         ((hash-has-key? df-table dx-name-str)
                          (let* ((der-bundle (hash-ref df-table dx-name-str))
                                 (dx-indexes
                                   (if (equal? 'dx-idxs (car der-bundle))
                                     (set->list (cadr der-bundle))
                                     (error "bug: car is not 'dx-idxs"))))
                            (length dx-indexes)))
                         (else 0)))))
                 ;(i (writeln (format "needed-dfdx-idxs-count: ~a" needed-dx-idx-count)))
                 (mapped-dx-size
                   (apply fxmax needed-dx-idx-count))
                 ;; NOTE: df-mappings is a vector size of df-array * dx-period. Left-aligned, padded right with -1
                 ;; df-mappings :: Vector (U (List dx-idx) Bool)
                 (df-mappings (make-fxvector (fx* n-mappings mapped-dx-size) -1))

                 (inv-mapping-period-value
                   (fx+ 1 (apply fxmax
                                 (for/list ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
                                   (let* ((df-table
                                            (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))
                                     (if (hash-has-key? df-table dx-name-str)
                                       (let* ((der-bundle (hash-ref df-table dx-name-str)))

                                         (apply fxmax
                                                (if (equal? 'dx-idxs (car der-bundle))
                                                  (set->list (cadr der-bundle))
                                                  (error "bug: car is not 'dx-idxs"))))
                                       0))))))

                 (dx-idx-mappings (make-fxvector (fx* n-mappings inv-mapping-period-value) -1)))
            ;(hash-set! mapped-dx-size-GLOBAL (make-mapped-dx-size-key var-name-key dx-name-str) mapped-dx-size)

            (if (hash-has-key? need-derivatives-table var-name-key)
              (let ((dx-names-sizes (hash-ref need-derivatives-table var-name-key)))

                (hash-set! dx-names-sizes dx-name-str (mapping-sizes mapped-dx-size inv-mapping-period-value)))
              (hash-set!
                need-derivatives-table
                var-name-key
                (make-hash (list (cons dx-name-str (mapping-sizes mapped-dx-size inv-mapping-period-value))))))

            (for ((var-array-idx (in-list (hash-ref grouped-keys-table var-name-key))))
              (let* ((df-table
                       (hash-ref der-table (ref-to-key (list 'array-ref var-name-key var-array-idx)))))

                (when (hash-has-key? df-table dx-name-str)
                  (let* ((der-bundle (hash-ref df-table dx-name-str))
                         (dx-indexes-sorted
                           ;; FIXME: try to cast to vector and sort it
                           (sort (if (equal? 'dx-idxs (car der-bundle))
                                   (set->list (cadr der-bundle))
                                   (error "bug: car is not 'dx-idxs")) <)))
                    ;(dx-idx-mapping (make-fxvector (fx+ (apply fxmax dx-indexes-sorted) 1) -1))
                    (for ((dx-idx (in-list dx-indexes-sorted))
                          (mapped-dx-index (in-naturals)))
                      (fxvector-set! dx-idx-mappings
                                     (fx+ (fx* var-array-idx inv-mapping-period-value) dx-idx) mapped-dx-index))
                    (for ((dx-idx (in-list dx-indexes-sorted))
                          (i (in-naturals)))
                      (fxvector-set! df-mappings (fx+ (fx* mapped-dx-size var-array-idx) i) dx-idx))))))

            (values
              (add-variable!
                current-variables
                (datum->syntax #f 
                               (make-var-mappings-name var-name-key dx-name-str))
                'mappings)
              df-mappings
              (add-variable!
                current-variables
                (datum->syntax #f 
                               (make-dx-idxs-mappings-name var-name-key dx-name-str))
                'dx-mappings)
              dx-idx-mappings)))))))

(define-for-syntax (fxvector-member val vec )
  (for/or ((vec-val (in-fxvector vec)))
    (equal? val vec-val)))

(define-for-syntax (include-minus-one? vec)
  (match (fxvector-member -1 vec)
    (#f #f)
    (_ #t)))

;; NOTE: Traverse the keys of all needed variables names and allocate mappings vectors if variable need derivatives
;;       Mappings and dx-idx-mappings are allocated for function return value and function arguments
;;       if thee belong to need-derivatives-set set         
;;       (e.g. df-table is empty) df-table :: Hash (String-key, (idx-type (Either (Set Int) Bool))
;; WARNING: Mutation of need-derivatives-set need-only-value-set mappings-table
(define/contract-for-syntax
  (make-mappings-decl-list-helper!
    dx-name-str
    grouped-keys-table
    real-vars-table
    der-table
    need-only-value-set
    need-derivatives-table
    current-variables
    mappings-table
    lang
    stx)
  (-> string?
      grouped-keys-table/c
      real-vars-table/c
      der-table/c
      need-only-value-set/c
      need-derivatives-table/c
      current-variables/c
      mappings-table/c
      lang/c
      (syntax/c any/c)

      (listof any/c))
  (begin
    ;  (println "grouped-keys-table")
    ;  (pretty-print grouped-keys-table)

    ;  (println "real-vars-table")
    ;  (pretty-print real-vars-table)

    ;  (println "need-derivatives-table")
    ;  (pretty-print need-derivatives-table)
    (with-syntax
      ((decl-list
         (for/list ((var-name-key (in-list (hash-keys grouped-keys-table))))
           (let* ((var-size (rvt-get-var-size real-vars-table var-name-key))
                  (n-mappings (if var-size var-size 1))
                  (is-basic-variable (not (rvt-var-is-array real-vars-table var-name-key)))
                  (derivatives-info-struct (derivatives-info der-table
                                                             need-only-value-set
                                                             need-derivatives-table)))
             (let-values
               (((mapping-var-symbol df-mappings-vec dx-idx-mappings-var-symbol dx-idx-mappings-vec)
                 (cond
                   ;; TODO: add flags to not generate inverse mappings
                   (is-basic-variable
                     ;; NOTE: non-array variable
                     (generate-mappings-for-basic-variable! derivatives-info-struct
                                                            current-variables  
                                                            var-name-key
                                                            dx-name-str))

                   ;; NOTE: array variable
                   (else
                     (generate-mappings-for-array-variable! derivatives-info-struct
                                                            current-variables 
                                                            var-name-key
                                                            dx-name-str 
                                                            n-mappings
                                                            grouped-keys-table)))))
               (when mapping-var-symbol
                (set-mappings! mappings-table dx-idx-mappings-var-symbol (if df-mappings-vec
                                                                          (mappings-info (include-minus-one? df-mappings-vec))
                                                                          #f)))
               
               (define mapping-is-needed? (and mapping-var-symbol dx-idx-mappings-var-symbol))
               (match lang
                 ;; NOTE: Racket lang
                 ('racket
                  (if mapping-is-needed?
                    (begin
                      (with-syntax* 
                        ((dx-idx-mappings-vec-syntax dx-idx-mappings-vec)
                         (mapping-var-symbol (datum->syntax stx mapping-var-symbol))
                         (dx-idx-mappings-var-symbol (datum->syntax stx dx-idx-mappings-var-symbol))
                         (dx-name-str dx-name-str))
                        (with-syntax
                          ((df-mappings-vec (fxvec->vec df-mappings-vec))
                           (dx-idx-mappings-vec (fxvec->vec dx-idx-mappings-vec)))
                          (quasisyntax/loc
                            stx
                            (begin (define mapping-var-symbol (vec->fxvec df-mappings-vec))
                                   (define dx-idx-mappings-var-symbol (vec->fxvec dx-idx-mappings-vec))
                                   )))))

                    #'(void)))

                 ;; TODO USE metalang
                 ;; NOTE: ANSI-C lang
                 ('ansi-c
                  (if mapping-is-needed? 
                    (with-syntax ((df-mappings-vec (vector->carray (fxvector->vector df-mappings-vec)))
                                  (df-mappings-vec-size (fxvector-length df-mappings-vec))
                                  (dx-idx-mappings-vec (vector->carray (fxvector->vector dx-idx-mappings-vec)))
                                  (dx-idx-mappings-vec-size (fxvector-length dx-idx-mappings-vec))
                                  (mapping-var-symbol (symbol->string mapping-var-symbol))
                                  (dx-idx-mappings-var-symbol (symbol->string dx-idx-mappings-var-symbol)))
                      (quasisyntax/loc stx
                                       (list
                                         (c-define-array "int" mapping-var-symbol
                                                         df-mappings-vec-size 
                                                         #,#'df-mappings-vec "const static")
                                         (c-define-array "int" dx-idx-mappings-var-symbol
                                                         dx-idx-mappings-vec-size 
                                                         #,#'dx-idx-mappings-vec "const static"))))

                    #'""))))))))

      (syntax-e #'decl-list))))
