#lang debug racket/base
;; TODO move the module to the syntax-analyser directory
(require 
  racket/base
  racket/contract
  racket/contract/region
  racket/fixnum
  racket/flonum
  racket/function
  racket/list
  racket/match
  racket/pretty
  racket/set
  racket/stxparam

  "environment.rkt"
  "type-utils.rkt"
  "common-for-syntax.rkt")

(provide process! splice-nested)


(define debug-backrun #t)

(define (log-debug msg)
  (when debug-backrun (pretty-print msg)))

(define (coerce-to-list any)
  (if (list? any)
    any
    (list)))

(define (splice-step x acc)
  (cond 
    ((list? x) (cond
                 ((empty? x) acc)
                 ((symbol? (first x))
                  (cons x acc))
                 (else
                   (append (splice-nested x) acc))))
    (else (error "not a list"))))

(define (splice-nested lst)
  (foldl splice-step '() lst))
      
(define (get-idx df)
  (if (equal? (car df) 'array-ref)
    (caddr df)
    #f))

(define (add-lval-ders-to-rval! der-table available-dx-table discard-table l-val-ref r-val-ref)

  (let* ((l-val-df-table (hash-ref! der-table (ref-to-key l-val-ref) (make-hash)))
         (r-val-df-table (hash-ref! der-table (ref-to-key r-val-ref) (make-hash)))
         (r-val-available-dx-table (if (hash-has-key? available-dx-table (ref-to-key r-val-ref))
                                     (hash-ref available-dx-table (ref-to-key r-val-ref))
                                     #f))
         (r-val-discard-table
           (if (hash-has-key? discard-table (ref-to-key r-val-ref))
             (hash-ref discard-table (ref-to-key r-val-ref))
             #f)))


    (when r-val-available-dx-table
      (for ((l-val-dx-name (in-list (hash-keys l-val-df-table))))
        ;; NOTE: If r-val has no such dx-name then skip it
        (when (hash-has-key? r-val-available-dx-table l-val-dx-name)
          (let* ((l-val-der-bundle (hash-ref l-val-df-table l-val-dx-name))
                 (l-val-der-bundle-type (car l-val-der-bundle))
                 (r-val-available-dx-idx (cadr (hash-ref r-val-available-dx-table l-val-dx-name)))
                 (r-val-available-dx-idx-immutable (list->set (set->list r-val-available-dx-idx)))
                 ;; maybe-dx-idxs :: Either (Set Int) (Box Bool)
                 (l-val-maybe-dx-idxs (cadr l-val-der-bundle))
                 (l-val-dx-idx-immutable (list->set (set->list l-val-maybe-dx-idxs)))
                 (r-val-der-discard-bundle
                   (cond
                     ((equal? r-val-discard-table #f) #f)
                     ((hash-has-key? r-val-discard-table l-val-dx-name)
                      (hash-ref r-val-discard-table l-val-dx-name))
                     (else #f))))

            (if (hash-has-key? r-val-df-table l-val-dx-name)
              (let* ((r-val-der-bundle (hash-ref r-val-df-table l-val-dx-name))
                     (r-val-der-bundle-type (car r-val-der-bundle))
                     (dx-idxs-to-discard (if r-val-der-discard-bundle (cadr r-val-der-discard-bundle) #f))
                     (r-val-maybe-dx-idxs (cadr r-val-der-bundle)))

                (cond
                  ((and (equal? l-val-der-bundle-type 'dx-idxs)
                        (equal? r-val-der-bundle-type 'dx-idxs))
                   ;; NOTE: if r-val already has such key, than union indexes

                   (begin
                     (set-union! r-val-maybe-dx-idxs
                                 (set-intersect
                                   (if dx-idxs-to-discard
                                     (set-subtract
                                       l-val-dx-idx-immutable
                                       (list->set (set->list dx-idxs-to-discard)))
                                     l-val-dx-idx-immutable)
                                   r-val-available-dx-idx-immutable))

                     (when (set-empty? r-val-maybe-dx-idxs)
                       (hash-remove! r-val-df-table l-val-dx-name))))
                  ((and (equal? l-val-der-bundle-type 'dx)
                        (equal? r-val-der-bundle-type 'dx))
                   ;; NOTE: Some code been there in 8110d17e717930726a296ba86c09873f9d1eb04f
                   (error (format "bug: (equal? l-val-der-bundle-type 'dx-idxs) failed. dx is not an array")))))

              ;; NOTE: if r-val had no such der-bundle just copy it
              ;; and then discard indexes if need
              (if r-val-der-discard-bundle
                ;; NOTE: subtract indexes that shoud be discarded
                (let* ((dx-idxs-to-discard (cadr r-val-der-discard-bundle))
                       (l-val-der-bundle-type (car l-val-der-bundle))

                       (idxs-after-discard (set-intersect
                                             r-val-available-dx-idx-immutable
                                             (set-subtract
                                               l-val-dx-idx-immutable
                                               (list->set (set->list dx-idxs-to-discard))))))

                  (unless (set-empty? idxs-after-discard)
                    (let ((l-val-der-bundle-discarded (list l-val-der-bundle-type
                                                            idxs-after-discard)))
                      (hash-set! r-val-df-table l-val-dx-name l-val-der-bundle-discarded))))

                (begin
                  (let* ((intersected-immutable-idxs-set (set-intersect r-val-available-dx-idx-immutable
                                                                        l-val-dx-idx-immutable))
                         (intersected-mutable-idxs-set
                           (for/mutable-set ((ix (in-set intersected-immutable-idxs-set)))
                                            ix)))
                    (unless (set-empty? intersected-immutable-idxs-set)
                      (hash-set! r-val-df-table l-val-dx-name
                                 (list l-val-der-bundle-type
                                       intersected-mutable-idxs-set)))))))))))))


;; NOTE: Traverse actions-list in forward direction and populate available-dx-table and deps-set
(define/contract
  (process-forward! forward-flat available-dx-table deps-set real-vars-table)
  (-> (listof any/c) any/c any/c real-vars-table/c any/c)
  (for ((action (in-list forward-flat)))
          (match action

            ((list 'var-decl name-vs type)
             (match type
               ((list 'real (list))
                (hash-set! real-vars-table name-vs #f))

               ((list 'real (list array-len))
                (hash-set! real-vars-table name-vs array-len))
               (_ (void))))

            ((list 'der-annot df-ref dx-ref)
             (let*-values
                 (((dx-key) (format "~a" (cadr dx-ref)))
                  ((df-idx) (get-idx df-ref))
                  ((dx-idx) (get-idx dx-ref))
                  ; ((i) (println (format "~a ~a" df-ref dx-ref)))
                  ;; NOTE: key for deps-set

                  ((df-key) (ref-to-key df-ref))
                  ((df-dx-table) (hash-ref! available-dx-table df-key (make-hash)))
                  ((derivatives-bundle)
                   (hash-ref! df-dx-table
                              dx-key
                              (if dx-idx
                                  ;; NOTE: if dx is array
                                  (list 'dx-idxs (mutable-set))
                                  ;; NOTE: if dx is basic variable
                                  (list 'dx (box #f)))))
                              
                  ((bundle-type maybe-dx-indexes)
                   (values (car derivatives-bundle) (cadr derivatives-bundle))))
               (match bundle-type
                 ('dx-idxs (set-add! maybe-dx-indexes dx-idx))
                 (error (format "bug: expect 'dx-idxs")))
               (set-add! deps-set df-key)))
               
            ((or (list 'assign l-val-ref refs-list) (list 'func-assign l-val-ref refs-list))
             (let* ((l-val-key (ref-to-key l-val-ref))
                    ;; NOTE: table with available derivatives. Exist only if l-val-ref been used in der-annot
                    (l-val-dx-table (hash-ref! available-dx-table (ref-to-key l-val-ref) (make-hash))))
               ;; NOTE: coerce-to-list used to ignore 'int and 'int-index values in the rhs 
               (for ((r-val-ref (in-list (coerce-to-list refs-list))))
                 (let ((r-val-dx-table (if (hash-has-key? available-dx-table (ref-to-key r-val-ref))
                                           (let
                                             ((tbl (hash-ref available-dx-table (ref-to-key r-val-ref))))
                                             tbl)
                                           #f)))
                   (when r-val-dx-table
                     (let ((r-val-dx-names (hash-keys r-val-dx-table)))
                       (for ((dx-name (in-list r-val-dx-names)))
                         (let ((r-val-dx-idxs (cadr (hash-ref r-val-dx-table dx-name)))
                               (l-val-dx-idxs (cadr (hash-ref! l-val-dx-table dx-name (list 'dx-idxs (mutable-set))))))
                           (set-union! l-val-dx-idxs r-val-dx-idxs)))))))))
            (_ (void)))))

;; NOTE: Traverse reverse actions list and populate deps-set der-table
;; available-dx-table and discard-table are read only
(define/contract
  (process-reverse! rev-flat available-dx-table deps-set der-table discard-table)
  (-> (listof any/c) any/c any/c any/c any/c any/c)
  (for ((action (in-list rev-flat)))
    (match action
      ((list 'var-decl _ _)
       ;; NOTE: Alreary done in process-forward!
       (void))

      ((list 'der-apply df dx)
       (let*-values
         (((dx-key) (format "~a" (cadr dx)))

          ((df-idx) (get-idx df))
          ((dx-idx) (get-idx dx))
          ;; NOTE: key for deps-set
          ((df-key) (ref-to-key df))
          ((df-derivatives-table) (hash-ref! der-table df-key (make-hash)))

          ((derivatives-bundle)
           (hash-ref! df-derivatives-table
                      dx-key
                      (if dx-idx
                        ;; NOTE: if dx is array
                        (list 'dx-idxs (mutable-set))
                        ;; NOTE: if dx is basic variable
                        (list 'dx (box #f)))))

          ((bundle-type maybe-dx-indexes)
           (values (car derivatives-bundle) (cadr derivatives-bundle))))
         (match bundle-type
           ('dx-idxs (set-add! maybe-dx-indexes dx-idx))
           ('dx (set-box! maybe-dx-indexes #t)))
         (set-add! deps-set df-key)))

      ((list 'discard df dx)
       (let*-values
         (((dx-key) (format "~a" (cadr dx)))

          ((df-idx) (get-idx df))
          ((dx-idx) (get-idx dx))
          ;; NOTE: key for deps-set
          ((df-key) (ref-to-key df))
          ;; NOTE: get discard-derivatives bundle
          ((df-derivatives-discard-table) (hash-ref! discard-table df-key (make-hash)))
          ((derivatives-discard-bundle)
           (hash-ref! df-derivatives-discard-table
                      dx-key
                      (if dx-idx
                        ;; NOTE: if dx is array
                        (list 'dx-idxs (mutable-set))
                        ;; NOTE: if dx is basic variable
                        (list 'dx (box #f)))))

          ((discard-bundle-type maybe-dx-indexes-to-discard)
           (values (car derivatives-discard-bundle) (cadr derivatives-discard-bundle)))
          ;; NOTE: get der-bundle
          ((df-derivatives-table)
           (cond
             ((hash-has-key? der-table df-key) (hash-ref der-table df-key))
             (else #f)))

          ((derivatives-bundle)
           (cond
             ((and df-derivatives-table
                   (hash-has-key? df-derivatives-table dx-key))
              (hash-ref df-derivatives-table dx-key))
             (else #f)))
          ((bundle-type maybe-dx-indexes)
           (cond
             (derivatives-bundle
               (values (car derivatives-bundle) (cadr derivatives-bundle)))
             (else
               (values #f #f)))))
         (match bundle-type
           ('dx-idxs
            (cond
              ((set-member? maybe-dx-indexes dx-idx)
               (begin
                 (set-remove! maybe-dx-indexes dx-idx)
                 (when (set-empty? maybe-dx-indexes) (hash-remove! df-derivatives-table dx-key))))

              ;; NOTE: idxs set haven't filled with discarded derivative (yet, or it is just not used)
              (else #f)))

           ('dx (hash-remove! df-derivatives-table dx-key))
           (else #f))
         ;; NOTE: fill discard table
         (match discard-bundle-type
           ('dx-idxs (set-add! maybe-dx-indexes-to-discard dx-idx))
           ('dx (set-box! maybe-dx-indexes-to-discard #t)))))

      ;; FIXME: maybe need to use discard table
      ((list 'func-assign l-val-ref refs-list)
       (let ((f-key (ref-to-key l-val-ref)))

         (set-add! deps-set f-key)
         ;; NOTE: coerce-to-list used to ignore 'int and 'int-index values in the rhs                
         (for ((ref (in-list (coerce-to-list refs-list))))
           (let ((r-val-key (ref-to-key ref)))
             (set-add! deps-set r-val-key)
             ;; NOTE: Add derivatives of reference (e.g. f[0] or just x) to
             ;;       l-val derivatives
             (add-lval-ders-to-rval!
               der-table
               available-dx-table
               discard-table
               l-val-ref
               ref)
             ))))

      ((list 'assign l-val-ref refs-list)
       (let ((l-val-key (ref-to-key l-val-ref)))
         (when (set-member? deps-set l-val-key)
           (set-remove! deps-set l-val-key)
           ;; NOTE: coerce-to-list used to ignore 'int and 'int-index values in the rhs
           (for ((ref (in-list (coerce-to-list refs-list))))
             (let ((r-val-key (ref-to-key ref)))
               (set-add! deps-set r-val-key)
               ;; NOTE: Add derivatives of reference (e.g. f[0] or just x) to
               ;;       l-val derivatives
               (add-lval-ders-to-rval!
                 der-table
                 available-dx-table
                 discard-table
                 l-val-ref
                 ref)
               ;  (displayln (format "~a <- ~a" l-val-ref ref))
               ;  (pretty-print der-table)
               )))))

      ((list 'der-annot ref-1 ref-2)
       (void))

      ;; TODO: 
      ;; t[0] -> (hash "x[0]" (list 'dx-idxs (mutable-set 0)))
      (_ (void)))))


;; NOTE: process actions-list (function actions)
(define/contract 
  (process! actions-list dx-names-set real-vars-table)
  (-> (listof any/c) (hash/c string? boolean?) real-vars-table/c 
      (values
       der-table/c
       dx-names-set/c
       real-vars-table/c))

  (begin
    (let ((rev-flat (splice-nested actions-list)))
      (let*
        ((der-table (make-hash))
         ;; NOTE: discard-table :: (var-ame, IntIndex) -> (table dx-name -> Set dx-indexes-to-discard)
         (discard-table (make-hash))
         (available-dx-table (make-hash))
         (forward-flat (reverse rev-flat))
         (deps-set (mutable-set)))
        ;; NOTE: populate table with function code for the future inlining, when this function is called
        (process-forward! forward-flat available-dx-table deps-set real-vars-table)
        (process-reverse! rev-flat available-dx-table deps-set der-table discard-table)
        (values der-table dx-names-set real-vars-table)))))


