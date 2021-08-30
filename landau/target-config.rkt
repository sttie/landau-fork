#lang racket
(require racket/base
         racket/contract
         json)

(provide target/c lang/c TARGET target-lang target-extfloat? target-debug? target-c-out-path)

(struct target
    (lang
     extfloat?
     debug?
     c-out-path)
    #:prefab)

(define lang/c (one-of/c 'racket 'ansi-c))
(define target/c (struct/c target 
                           lang/c 
                           boolean?
                           boolean?
                           path-for-some-system?))

(define (get-source-name)
  (let ((full-name-lst (string-split (format "~a" (find-system-path 'run-file)) ".")))
    (when (empty? full-name-lst)
      (error "Empty Landau source name"))
    (cond
      ;; FIXME: when .dau is called with drracket full-name-lst is drracket's path
      ((absolute-path? (car full-name-lst)) "landau_output")
      (else (file-name-from-path (car full-name-lst))))))

(define (path-only_ path)
  (define dir (path-only path))
  (if dir dir "./"))

(define (resolve-path relative-path)
  (let* ((source-name (get-source-name))
         (target-name (format "~a.c" source-name))
         (out-path relative-path)
         (resolved-path 
          (cond
            ((absolute-path? out-path)
             (build-path out-path))
            (else
             (let*
              ((root-path (path-only_ (find-system-path 'run-file)))
               (out-dir-resolved (build-path root-path relative-path)))
              (simplify-path out-dir-resolved))))))
        
        (unless (directory-exists? resolved-path)
                (error (format "output_directory: directory: ~a does not exist" resolved-path)))
        (build-path resolved-path target-name)))

(define (parse-target-cfg inp-path)
  (let* ((inp (open-input-file inp-path))
         (parsed (read-json inp))
         (lang (hash-ref parsed 'target_language))
         (extfloat (hash-ref parsed 'use_extfloat))
         (path (hash-ref parsed 'output_directory)))
        (close-input-port inp)
        (unless (string? lang)
                (raise-type-error 'target_language "string" lang))
        (unless (lang/c (string->symbol lang))
                (error (format
                        "target_language should be one of: \"racket\", \"ansi-c\". Given: \"~a\""
                        (string->symbol lang))))

        (unless (boolean? extfloat)
                (raise-type-error 'use_extfloat "boolean" extfloat))

        (unless (string? path)
                (raise-type-error 'output_directory "path" path))

        (define c-output-resolved-path (resolve-path path))

        (target (string->symbol lang)
                extfloat
                #f
                c-output-resolved-path)))

(define/contract (override-with-command-line cfg-target)
  (-> target/c
      target/c)
  (let* ((args (current-command-line-arguments))
         (config-target-lang (target-lang cfg-target))
         (config-target-extfloat (target-extfloat? cfg-target))
         ;; FIXME: do not affect if .dau is called not from ./
         (final-target-lang (if (vector-member "-c" args)
                              'ansi-c
                              config-target-lang))
         (final-target-extfl (if (vector-member "-extfl" args)
                               #t
                                config-target-extfloat))
         (final-target-debug? (if (vector-member "-debug" args)
                                #t
                                #f)))
        (target final-target-lang
                final-target-extfl
                final-target-debug?
                (target-c-out-path cfg-target))))

(define/contract
  TARGET
  target/c
  (override-with-command-line
   (parse-target-cfg 
    (collection-file-path "../config.json" "landau"))))
