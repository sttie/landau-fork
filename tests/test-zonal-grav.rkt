#lang racket/base

(require racket/flonum math/flonum racket/extflonum)
(require "zonal_grav.dau")

(define jupiter-grav
  (flvector 0.0 0.0 -0.006572081058301649 0.0 0.00019571333333323615 0.0 -9.499240999993885e-06))

(zonal_grav .1 .2 .3 jupiter-grav)
