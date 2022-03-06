;;;; Copyright (c) Alan M. Watson 2016.

(define-library (grey-atmosphere)

  (export q-infinity q H H_nu H_alpha_hat J_alpha B_alpha_hat theta)

  (import (scheme base))
  (import (scheme inexact))
  (import (mathematics))
  (import (integral))
  (import (exponential-integral))
  (import (astrophysics))

  (begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (define (Z u)
      (+ (expt (- 1 (* 0.5 u (log (/ (+ 1 u) (- 1 u))))) 2)
         (expt (* 0.5 pi u) 2)))

    (define (H mu)
      (/ (exp (/ (integral 0.0 (* 0.5 pi)
                           (lambda (theta)
                             (/ (* theta (atan (* mu (tan theta))))
                                (- 1 (/ theta (tan theta))))))
                 pi))
         (sqrt (+ mu 1))))

    ;;(display (H 0.0))
    ;;(newline)
    ;;(display (H 1.0))
    ;;(newline)

    (define q-infinity (/ (integral 0.0 1.0 (lambda (mu) (* (H mu) mu mu)))
                          (integral 0.0 1.0 (lambda (mu) (* (H mu) mu)))))

    ;;(display q-infinity)
    ;;(newline)

    (define (q tau)
      ;; The Hofp function. See Mihalas (1978, pp. 71-73)
      (- q-infinity
         (/ (integral 0.0 1.0
                      (lambda (u)
                        (/ (exp (- (/ tau u))) (H u) (Z u))))
            (* 2 (sqrt 3)))))

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    (define (theta tau)
      (expt (* 0.75 (+ tau (q tau))) 0.25))
      
    (define (J_alpha alpha tau)
      (if (= alpha 0.0)
        0.0
        (let ((integrand (lambda (tau-prime)
                           (/ (E_1 (abs (- tau-prime tau)))
                              (- (exp (/ alpha (theta tau-prime))) 1)))))
          (* 0.25 0.3079
             (+ (integral (+ tau 2)          (+ tau 4)        integrand)
                (integral tau                (+ tau 2)         integrand)
                (integral (max 0  (- tau 2)) tau               integrand)
                (integral 0                  (max 0 (- tau 4)) integrand))
             (expt alpha 3)
             (expt (theta tau) -4)))))

    (define (H_alpha_hat alpha tau)
      (if (= alpha 0.0)
        0.0
        (let ((integrand (lambda (tau-prime)
                           (/ (E_2 (abs (- tau-prime tau)))
                              (- (exp (/ alpha (theta tau-prime))) 1)))))
          (* (/ 30 (expt pi 4))
             (expt alpha 3)
             (- (integral tau (+ tau 20) integrand)
                (integral 0.0 tau integrand))))))


    (define (H_nu T-eff nu tau)
      (let* ((C (/ h k T-eff))
             (alpha (* C nu))
             (H (/ (* sigma (expt T-eff 4)) (* 4 pi)))
             (H_alpha (* (H_alpha_hat alpha tau) H))
             (H_nu (* C H_alpha)))
        H_nu))

    (define (B_alpha_hat alpha theta)
      (if (= alpha 0.0)
        0.0
        (* (/ 15 (expt pi 3)) (/ (expt alpha 3) (expt theta 4) (- (exp (/ alpha theta)) 1)))))
             
  ))