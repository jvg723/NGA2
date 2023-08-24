module myrandom
   use random, only: random_normal
   use precision, only: WP
   implicit none
   private

   REAL(WP), PRIVATE :: zero = 0.0_WP, half = 0.5_WP, one = 1.0_WP,  &
   vsmall = TINY(1.0_WP)

   ! Function and subroutines
   public :: random_gamma

contains

   FUNCTION random_gamma(s, first) RESULT(fn_val)

      ! Adapted from Fortran 77 code from the book:
      !     Dagpunar, J. 'Principles of random variate generation'
      !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

      !     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
      !     CALLS EITHER random_gamma1 (S > 1.0_WP)
      !     OR random_exponential (S = 1.0_WP)
      !     OR random_gamma2 (S < 1.0_WP).

      !     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL(WP)).

      REAL(WP), INTENT(IN)    :: s
      LOGICAL, INTENT(IN) :: first
      REAL(WP)                :: fn_val

      IF (s <= zero) THEN
         WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
         STOP
      END IF

      IF (s > one) THEN
         fn_val = random_gamma1(s, first)
      ELSE IF (s < one) THEN
         fn_val = random_gamma2(s, first)
      ELSE
         fn_val = random_exponential()
      END IF

      RETURN
   END FUNCTION random_gamma


   FUNCTION random_gamma1(s, first) RESULT(fn_val)

      ! Uses the algorithm in
      ! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
      ! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

      ! Generates a random gamma deviate for shape parameter s >= 1.

      REAL(WP), INTENT(IN)    :: s
      LOGICAL, INTENT(IN) :: first
      REAL(WP)                :: fn_val

      ! Local variables
      REAL(WP), SAVE  :: c, d
      REAL(WP)        :: u, v, x

      IF (first) THEN
         d = s - one/3.
         c = one/SQRT(9.0_WP*d)
      END IF

      ! Start of main loop
      DO

         ! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

         DO
            x = random_normal(0.0_WP,1.0_WP)
            v = (one + c*x)**3
            IF (v > zero) EXIT
         END DO

         ! Generate uniform variable U

         CALL RANDOM_NUMBER(u)
         IF (u < one - 0.0331_WP*x**4) THEN
            fn_val = d*v
            EXIT
         ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
            fn_val = d*v
            EXIT
         END IF
      END DO

      RETURN
   END FUNCTION random_gamma1



   FUNCTION random_gamma2(s, first) RESULT(fn_val)

      ! Adapted from Fortran 77 code from the book:
      !     Dagpunar, J. 'Principles of random variate generation'
      !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

      ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
      ! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
      ! GAMMA2**(S-1) * EXP(-GAMMA2),
      ! USING A SWITCHING METHOD.

      !    S = SHAPE PARAMETER OF DISTRIBUTION
      !          (REAL(WP) < 1.0_WP)

      REAL(WP), INTENT(IN)    :: s
      LOGICAL, INTENT(IN) :: first
      REAL(WP)                :: fn_val

      !     Local variables
      REAL(WP)       :: r, x, w
      REAL(WP), SAVE :: a, p, c, uf, vr, d

      IF (s <= zero .OR. s >= one) THEN
         WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
         STOP
      END IF

      IF (first) THEN                        ! Initialization, if necessary
         a = one - s
         p = a/(a + s*EXP(-a))
         IF (s < vsmall) THEN
            WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
            STOP
         END IF
         c = one/s
         uf = p*(vsmall/a)**s
         vr = one - vsmall
         d = a*LOG(a)
      END IF

      DO
         CALL RANDOM_NUMBER(r)
         IF (r >= vr) THEN
            CYCLE
         ELSE IF (r > p) THEN
            x = a - LOG((one - r)/(one - p))
            w = a*LOG(x)-d
         ELSE IF (r > uf) THEN
            x = a*(r/p)**c
            w = x
         ELSE
            fn_val = zero
            RETURN
         END IF

         CALL RANDOM_NUMBER(r)
         IF (one-r <= w .AND. r > zero) THEN
            IF (r*(w + one) >= one) CYCLE
            IF (-LOG(r) <= w) CYCLE
         END IF
         EXIT
      END DO

      fn_val = x
      RETURN

   END FUNCTION random_gamma2


   FUNCTION random_exponential() RESULT(fn_val)

      ! Adapted from Fortran 77 code from the book:
      !     Dagpunar, J. 'Principles of random variate generation'
      !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

      ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
      ! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
      ! TO EXP(-random_exponential), USING INVERSION.

      REAL(WP)  :: fn_val

      !     Local variable
      REAL(WP)  :: r

      DO
         CALL RANDOM_NUMBER(r)
         IF (r > zero) EXIT
      END DO

      fn_val = -LOG(r)
      RETURN

   END FUNCTION random_exponential

end module myrandom
