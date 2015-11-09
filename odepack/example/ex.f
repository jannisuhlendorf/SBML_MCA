

      EXTERNAL FEX
      DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
      DIMENSION Y(9), ATOL(9), RWORK(184), IWORK(29)
      NEQ = 9
      y(1) = 0.
      y(2) = 0.
      y(3) = 0.
      y(4) = 0.
      y(5) = 0.
      y(6) = 0.
      y(7) = 0.
      y(8) = 0.
      y(9) = 0.

      T = 0.
      TOUT = .4
      ITOL = 1
      RTOL = 1.D-4
      ATOL = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = 184
      LIW = 29
      JT = 2
      DO 40 IOUT = 1,12
        CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
        WRITE(*,*)T,Y
        IF (ISTATE .LT. 0) GO TO 80
 40     TOUT = TOUT + 0.1
      STOP
 80   WRITE(6,90)ISTATE
 90   FORMAT(///' Error halt.. ISTATE =',I3)
      STOP
      END
 
      SUBROUTINE FEX (NEQ, T, Y, YDOT)
      DOUBLE PRECISION T, Y(9), YDOT(9), V(2)
      V(1) = Y(8) * (4.84460624534 * (Y(1) / 0.16999) * (Y(2) / 0.16999)
     & - Y(3) / 0.16999 * 5.06381896849) / ((1 + Y(1) / 0.16999) * (1 + 
     &Y(2) / 0.16999) + 1 + Y(3) / 0.16999 - 1)
      V(2) = Y(9) * (0.73971823737 * (Y(4) / 0.21189) * (Y(5) / 0.21189)
     & - Y(6) / 0.13638 * (Y(7) / 0.13638) * 5.8960484407) / ((1 + Y(4) 
     &/ 0.21189) * (1 + Y(5) / 0.21189) + (1 + Y(6) / 0.13638) * (1 + Y(
     &7) / 0.13638) - 1)
      ydot(1) = 0.  + -1.0 * v(1)
      ydot(2) = 0.  + -1.0 * v(1)
      ydot(3) = 0.  + 1.0 * v(1)
      ydot(4) = 0.  + -1.0 * v(2)
      ydot(5) = 0.  + -1.0 * v(2)
      ydot(6) = 0.  + 1.0 * v(2)
      ydot(7) = 0.  + 1.0 * v(2)
      ydot(8) = 0. 
      ydot(9) = 0. 

      RETURN
      END

