c F      = name of subroutine for right-hand side vector f.
c          This name must be declared External in calling program.
c NEQ    = number of first order ODEs.
c Y      = array of initial values, of length NEQ.
c T      = the initial value of the independent variable.
c TOUT   = first point where output is desired (.ne. T).
c ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
c RTOL   = relative tolerance parameter (scalar).
c ATOL   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be less than
c             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
c             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
c          Thus the local error test passes if, in each component,
c          either the absolute error is less than ATOL (or ATOL(i)),
c          or the relative error is less than RTOL.
c          Use RTOL = 0.0 for pure absolute error control, and
c          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
c          control.  Caution: actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c ITASK  = 1 for normal computation of output values of y at t = TOUT.
c ISTATE = integer flag (input and output).  Set ISTATE = 1.
c IOPT   = 0 to indicate no optional inputs used.
c RWORK  = real work array of length at least:
c             22 + NEQ * MAX(16, NEQ + 9).
c          See also Paragraph E below.
c LRW    = declared length of RWORK (in user's dimension).
c IWORK  = integer work array of length at least  20 + NEQ.
c LIW    = declared length of IWORK (in user's dimension).
c JAC    = name of subroutine for Jacobian matrix.
c          Use a dummy name.  See also Paragraph E below.
c JT     = Jacobian type indicator.  Set JT = 2.
c          See also Paragraph E below.
 

      external fex
      integer neq, itol, itask, istate, iopt, lrw, iwork(22), liw, jt,
     1    steps
      double precision y(2), t, tout, rtol, atol(2), rwork(54), stpsize


      neq    =  2
      y(1)   =  1.
      y(2)   =  0.
      t      =  0
      tout   =  0.
      itol   =  1
      rtol   =  1.D-4
      atol   =  1.D-6
      itask  =  1
      istate =  1
      iopt   =  0
      lrw    =  54
      ltw    =  22
      jt     =  2

      steps    =  100
      stpsize  =  0.1



      do 10 iout = 1,steps
         CALL dlsoda(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1        iopt,rwork,lrw,iwork,liw,jdum,jt)

         write(*,*) t, y(1), y(2)
 10      tout = tout+stpsize
         
      stop


 



      stop
      end


      subroutine fex (neq, t, y, ydot)
      integer neq
      double precision t, y(2), ydot(2)
      ydot(1) =   y(2)
      ydot(2) = - y(1)
      return
      end
    
