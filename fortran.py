    def _get_fortran( self, time, steps ):
        def str2f(str):
            ret='      %s\n' %(str)
            if len(str) > 61:
                ret=''
                pos=0
                cont=' '
                while ( pos < len(str) ):
                    if pos>0:
                        cont = '&'
                    ret += '     %s%s\n' %( cont,str[ pos: pos+66 ])
                    pos += 66
            return ret
        s0 = self.get_initial_conc()
        ode_func=''
        initial_values=''
        # reaction velocities
        for i,kl in enumerate( [ r.getKineticLaw()for r in self._model.getListOfReactions() ] ):
            math=kl.getMath()
            formula = self._ast_to_string( math)
            ode_func += str2f( 'V(%s) = %s' %(i+1,formula) )
            

        # dS/dt
        for s in range(self._model.getNumSpecies()):
            ode = 'ydot(%s) = 0. ' %(s+1)
            for r in range(self._model.getNumReactions()):
                if self._N[s,r] != 0:
                    ode += ' + %s * v(%s)' %(self._N[s,r], r+1)

            ode_func        += str2f( ode )
            initial_values  += str2f( 'y(%s) = %s' %( s+1, s0[s])  )

        
        num_sp     = self._model.getNumSpecies()
        stepsize   = float(time)/steps

        config_dict = { 'num_sp'          :  num_sp,\
                        'num_reac'        : self._model.getNumReactions(),\
                        'lrw'             :  22 + ( num_sp * max(16, 9+num_sp ) ),\
                        'liw'             :  20 + num_sp,\
                        'initial_values'  :  initial_values,\
                        'ode_func'        :  ode_func,\
                        'stepsize'        :  stepsize,\
                        'steps'           :  steps }

        f_code="""
      SUBROUTINE GET_RESULT (RET)
      EXTERNAL FEX
      DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y, RET(%(steps)s,%(num_sp)s)
cf2py intent(in,out) ret
      DIMENSION Y(%(num_sp)s), RWORK(%(lrw)s), IWORK(%(liw)s)
      NEQ = %(num_sp)s
%(initial_values)s
      T = 0.
      TOUT = 0.
      ITOL = 1
      RTOL = 1.D-6
      ATOL = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = %(lrw)s
      LIW = %(liw)s
      JT = 2
      DO 40 IOUT = 1,%(steps)s
        CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
        DO 50 COUNT = 1,%(num_sp)s
 50       RET(IOUT,COUNT) = Y(COUNT)

 40     TOUT = TOUT + %(stepsize)s
      END

 
      SUBROUTINE FEX (NEQ, T, Y, YDOT)
      DOUBLE PRECISION T, Y(%(num_sp)s), YDOT(%(num_sp)s), V(%(num_reac)s)
%(ode_func)s
      RETURN
      END
""" %config_dict
        return f_code


    def compile_model(self, time, steps):
        code = self._get_fortran(time,steps)
        f = open('tmp/out2.f', 'w')
        f.write(code)
        f.close()

        #cmd = 'gfortran -Wl,-rpath,odepack/ -Lodepack/ -lodepack -o tmp/run_model tmp/out2.f'
        cmd = 'f2py -c -m run_sbml -lodepack -L. tmp/out2.f'
        print os.popen(cmd).read()

    def integrate_f(self, steps):

        a = numpy.zeros( (steps,self._model.getNumSpecies()) )
        return run_sbml.get_result(a)
