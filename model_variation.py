#!/usr/bin/env python
import numpy
import sbml_mca, misc

class model_variation:
    """ Class to compute the steady state distribution from parameter distribution """
    
    def __init__(self, model, cov_p, mu_p=None, param_names=None, distr='normal', mode='taylor2' ):
        """
        mode   = ( taylor1 | taylor2 | carlo )
        distr  = ( normal | lognormal )
        params = sbml ids of parameters
        """
        self._model      =  model
        self._cov_p      =  cov_p
        self._model_sim  =  sbml_mca.sbml_mca(model)
        self._param_names=param_names
        if param_names==None:
            self._param_names = self._model_sim._external_species_ids
        if mu_p==None:
            self._mu_p = numpy.array([ misc.get_parameter_value(self._model_sim._model, p_name) for p_name in self._param_names ])
        if cov_p==None:
            self._cov_p = numpy.eye(len(self._param_names))
        
        params = ( self._model_sim, self._cov_p, self._mu_p, self._param_names, distr )
        if mode=='carlo':
            self._variator = apply( carlo_variation, params )
        elif mode=='taylor1':
            self._variator = apply( taylor_variation, params + (1,) )
        elif mode=='taylor2':
            self._variator = apply( taylor_variation, params + (2,) )

    def get_flux_distr(self):
        """ get (mean,cov) of flux distribution """
        return self._variator.get_flux_distr()

    def get_conc_distr(self):
        """ get (mean,cov) of concentration distribution """
        return self._variator.get_conc_distr()


class carlo_variation:
    """ Class to compute steady state discribution by monte carlo method """
    def __init__(self, model_sim, cov_p, mu_p, params, distr, order=None):
        """
        @type   model_sim:  sbml_mca
        @param  model_sim:  Instance of sbml_mca object
        @type   cov_p:      numpy.array
        @param  cov_p:      covariance matrix of parameters
        @type   mu_p:       numpy.array
        @param  mu_p:       mean of parameter values
        @type   params:     list
        @param  params:     list of parameter ids
        @type   distr:      string
        @param  distr:      specify the type of distribution (normal or lognormal)
        """
        self._distr     = multivariate_distr( mu_p, cov_p, distr )
        self._params    = params
        self._model_sim = model_sim
        self._cov_conc = self._mu_conc = self._cov_flux = self._mu_flux =None
        
    def get_flux_distr(self):
        """ Get distribution of fluxes """
        if self._cov_flux==None:
            self._get_distr()
        return (self._mu_flux, self._cov_flux)

    def get_conc_distr(self):
        """ Get the distribution of steady state concentrations"""
        if self._cov_conc==None:
            self._get_distr()
        return (self._mu_conc, self._cov_conc)
            
    def _get_distr( self, num_iter=100 ):
        nep = self._model_sim._get_not_enzyme_positions()
        data_conc = numpy.zeros((len(nep),num_iter))
        data_flux = numpy.zeros((self._model_sim._model.getNumReactions(),num_iter))
        for i in range(num_iter):
            p = self._distr.draw()
            #self._model_sim.set_external_metabolite_concentrations(p)
            self._model_sim.set_parameter_values( self._params, p )
            ss=self._model_sim.get_steady_state()
            data_conc.T[i]=ss[nep]
            data_flux.T[i]=self._model_sim._v(ss,1)    
        mean_conc         =  numpy.dot( data_conc, numpy.ones( num_iter ) ) / num_iter
        data_conc_center  =  (data_conc.T-mean_conc).T
        cov_conc          =  numpy.dot( data_conc_center,data_conc_center.T ) / (num_iter-1)
        mean_flux         =  numpy.dot( data_flux, numpy.ones( num_iter ) ) / num_iter
        data_flux_center  =  (data_flux.T-mean_flux).T
        cov_flux          =  numpy.dot( data_flux_center,data_flux_center.T ) / (num_iter-1)
        self._cov_conc  = cov_conc
        self._mu_conc   = mean_conc
        self._cov_flux  = cov_flux
        self._mu_flux   = mean_flux

class taylor_variation:
    """ Class to approximate steady state discribution by taylor expansion"""
    def __init__(self, model_sim, cov_p, mu_p, params, distr, order=None):
        """
        @type   model_sim:  sbml_mca
        @param  model_sim:  Instance of sbml_mca object
        @type   cov_p:      numpy.array
        @param  cov_p:      covariance matrix of parameters
        @type   mu_p:       numpy.array
        @param  mu_p:       mean of parameter values
        @type   params:     list
        @param  params:     list of parameter ids
        @type   distr:      string
        @param  distr:      specify the type of distribution (normal or lognormal)
        """
        self._distr     = multivariate_distr( mu_p, cov_p, distr )
        self._params    = params
        self._model_sim = model_sim
        if distr=='lognormal':
            self._normalize='right'
        else:
            self._normalize=None
        self._order = order

    def get_flux_distr(self):
        """ Get distribution of fluxes """
        RJ   =  self._model_sim.get_custom_flux_resp( self._params, normalize=self._normalize )
        R2J  =  self._model_sim.get_2nd_custom_flux_resp( self._params, normalize=self._normalize )
        ss   =  self._model_sim.get_steady_state()
        mu   =  self._model_sim._v(ss,1)
        return self._taylor_expansion( RJ, R2J, mu, self._order )
    
    def get_conc_distr(self):
        """ Get the distribution of steady state concentrations"""
        nep  =  self._model_sim._get_not_enzyme_positions()
        RS   =  self._model_sim.get_custom_conc_resp( self._params, normalize=self._normalize )[nep]
        R2S  =  self._model_sim.get_2nd_custom_conc_resp( self._params, normalize=self._normalize )[nep]        
        mu   =  self._model_sim.get_steady_state()[nep]
        return self._taylor_expansion( RS, R2S, mu, self._order )

        
    def _taylor_expansion( self, R, R2, mean_y, order):
        C     =  self._distr._cov2
        cov1  =  numpy.dot( numpy.dot(R, C), R.T )
        u     =  numpy.tensordot( R2, C, [1,0] )
        v     =  numpy.tensordot( R2, C, [2,1] )
        cov2  =  .5 * numpy.tensordot( u, v, [(1,2),(2,1)] )
        mu1   =  mean_y
        mu2   =  .5 * numpy.tensordot( R2, C, [(1,2),(0,1)] )
        
        if order==1:
            return (mu1,cov1)
        else:
            return (mu1+mu2,cov1+cov2)

        
class multivariate_distr:
    """ Class providing a multivariate (log)-normal distribution """
    def __init__(self, mu, cov, mode='normal'):
        """
        @type  mu:   numpy.array
        @param mu:   mean vector
        @type  cov:  numpy.array
        @param cov:  covariance matrix
        @type  mode: string
        @param mode: specify the type of distribution (normal or lognormal)
        """
        self._mode=mode
        if mode=='lognormal':
            if numpy.linalg.norm( cov - numpy.diag( cov.diagonal() ))>1e-16:
                raise Exception('Only diagonal covariance matrices supported. If you know the formula tell me ...')
            m = numpy.dot
            d = numpy.diag
            self._cov2 = numpy.log( numpy.ones((mu.size,mu.size))  + m( m( d(1./mu), cov ), d(1./mu) ) )
            self._mu2 = numpy.log(mu) - .5 * self._cov2.diagonal()
        elif mode=='normal':
            self._cov2 = cov
            self._mu2 = mu
        else:
            raise Exception('Unknown mode')
        
    def draw(self):
        """
        Draw random values from the distribution
        @rtype:  numpy.array
        @return: vector
        """
        if self._mode=='lognormal':
            return numpy.exp( numpy.random.multivariate_normal( self._mu2,self._cov2 ) )
        else:
            return numpy.random.multivariate_normal( self._mu2,self._cov2 )




if __name__=='__main__':
    import optparse, sys
    parser = optparse.OptionParser(usage='usage: %prog [options] model')
    parser.add_option('-d', '--distribution', metavar='(normal | lognormal)', dest='distribution', default='normal', help='Specify distribution type of parameters')
    parser.add_option('-m', '--mode', metavar='(taylor1 | taylor2 | carlo)', dest='mode', default='taylor2', help='Approximation method')
    parser.add_option('--cov', dest='cov', help='File with parameter covariance')
    parser.add_option('--mean', dest='mean', help='File with parameter mean' )
    parser.add_option('-c', action='store_true', dest='conc', help='Get distribution of steady state inner metabolites')
    parser.add_option('-j', action='store_true', dest='flux', help='Get distribution of steady state fluxes')
    parser.add_option('-s', '--string_output', action='store_true', dest='string_output', help='Write nice output')
    parser.add_option('-p', '--param_names', dest='params', help='Specify parameters')
                      
    (options,args) = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
        
    cov = mean = params = None
    if options.cov:
        cov = misc.load_matrix_from_file(options.cov)
    if options.mean:
        mean = misc.load_matrix_from_file(options.mean)[0]

    if options.params:
        params=options.params.split()
        
    mv = model_variation( sys.argv[-1],\
                          cov_p=cov,\
                          mu_p=mean,\
                          param_names=params,\
                          distr=options.distribution,\
                          mode=options.mode )


    print 'Parameters:'
    print mv._param_names
    print 'Covariance Parameters:'
    print misc.matrix2string( mv._cov_p, mv._param_names, mv._param_names )
        
    if options.conc:
        mean, cov = mv.get_conc_distr()
        print 'Concentration distribution:'
        species = [ mv._model_sim._species_ids[i] for i in mv._model_sim._get_not_enzyme_positions() ]
        print 'Species: ', species
        print 'Mean:'
        print mean
        print 'Covariance:'
        if options.string_output:
            print misc.matrix2string(cov, species, species )
            #print misc.matrix2string( misc.cov2cor(cov), delim='\t' )
        else:
            print '\t'.join( species )
            print misc.matrix2string( cov, delim='\t' )
        print
        
    if options.flux:
        mean, cov = mv.get_flux_distr()
        print 'Flux distribution:'
        reacs = [ r.getId() for r in mv._model_sim._model.getListOfReactions() ]
        print 'Reactions: ', reacs
        print 'Mean:'
        print mean
        print 'Covariance: ' 
        if options.string_output:
            print misc.matrix2string(cov, reacs, reacs, justify='left')
            #print misc.matrix2string(misc.cov2cor( cov ), justify='left', delim='\t')
            #print misc.matrix2string( cov, justify='left', delim='\t')
        else:
            print cov
        print
        
