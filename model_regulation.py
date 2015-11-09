#!/usr/bin/env python

import numpy, scipy, pp, libsbml, sys, pylab, matplotlib
import sbml_mca, regulator, misc, model_variation


_hosts                  = ()#('fish', 'barbabea', 'marvin', 'trillian', 'ford')
_password               = 'katze'
#_monte_carlo_iterations = 10


class model_regulation:
    """ Class to examine the influence of regulations """
    def __init__(self, model, evaluator_class, evaluator_params, dev, mod_type, k_mod, params=None, server=None ):
        """
        @type    model:               libsbml.model
        @param   model:               The libsbml model
        @type    evaluator_class:     evaluator_base
        @param   evaluator_class:     reference to the class to be used for evaluation
        @type    evaluator_params:    tuple
        @param   evaluator_params:    parameters to be passed to the evaluator constructor
        @type    dev:                 number or numpy.array
        @param   dev:                 deviation (if number then equal for all parameters) otherwise covariance matrix
        @type    mod_type:            string
        @param   mod_type:            indicate wheter inhibition ('inh') or activations ('act')
        @type    k_mod:               number
        @param   k_mod:               inhibition or activation constant
        @type    params:              list
        @param   params:              list of parameter ids to be varied, if not specified the external metaboltiexs are being varied
        @type    server:              pp.Server
        @param   server:              parallel python server object (optional)
        """
        self._model                 =  model
        self._sbml_mca              =  sbml_mca.sbml_mca(model)
        self._eval_class            =  evaluator_class
        self._eval_params           =  evaluator_params
        self._dev                   =  dev
        self._mod_type              =  mod_type
        self._k_mod                 =  k_mod
        if params:
            for p in params:
                if misc.get_parameter_value(self._model, p )==None:
                    raise Exception('Error: parameter %s not found' %p)
            self._params            =  params
        else:
            self._params            =  self._sbml_mca._external_species_ids
        if server:
            self._server            =  server
        else:
            s=self._server            =  pp.Server( ppservers=_hosts, secret=_password )

        self._pos_b = [r.getId() for r in self._model.getListOfReactions()].index('reaction_8')
        print self._model.getReaction(self._pos_b).getName()
            
        
    def eval_flux_curvature( self ):
        """ Get the curvature of the flux w.r. to parameters"""
        rj2 = self._sbml_mca.get_2nd_custom_flux_resp( self._params )
        print 'Curvature of flux w.r.t to external parameters:'
        for pos,r in enumerate([r.getId() for r in self._model.getListOfReactions()]):
            print r
            for p_pos,p in enumerate(self._params):
                print '%s: %s' %(p, rj2[pos][p_pos][p_pos])
            print
        print

    def eval_variance_change_for_regulations(self,flux=1,conc=1):
        """ Evaluate how the different regulations influence the variance of the steady state flux and metabolite concnetrations """
        import heapq
        
        rp = misc.get_positions_enzyme_catalysed_reactions( self._model )
        nep = self._sbml_mca._get_not_enzyme_positions()
        reaction2pos = dict( zip( [self._model.getReaction(i).getId() for i in rp], range(len(rp)) ) )
        species2pos  = dict( zip([ self._sbml_mca._species_ids[x] for x in nep ], range(len(nep))) )


        def interpred( unreg, variants ):
            reaction_q = []
            regulator_q = []
            comb_q = []
            sp_count_reduced={}
            for r_id in reaction2pos:
                r_count_reduced=0.
                for s_id in species2pos:
                    comb_count_reduced=0.
                    #print r_id, 'regulated by', s_id, ':'
                    diff = abs(variants[species2pos[s_id]][reaction2pos[r_id]]) - abs(unreg)
                    #print diff
                    
                    #red = sum(sum( diff < 0 ))
                    #r_count_reduced += red
                    #comb_count_reduced += red
                    #try:
                    #    sp_count_reduced[s_id]+=red
                    #except:
                    #    sp_count_reduced[s_id]=red
                        
                    for pos,p in enumerate(self._params):
                        #print p, diff.T[pos]
                        red = diff.T[pos][self._pos_b]
                        r_count_reduced += red
                        comb_count_reduced += red
                        try:
                            sp_count_reduced[s_id]+=red
                        except:
                            sp_count_reduced[s_id]=red
                        #if (diff.T[pos]<0).all():
                        #    print '  -> all ... have reduced variability w.r.t %s for this regulation' %p
                        #    pass
                        #if (diff.T[pos]>0).all():
                        #    print '  -> all ... have enhanced variability w.r.t %s for this regulation' %p
                        #    pass
                        #print
                        
                    #heapq.heappush( comb_q, (comb_count_reduced, r_id, s_id) )
                    comb_q.append( (comb_count_reduced, r_id, s_id) )
                #heapq.heappush( reaction_q, (r_count_reduced, r_id) )
                reaction_q.append( (r_count_reduced, r_id) )
            for sp in sp_count_reduced: 
                #heapq.heappush( regulator_q, (sp_count_reduced[sp],sp) )
                regulator_q.append( (sp_count_reduced[sp],sp) )

            comb_q.sort()
            reaction_q.sort()
            regulator_q.sort()

            print '*******************************'
            print 'ranking'
            print 'combinations:'
            for i, (rank, r, s) in enumerate(comb_q):
                print i, rank, r, s
            print
            print 'reactions:'
            for i, (rank, r) in enumerate(reaction_q):
                print i, rank, r
            print
            print 'regulators:'
            for i, (rank, s) in enumerate(regulator_q):
                print i, rank, s
            print
            
        
        if flux:
            print 'Influence of regulation on flux variability:'
            print [ r.getId() for r in self._model.getListOfReactions() ]
            unreg,variants = self._apply_func_to_model_variants('get_custom_flux_resp', None, (self._params,))
            interpred( unreg, variants)

        
        if conc:
            print 'Influence of regulation on concentration variability:'
            print self._sbml_mca._species_ids
            unreg,variants = self._apply_func_to_model_variants('get_custom_conc_resp', None, (self._params,))
            interpred( unreg, variants )

        print
        print

        

    def eval_regulation_patterns(self, output='relative'):
        """ generate model variants and evaluate the fitness model for them """
        nep = self._sbml_mca._get_not_enzyme_positions()
        #m = self._eval_model_variants(mode=output)
        eval_class_wrap = lambda model: self._eval_class( model, self._eval_params, self._params, self._dev )
        ref, mat = self._apply_func_to_model_variants( 'evaluate_model', eval_class_wrap, () )
        mat = numpy.array(mat)
        print 
        print 'unregulated: ', ref
        if output=='relative':
            mat = ref-mat
        print misc.matrix2string( mat,\
                                  [r.getId() for r in misc.get_enzyme_catalysed_reactions(self._model)],\
                                  [ self._sbml_mca._species_ids[x] for x in nep ] )
        return mat

    def _apply_func_to_model_variants(self, func, obj_class=None, params=() ):
        # evaluate func for every regulation pattern
        nep           =  self._sbml_mca._get_not_enzyme_positions()
        variants      =  self._generate_model_variants()
        rp            =  misc.get_positions_enzyme_catalysed_reactions( self._model )
        reaction2pos  =  dict( zip( [self._model.getReaction(i).getId() for i in rp], range(len(rp)) ) )
        species2pos   =  dict( zip([ self._sbml_mca._species_ids[x] for x in nep ], range(len(nep))) )
        nr            =  len(misc.get_enzyme_catalysed_reactions(self._model))
        result        =  [[None]*nr for i in range(len(nep))]
        jobs=[]
        for (reaction,species,model_sim) in variants:
            if obj_class:
                o = obj_class(model_sim)
            else:
                o = model_sim
            sys.stderr.write('subm')
            j = self._server.submit( getattr(o,func), params, (), ('numpy', 'scipy' ), globals=globals() )
            jobs.append((reaction,species,j))
            
        #for i,(reaction,species,model_sim) in enumerate(variants):
        for reaction, species, j in jobs:
            print reaction, species
            if species=='unregulated':
                #unreg=jobs[i]()
                unreg=j()
                continue
            sp_pos = species2pos[species]
            r_pos  = reaction2pos[reaction]
            result[sp_pos][r_pos]=j() #obs[i]()
        return (unreg,result)
    
    def _generate_model_variants(self, regulators=None):
        ss_ref         =  self._sbml_mca.get_steady_state()
        if regulators==None:
            regulators = self._sbml_mca._species_ids
        yield ('unregulated', 'unregulated', self._sbml_mca)
        for reac in self._model.getListOfReactions():
            if reac.getNumModifiers()==0:   # non enzymatic reactions can not be regulated
                continue
            for species in regulators:
                if 'enzyme' in species:  continue
                m = regulator.add_regulation( self._model, reac.getId(), species, self._mod_type, self._k_mod,  ss_ref)
                yield ( reac.getId(), species, sbml_mca.sbml_mca(m) ) 



class evaluator_base(object):
    """ Base class for evaluators """
    def __init__(self, sbml_mca, fitness_params, parameter_ids, cov ):
        raise Exception('Implement me')
    
    def evaluate_model(self):
        raise Exception('Implement me')


class evaluator_std_taylor:
    """ Implementation of the standard fitness models with taylor approx. """
    def __init__(self, sbml_mca, fitness_params, parameter_ids, cov ):
        self._sbml_mca                                             =  sbml_mca
        ( self._y_id, self._h, self._alpha, self._ep, self._nep )  =  fitness_params
        self._param_ids                                            =  parameter_ids
        if type(cov)==type(numpy.array(1)):
            self._cov = cov
        else:
            self._cov = numpy.eye(len(parameter_ids)) * float(cov)
            
    def evaluate_model(self):
        """ return linear approximated fitness value """
        ss            =  self._sbml_mca.get_steady_state()
        flux          =  self._sbml_mca._v(ss,1)
        enzyme_conc   =  self._sbml_mca.get_initial_conc()[self._ep]
        #p             =  self._sbml_mca.get_external_metabolite_concentrations()
        p             =  self._sbml_mca.get_parameter_values( self._param_ids )
        F             =  self._f( ss, flux, enzyme_conc )
        Fpp           =  self._dpp_f( ss, flux, enzyme_conc, p )
        return F + numpy.trace( numpy.dot( Fpp, self._cov ) )
        
    def _f( self, ss, flux, enzyme_conc ):
        """ standard fitness function """
        y = numpy.hstack((ss[self._nep], flux))
        v = numpy.dot( numpy.sqrt(self._h), y - self._y_id )
        fp = -.5 * numpy.dot( v, v )
        fm = numpy.dot( self._alpha, enzyme_conc )        
        return fp - fm

    
    def _dpp_f( self, ss, flux, enzyme_conc, p ):
        """ get second defivative w.r.t external metabolite concentrations """
        #s_ext_names    =  self._sbml_mca._external_species_ids
        #p_old          =  self._sbml_mca.get_external_metabolite_concentrations()
        p_old          =   self._sbml_mca.get_parameter_values(self._param_ids)
        
        RS             =   self._sbml_mca.get_custom_conc_resp( self._param_ids )[self._nep]
        RJ             =   self._sbml_mca.get_custom_flux_resp( self._param_ids )
        RY             =   numpy.concatenate((RS,RJ))
        size_p         =   len(self._param_ids)
        R2S            =   self._sbml_mca.get_2nd_conc_resp()[self._nep, :size_p, :size_p]
        R2J            =   self._sbml_mca.get_2nd_flux_resp()[:, :size_p, :size_p]
        R2Y            =   numpy.concatenate( (R2S,R2J) )
        # evalutate fitness func + derivatives
        F              =   self._f(ss,flux,enzyme_conc)
        Fy             =   self._dy_f(ss,flux,enzyme_conc)
        Fyy            =   self._dyy_f( ss, flux, enzyme_conc )
        # calculate 2nd parameter derivative from 1st and 2nd Y derivative
        Fpp   = numpy.dot( numpy.dot( RY.T, Fyy ), RY ) + numpy.tensordot( Fy, R2Y, [0,0] )
        #self._sbml_mca.set_external_metabolite_concentrations( p_old )
        self._sbml_mca.set_parameter_values( self._param_ids, p_old )
        return Fpp

    def _dy_f( self, ss, flux, enzyme_conc ):
        y = numpy.hstack((ss[self._nep], flux))
        """ derivative w.r.t Y of standard fitness function"""
        return -numpy.dot( self._h, y - self._y_id )
    
    def _dyy_f( self, ss, flux, enzyme_conc ):
        return -self._h

class evaluator_std_monte_carlo:
    """ Implementation of the standard fitness model with monte carlo """
    def __init__(self, sbml_mca, fitness_params, parameter_ids, cov ):
        self._sbml_mca                                             =  sbml_mca
        ( self._y_id, self._h, self._alpha, self._ep, self._nep )  =  fitness_params
        self._param_ids                                            =  parameter_ids
        if type(cov)==type(numpy.array(1)):
            self._cov = cov
        else:
            self._cov = numpy.eye(len(parameter_ids)) * float(cov)

    def evaluate_model(self):
        _monte_carlo_iterations=100
        enzyme_conc   =  self._sbml_mca.get_initial_conc()[self._ep]
        mu            =  self._sbml_mca.get_parameter_values( self._param_ids )
        rand          =  model_variation.multivariate_distr(mu, self._cov, mode='normal')
        result        =  numpy.zeros(_monte_carlo_iterations)
        for i in range(_monte_carlo_iterations):
            p = rand.draw()
            self._sbml_mca.set_parameter_values(self._param_ids, p)
            ss            =  self._sbml_mca.get_steady_state()
            flux          =  self._sbml_mca._v(ss,1)
            result[i]     =  self._f( ss, flux, enzyme_conc )
        return sum(result)/_monte_carlo_iterations


    def _f( self, ss, flux, enzyme_conc ):
        """ standard fitness function """
        y = numpy.hstack((ss[self._nep], flux))
        v = numpy.dot( numpy.sqrt(self._h), y - self._y_id )
        fp = -.5 * numpy.dot( v, v )
        fm = numpy.dot( self._alpha, enzyme_conc )        
        return fp - fm



class params_eval_std:
    def get_parameters( self, model, factor_y_opt, h, force=False):
        model_sim            =  sbml_mca.sbml_mca(model)
        ep                   =  model_sim._get_enzyme_positions()
        nep                  =  model_sim._get_not_enzyme_positions()
        self._force          = force
        if h==None:
            sys.stderr.write('No weight vector given, using unitary\n')
            h  =  numpy.eye( len(model_sim._species_ids) - len(ep) + len(model_sim._kinetic_laws) )
        else:
            h  =  numpy.diag( h )
        ss                   =  model_sim.get_steady_state()
        y_id                 =  numpy.hstack((ss[nep], model_sim._v(ss,1) )) * factor_y_opt
        alpha                =  self._get_alpha( model_sim, ep, nep, h, y_id )
        return( y_id, h, alpha, ep, nep )
    
    def _get_alpha( self, model_sim, ep, nep, h, y_id ):
        """ compute weights for the enzymes so that fitness function is in optimal state at reference state"""
        enzyme_names = [ model_sim._species_ids[x] for x in ep ]
        RS          =  model_sim.get_custom_conc_resp( enzyme_names )[nep]
        RJ          =  model_sim.get_custom_flux_resp( enzyme_names )
        RY          =  numpy.concatenate( (RS, RJ) )
        ss          =  model_sim.get_steady_state()
        ss_reduced  =  ss[nep]
        flux        =  model_sim._v(ss,1)
        y           =  numpy.concatenate( (ss_reduced,flux) )
        try:
            alpha = numpy.dot( numpy.dot( - h, (y - y_id) ), RY )
        except:
            print RY.shape
            sys.exit()            
        if sum(x<0 for x in alpha)>0 and not self._force:
            print alpha
            #DEGUG
            raise Exception('Error: negative weights for enzymes.')        
        return alpha



if __name__=='__main__':
    import optparse
    parser = optparse.OptionParser(usage='usage: %prog [options] model')
    parser.add_option('-m', '--mode', metavar='(taylor | carlo)', dest='mode', default='taylor', help='Evaluation method')
    parser.add_option('-t', '--modifier_type', metavar='(inh | act)', dest='modifier_type', help='specify modifier type')
    parser.add_option('-k', '--k_mod', dest='k_mod', help='specify inhibition / activation constant')
    parser.add_option('-y', '--factor_y_opt', dest='fac_y_opt', help='factor for optimal ss values')
    parser.add_option('--cov', dest='cov', help='File with parameter covariance')
    parser.add_option('--mean', dest='mean', help='File with parameter mean' )
    parser.add_option('--dev', dest='dev', help='Same devation for all parameters')
    parser.add_option('-p', '--param_names', dest='params', help='Specify parameters')
    parser.add_option('-H', '--weight_vector', dest='h', help='Specify weight vector for the steady state values')
    parser.add_option('-f', '--force', action='store_true', dest='force',  help='continue with negative alphas')
    parser.add_option( '--flux_curv', action='store_true', dest='flux_curv', help='Compute flux curvature w.r.t parameters' )
    parser.add_option( '--flux_var', action='store_true', dest='flux_var', help='Compute change of flux variance for the different regulations' )
    parser.add_option( '--conc_var', action='store_true', dest='conc_var', help='Compute change of concentration variance for the different regulations' )
    parser.add_option( '--fitness', action='store_true', dest='fitness', help='Evaluate fitness modell for different regulations' )
    
    
                      
    (options,args) = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    h = cov = mean = params = None
    if options.dev!=None:
        cov = options.dev
    if options.cov:
        cov = misc.load_matrix_from_file(options.cov)
    if options.mean:
        mean = misc.load_matrix_from_file(options.mean)[0]
    if options.params:
        params=options.params.split()

    if options.mode=='taylor':
        eval_class = evaluator_std_taylor
    elif options.mode=='carlo':
        eval_class = evaluator_std_monte_carlo
    else:
        raise Exception('Unknown mode')

    if options.k_mod:
        k_mod=float(options.k_mod)
    else:
        print 'Error: specify k_mod'
        sys.exit()
    if not options.modifier_type in ['inh','act']:
        print 'Error: specify modifier type'
        sys.exit()
    if options.fac_y_opt!=None:
        fac_y_opt=options.fac_y_opt
    else:
        print 'Setting factor y opt to 1.2'
        fac_y_opt=1.2

    if cov==None:
        print 'Error: specify covariance matrix or deviation'
        sys.exit()
        
        
    if options.h:
        h = [ float(x) for x in options.h.split() ]
        
    doc           =  libsbml.readSBML(sys.argv[-1])
    model         =  doc.getModel()
    pe            =  params_eval_std()
    eval_params   = pe.get_parameters( model, fac_y_opt, h, options.force )
    mr            = model_regulation( model, eval_class, eval_params, cov, options.modifier_type, k_mod, params )

    if options.flux_curv:
        mr.eval_flux_curvature()
    if options.flux_var:
        mr.eval_variance_change_for_regulations(flux=1, conc=0)
    if options.conc_var:
        mr.eval_variance_change_for_regulations(flux=0, conc=1)
    if options.fitness:
        mr.eval_regulation_patterns()



