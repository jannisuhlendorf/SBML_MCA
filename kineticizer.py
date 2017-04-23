#!/usr/bin/env python
import misc
import Thermodynamics
import sys
import libsbml
import numpy


def myeval(formula, d1, d2=None):
    """ replace libsbml power ^ with python power symbol **   """
    return eval(formula.replace('^', '**'), d1, d2)


class Kineticizer(object):
    """ Base class for the assignment of kinetics to SBML models"""
    
    def __init__(self,\
                 model,\
                 ss_conc=None,\
                 ss_flux=None,\
                 vec_vmax=None,\
                 vec_keq=None,\
                 vec_km=None,\
                 vec_kforw=None,\
                 vec_kback=None,\
                 reactions_with_stoich_one=[],\
                 lazy_parameters=False):
        """
        @type     model:                      libsbml.model
        @param    model:                      SBML-model
        @type     ss_conc:                    numpy.array
        @param    ss_conc:                    Vector of metabolite concentrations (without enzymes)
        @type     ss_flux:                    numpy.array
        @param    ss_flux:                    Vector of steady state fluxes
        @type     vec_vmax:                   numpy.array
        @param    vec_vmax:                   Vector of vmax values
        @type     vec_keq:                    numpy.array
        @param    vec_keq:                    Vector of equilibrium constants
        @type     vec_km:                     numpy.arary
        @param    vec_km:                     Vector of KM values
        @type     vec_kforw:                  numpy.array
        @param    vec_kforw:                  Vector of forward rate constants
        @type     vec_kback:                  numpy.array
        @param    vec_kback:                  Vector of backward rate constants
        @type     reactions_with_stoich_one:  list
        @param    reactions_with_stoich_one:  List of indices of reactions where the stoichiometry of all participants should be 1 in the kl
        """
                 
        self._model = model
        misc.add_enzymes_to_reactions(self._model)

        # check input parameters
        if ss_conc != None and ss_flux != None and vec_vmax != None and vec_keq != None:
            params = self._find_parameters(ss_conc, ss_flux, vec_vmax, vec_keq, vec_km)
        elif vec_kforw != None and vec_kback != None: # all parameters specified
            params = self._pack_parameters(vec_kforw, vec_kback, vec_km)
        elif lazy_parameters:
            vec_kforw = numpy.ones(model.getNumReactions())
            vec_kback = numpy.ones(model.getNumReactions())
            vec_km = numpy.ones(sum([r.getNumReactants()+r.getNumProducts() for r in model.getListOfReactions()]))
            params = self._pack_parameters(vec_kforw, vec_kback, vec_km)
        else:
            raise Exception('Not enought parameters specified')

        # define a set of reactions where all stoichiometries are set to 1 in the kinetic law
        self._reactions_with_stoich_one = [r for i,r in enumerate(model.getListOfReactions()) if i in reactions_with_stoich_one]
        if ss_conc != None:
            # set the model to the steady state
            ss_dict = dict(zip([s.getId() for s in misc.get_species_wo_enzymes(self._model)], ss_conc))
            for s in model.getListOfSpecies():
                if ss_dict.has_key(s.getId()):
                    s.setInitialConcentration(ss_dict[s.getId()])
        
        # assign values to empty concentrations
        for s in model.getListOfSpecies():
            if not s.isSetInitialConcentration():
                s.setInitialConcentration(1)
                sys.stderr.write('Warning: Species %s has not set initial concentration. Setting to 1.\n' %s.getId())

        # assign the kinetics
        for r,param in zip(model.getListOfReactions(),params):
            self._assign_kinetic(r, param)

    def _get_stoichiometry(self, species_ref):
        """ use this func instead of the libsbml one if you want to support the funky reaction_with_stoich_one thing """
        r = species_ref.getParentSBMLObject().getParentSBMLObject()
        if r in self._reactions_with_stoich_one:
            return 1
        return species_ref.getStoichiometry()

    def _set_enzyme_concentrations(self, enzyme_conc):
        """ set enzyme concentrations with vector """
        for i,reaction in enumerate(self._model.getListOfReactions()):
            e=misc.get_enzyme_for_reaction(self._model, reaction)
            e.setInitialConcentration(enzyme_conc[i])

    def _set_external_metabolite_conc(self,ss_conc):
        """ set external metabolite concentrations with vector """
        for i,s in enumerate(misc.get_species_wo_enzymes(self._model)):
            if s.getConstant() or s.getBoundaryCondition():
                s.setInitialConcentration(ss_conc[i])

    def _get_enzyme(self, reaction):
        """ get an enzyme for a reaction or make one"""
        e = misc.get_enzyme_for_reaction(self._model, reaction)
        if not e:
            s = self._model.createSpecies()
            s.setId('enzyme_'+reaction.getId())
            e.setInitialConcentration(1)
            mod_ref = reaction.createModifier()
            e = s.getId()
            mod_ref.setSpecies(e)
        return e

    def _test_thermodynamics(self, ss, flux, keq, enzymes):
        td = Thermodynamics.thermodynamics(self._model)
        nep = misc.get_not_enzyme_positions(self._model)
        if not td.check_flux_signs(ss[nep], flux, keq):
            raise Exception('Flux signs are thermodynamically not feasible')

    def _get_standard_kms(self):
        kms = numpy.ones(sum([r.getNumReactants()+r.getNumProducts() for r in model.getListOfReactions()]))
        return kms

    def _pack_parameters(self, vec_kforw, vec_kback, vec_km):
        """ Abstract method to bring input parameters to a structured form """
        raise Exception('Implement me')
        
    def _find_parameters(self, ss_conc, fluxes, eq, kms=None):
        """ Abstract method to find a suitable parameter set. (returns list of parameters for each reaction)"""
        raise Exception('Implement me')

    def _assign_kinetic(self, reaction, params):
        """ Abstract method to assign a new kinetic to the reaction with parameters  """
        raise Exception('Implement me')


class MassAction(Kineticizer):
    """ Class to assign mass actions kinetics to SBML models """
    
    def _find_parameters(self, ss_conc, fluxes, kvs, eq, kms):
        species = misc.get_species_wo_enzymes(self._model)
        species2pos = dict( zip([s.getId() for s in species], range(len(species)) ) )
        params=[]
        for i,reaction in enumerate(self._model.getListOfReactions()):
            v = fluxes[i]
            e = self._get_enzyme(reaction).getInitialConcentration()
            q = eq[i]
            f_term = b_term =0.
            for s in reaction.getListOfReactants():
                f_term += ss_conc[species2pos[s.getSpecies()]] ** self._get_stoichiometry(s)
            for s in reaction.getListOfProducts():
                b_term += ss_conc[species2pos[s.getSpecies()]] ** self._get_stoichiometry(s)
            k_plus = v/(e*( f_term - (b_term/q) ))
            k_minus = k_plus/q
            params.append((k_plus,k_minus))
        return params

    def _assign_kinetic(self, reaction, params):
        (kf, kb) = params
        r_id = reaction.getId()
        enzyme = self._get_enzyme(reaction).getId()
            
        term = {'Reactants': '( kf_%s ' %r_id, 'Products': '( kb_%s '%r_id}
        for type in ['Reactants', 'Products']:
            for species in getattr(reaction, 'getListOf' + type)():
                term[type] += ' * %s' %(species.getSpecies())
                if self._get_stoichiometry(species) != 1:
                    term[type]+='^%s' %(self._get_stoichiometry(species))
        kin = '%s * ( %s ) - %s ) )' %(enzyme, term['Reactants'], term['Products'])
        kl = libsbml.KineticLaw()
        kl.setFormula(kin)
        k_f = libsbml.Parameter('kf_'+r_id, kf)
        kl.addParameter(k_f)
        k_b = libsbml.Parameter('kb_'+r_id, kb)
        kl.addParameter(k_b)
        reaction.setKineticLaw(kl)
        
    def _pack_parameters(self, vec_kforw, vec_kback, vec_km):
        return zip(vec_kforw,vec_kback)


class CompleteRandomOrder(Kineticizer):
    """ Class to assign complete random order kinetics to SBML models """
    
    def _find_parameters(self, ss_conc, fluxes, kvs, eqs, kms):
        ss_dict = dict(zip([s.getId() for s in misc.get_species_wo_enzymes(self._model)], [float(x) for x in ss_conc]))
        params=[]
        value_dict={}
        count=0
        # find kf,kb with formula (4), page 5 in manuscript
        for i, reaction in enumerate(self._model.getListOfReactions()):
            r_id=reaction.getId()
            H=1.
            km_sub=[]
            km_prod=[]
            for s in reaction.getListOfReactants():
                km_sub.append(kms[count])
                H *= kms[count]** self._get_stoichiometry(s)
                name = 'kM_%s_%s' %(r_id,s.getSpecies())
                value_dict[name]= kms[count]
                count+=1
                
            for s in reaction.getListOfProducts():
                km_prod.append(kms[count])
                H *= kms[count]**-self._get_stoichiometry(s)
                name = 'kM_%s_%s' %(r_id,s.getSpecies())
                value_dict[name]= kms[count]
                count+=1
                
            kf = kvs[i]*(eqs[i]*H)**.5
            kb = kvs[i]*(eqs[i]*H)**-.5

            # find enzyme conc
            forward,backward,denominator  = self._get_formula_factors(reaction, 1)
            print 
            print r_id
            print forward
            print backward
            print denominator
            forward_v = myeval(forward, value_dict, ss_dict)
            backward_v = myeval(backward, value_dict, ss_dict)
            denominator_v = myeval(denominator, value_dict, ss_dict)
            
            v = (kf*forward_v - kb*backward_v)/denominator_v
            enzyme_conc = fluxes[i]/v
            
            if enzyme_conc<0:
                raise Exception('Error: Flux, concnetrations of equilibrium constants do not fit.')
            enzyme = misc.get_enzyme_for_reaction(self._model,reaction)
            enzyme.setInitialConcentration(enzyme_conc)
            params.append((km_sub,km_prod, kf, kb))
        return params

    def _assign_kinetic(self, reaction, params):
        (km_sub, km_prod, kf, kb) = params
        r_id = reaction.getId()
        enzyme = self._get_enzyme(reaction).getId()

        # assemble the formula
        kl = libsbml.KineticLaw()
        (forward,backward,denominator) = self._get_formula_factors( reaction, enzyme )
        formula = ' %s * ( kf_%s * %s - kb_%s * %s) / ( %s )' %(enzyme, r_id, forward, r_id, backward, denominator)
        kl.setFormula(formula)
        
        # set the parameters
        k_f = libsbml.Parameter('kf_'+r_id, kf)
        kl.addParameter(k_f)
        k_b = libsbml.Parameter('kb_'+r_id, kb)
        kl.addParameter(k_b)
        
        for i,s in enumerate(reaction.getListOfReactants()):
            pid = 'kM_%s_%s' %(r_id, s.getSpecies())
            km = libsbml.Parameter(pid, km_sub[i])
            kl.addParameter(km)
        for i,s in enumerate(reaction.getListOfProducts()):
            pid = 'kM_%s_%s' %(r_id, s.getSpecies())
            km = libsbml.Parameter(pid, km_prod[i])
            kl.addParameter(km)        
        reaction.setKineticLaw(kl)

    def _get_formula_factors(self, reaction, enzyme):
        r_id = reaction.getId()
        denom=[]
        term = {'Reactants': [], 'Products': []}
        for type in ['Reactants', 'Products']:
            for species in getattr(reaction, 'getListOf' + type)():
                s_id = species.getSpecies()
                km  = 'kM_%s_%s'  %(r_id,s_id)
                x   = '(%s/%s)'   %(s_id, km)
                d   = '(1 + %s)'  %x
                if self._get_stoichiometry(species) != 1:
                    stoich = '^%s' %(self._get_stoichiometry(species))
                    x +=  stoich
                    d +=  stoich
                term[type].append(x)
                denom.append(d)

        forward = ' * '.join(term['Reactants'])
        backward = ' * '.join(term['Products'])
        denominator = ' * '.join(denom)
        return (forward, backward, denominator)

    def _pack_parameters(self, vec_kforw, vec_kback, vec_km):
        if type(vec_km) != list:
            vec_km = vec_km.tolist()
        vec_km.reverse()
        km_substr = []
        km_prod = []
        # sorry for this
        [(km_substr.append([vec_km.pop() for i in range(x)]),\
          km_prod.append([vec_km.pop() for i in range(y)]))\
                    for x,y in [ (r.getNumReactants(),r.getNumProducts())\
                                 for r in self._model.getListOfReactions()]]
        return zip(km_substr, km_prod, vec_kforw, vec_kback)


class ConvenienceKinetics(CompleteRandomOrder):

    """ Class to assign convenience kinetics to SBML models """
    def _get_formula_factors(self, reaction, enzyme):
        r_id = reaction.getId()
        denom = []
        term = {'Reactants': [], 'Products': []}
        for type in ['Reactants', 'Products']:
            for species in getattr(reaction, 'getListOf' + type)():
                s_id = species.getSpecies()
                km = 'kM_%s_%s'  %(r_id,s_id)
                x = '(%s/%s)'   %(s_id, km)
                denom_high_stoich = ''
                add_stoich = ''
                if self._get_stoichiometry(species) != 1:
                    denom_high_stoich = ' + ' + ' + '.join(
                        [x+'^'+str(i) for i in range(2,int(self._get_stoichiometry(species)+1))])
                    add_stoich = '^%s' %(self._get_stoichiometry(species))
                d = '( 1 + %s %s )' %(x,denom_high_stoich)
                term[type].append(x + add_stoich)
                denom.append(d)
        forward = ' * '.join( term['Reactants'] )
        backward = ' * '.join( term['Products'] )
        denominator = ' * '.join( denom )
        return (forward, backward, denominator)


def parse_input_file(filename):
    str2arr = lambda str: numpy.array([float(x) for x in str.split()])
    fluxes = ss_conc = vmax = eq_const = kms = None
    f = open(filename, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        if line.lower().startswith('fluxes'):
            fluxes = str2arr( line[6:] )
        elif line.lower().startswith('ss_conc'):
            ss_conc = str2arr( line[7:] )
        elif line.lower().startswith('vmax'):
            vmax = str2arr( line[4:] )
        elif line.lower().startswith('eq_const'):
            eq_const = str2arr( line[8:] )
        elif line.lower().startswith('kms'):
            kms = str2arr( line[3:] )
    return (fluxes, ss_conc, vmax, eq_const, kms)


if __name__=='__main__':
    import optparse
    parser = optparse.OptionParser(usage='usage: %prog [options] model')
    parser.add_option('-p','--param_file', metavar='file', dest='param_file', help='file with kinetic constants')
    parser.add_option('-t', '--type', metavar='kinetic', dest='type',
                      help='type of kinetic to us: ck: convenience kinetics | cro: complete random order | ma mass action')
    parser.add_option('-f', '--fantasy_params', action='store_true', dest='fantasy_params',
                      help='user fantasy parameters for missing values')

    (options,args) = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    doc = libsbml.readSBML(sys.argv[-1])
    model = doc.getModel()

    if options.param_file:
        flux, ss, vmax, keq, km = parse_input_file(options.param_file)
    else:
        flux = ss = vmax = keq = km = None
    
    type2amount = {'flux': model.getNumReactions(),\
                   'ss': len( misc.get_species_wo_enzymes(model) ),\
                   'vmax': model.getNumReactions(),\
                   'keq': model.getNumReactions(),\
                   'km': sum([r.getNumReactants()+r.getNumProducts() for r in model.getListOfReactions()]) }

    for key in type2amount:
        if locals()[key]==None:
            if options.fantasy_params:
                sys.stderr.write('%s parameters missing, setting to 1\n' %key)
                locals()[key] = numpy.ones( type2amount[key] )
            else:
                if key=='km' and options.type=='ma':
                    continue
                sys.stderr.write('%s parameters missing\n' %key)
                sys.exit()
        elif len(locals()[key]) != type2amount[key]:
            sys.stderr.write('Wrong number of parameters for %s\n' %key)
            sys.exit()

    d= {'ck': ConvenienceKinetics,\
        'cro': CompleteRandomOrder,\
        'ma': MassAction}
    if not options.type in d:
        sys.stderr.write('Error: unknown kinetics type\n')
        sys.exit()
    d[options.type]( model, ss_conc=ss, ss_flux=flux, vec_vmax=vmax, vec_keq=keq, vec_km=km, reactions_with_stoich_one=[] )
    
    print '<?xml version="1.0" encoding="UTF-8"?>'
    print doc.toSBML()
