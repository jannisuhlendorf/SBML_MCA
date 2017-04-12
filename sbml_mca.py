#!/usr/bin/env python
import sys
import scipy
import scipy.integrate
import numpy
import libsbml
import pylab
import warnings
import time
import misc
import math
import scipy.optimize
import sympy

sbo_type2id = {'activator':           21,
               'inhibitor':           20,
               'enzyme':              14}


class Noncritical_error(Exception):
    """ error-class for errors that are non-critical (computations can continue) """
    pass


class Critical_error(Exception):
    """ error-class for errors that are critical (computations should be stopped) """
    pass


class sbml_mca:
    """ Class providing MCA functionality """

    _time_variable = 'SBML_MCA_TIME_VARIABLE'
    # maximal number for product of parameternumber and species number (for
    # time-varying coefficients)
    _max_max_model_size = 2500
    _small_factor = 1e-12  # used to add to the trace of a singular jacobian

    def __init__(self, model):
        """
        @type model:  libsbml.model or string
        @param model: SBML model, or filename
        """
        if type(model) == str:
            self._doc = libsbml.readSBML(model)
            self._model = self._doc.getModel()
        else:
            self._model = model.clone()
            self._doc = libsbml.SBMLDocument(
                model.getLevel(), model.getVersion())
            self._doc.setModel(self._model)

        self._check_not_supported(self._model)
        misc.make_unique_local_parameters(self._model)
        self._species_changed_by_rule = self._get_species_changed_by_rule()
        self._species_ids = [s.getId() for s in filter(
            self._is_species_not_constant, self._model.getListOfSpecies())]
        self._species2pos = dict(
            zip(self._species_ids, range(self._species_ids.__len__())))
        self._N = self._get_stoich_mat(self._model)
        self._species_volume_conversion = self._get_species_volume_conversion()
        self._external_species_conc = self._get_external_species_conc(
            self._model)
        self._external_species_ids = [s.getId()
                                      for s in self._get_constant_species()]
        self._parameter_ids = None
        self._enzyme_positions = self._not_enzyme_positions = None
        # memorize parameters (e.g.) steady state conc. for the computation of
        # elasticities 
        self._reference_params = {}
        self._reset_calculated_parameters()
        self._build_math()

        self._no_steady_state = False
        [self._L_inv, self._L, self._Nr] = self._partition_stoich_matrix(self._N)
        # consider RuntimeWarnings (e.g. zero division) as errors
        warnings.filterwarnings("error", '.*', RuntimeWarning)

    def _is_species_not_constant(self, s):
        if s.getId() in self._species_changed_by_rule:
            return False
        return not(s.getConstant() or s.getBoundaryCondition())

    def _get_species_changed_by_rule(self):
        species_changed_by_rule = []
        for rule in self._model.getListOfRules():
            var = rule.getVariable()
            if self._model.getSpecies(var) is not None:
                species_changed_by_rule.append(var)
        return species_changed_by_rule

    def _build_math(self):
        """ construct the kinetic laws and rules """
        self._replacements = self._get_replacements(self._model)
        self._rate_rules = {}
        self._assignment_rules = {}
        self._handle_rules(self._model)
        self._kinetic_laws = self._get_kinetic_laws(self._model)
        self._rate_rules = self._replace_flux_symbols(self._rate_rules)

    def _reset_calculated_parameters(self):
        """ reset/initialize all parameters that are calculated and depend on initial coniditions / parameters """
        self._ss = self._ss_s0 = self._p_ela = self._e_ela = self._ccc = self._fcc = self._crc = self._frc = self._2nd_crc = self._2nd_frc = None

    def _get_stoich_mat(self, model):
        """ get the stoichiometric matrix (not including constant and boundary condition species) """
        N = numpy.zeros((self._species_ids.__len__(), model.getNumReactions()))
        for i, r in enumerate(model.getListOfReactions()):
            modes = [(-1, 'Reactants'), (1, 'Products')]
            for sign, direction in modes:
                for sr in getattr(r, 'getListOf' + direction)():
                    s = model.getSpecies(sr.getSpecies())
                    if s.getBoundaryCondition() \
                            or s.getConstant() \
                            or s.getId() in self._species_changed_by_rule:  # we dont want no constant species in our stoich matrix
                        continue
                    j = self._species2pos[sr.getSpecies()]
                    N[j, i] += sign * sr.getStoichiometry()
        if len(N) == 0:
            raise Critical_error('Empty stoichiometric matrix.')
        return N

    def _get_species_volume_conversion(self):
        comp_size = {comp.getId(): comp.getSize()
                     for comp in self._model.getListOfCompartments()}
        for c in comp_size:
            if math.isnan(comp_size[c]):
                comp_size[c] = 1.

        factors = []
        for s in [self._model.getSpecies(id) for id in self._species_ids]:
            if s.getHasOnlySubstanceUnits():
                factors.append(1.)
                continue
            factors.append(1. / comp_size[s.getCompartment()])
        return numpy.array(factors)

    def _get_kinetic_laws(self, model, compile_formulas=True):
        """ gather string representation of the kinetic laws """
        formulas = []
        self.kl = []
        for pos, kl in enumerate([r.getKineticLaw() for r in model.getListOfReactions()]):
            formula = self._ast_to_string(kl.getMath(), mode='python')
            formula = formula.replace(' ', '')
            # print model.getReaction(pos).getId(),  formula
            if compile_formulas:
                formula = compile(formula, 'kl_compile', 'eval')
            formulas.append(formula)
        return formulas

    def _get_replacements(self, model):
        """ get dictionary of parameter values and compartmets"""
        # do not take include parameters that are modified by a rule
        params_changed_by_rule = [r.getVariable()
                                  for r in self._model.getListOfRules()]
        params = {}
        for (pos, base) in enumerate([model] + [r.getKineticLaw() for r in model.getListOfReactions()]):
            for p in base.getListOfParameters():
                if not p.getId() in params_changed_by_rule:
                    params[p.getId()] = p.getValue()
        for comp in model.getListOfCompartments():
            s = 1.
            if comp.isSetSize():
                s = comp.getSize()
            elif comp.isSetVolume():
                s = comp.getVolume()
            params[comp.getId()] = s

        # handle initial assignments (they might be dependent on each other,
        # therefore try 50 evaluations)
        max_iter = 50
        # this has to be done manually, because get_initial_conc method depends
        # on parameters which are not yet known
        params_with_species = {}
        for s in model.getListOfSpecies():
            if s.isSetInitialConcentration():
                params_with_species[s.getId()] = s.getInitialConcentration()
            elif s.isSetInitialAmount():
                params_with_species[s.getId()] = s.getInitialAmount()

        while (any([math.isnan(params[x]) for x in params]) and max_iter > 0):
            params_with_species.update(params)
            ass = self._evaluate_initial_assignments(
                model, params_with_species)
            for p in ass:
                if params.has_key(p) and math.isnan(params[p]):
                    params[p] = float(ass[p])
                    # print p
                    # print ass[p]
                max_iter -= 1
        return params

    def _handle_rules(self, model):
        """ get expressions for rules in the model """
        # assignment rules have to be handeled before rate rules
        is_loop = True
        while is_loop:  # loop until no assignment rule is dependent on another assignment
            for rule in model.getListOfRules():
                if rule.isAssignment():
                    var = rule.getVariable()
                    formula = self._ast_to_string(
                        rule.getMath(), mode='', replace=True)
                    formula_wo_replace = self._ast_to_string(
                        rule.getMath(), mode='', replace=False)
                    self._assignment_rules[var] = {
                        True: formula, False: formula_wo_replace}
            # check dependencies
            is_loop = False
            for var1 in self._assignment_rules:
                for var2 in self._assignment_rules:
                    if var2 in self._assignment_rules[var1][True]:
                        is_loop = True

        for rule in model.getListOfRules():
            var = rule.getVariable()
            formula = self._ast_to_string(rule.getMath(), mode='')
            if rule.isRate():
                # print formula
                self._rate_rules[var] = formula
            elif rule.isAlgebraic():
                raise Critical_error('Algebraic rules not supported')

    def _compute_assignment_rules(self, simulation_result, species_ids):
        """ compute assignment rules for species after a simulation """
        simulation_result_dict = {}
        for pos, id in enumerate(species_ids):
            simulation_result_dict[id] = simulation_result[:, pos]
        result = {}
        for ar in self._assignment_rules:
            result[ar] = eval(self._assignment_rules[ar][True],
                              self._external_species_conc, simulation_result_dict)
        return result

    def _replace_flux_symbols(self, formula_dict):
        """ replace references to reaction fluxes in rate rules """
        r_ids = [r.getId() for r in self._model.getListOfReactions()]
        kinetic_laws = self._get_kinetic_laws(
            self._model, compile_formulas=False)
        for id in formula_dict:
            for (r_id, kl) in zip(r_ids, kinetic_laws):
                formula_dict[id] = formula_dict[id].replace(r_id, kl)
        return formula_dict

    def _evaluate_initial_assignments(self, model, params):
        params = params.copy()
        assignments = {}
        for ass in model.getListOfInitialAssignments():
            var = ass.getSymbol()
            formula = self._ast_to_string_libsbml(ass.getMath())
            # print formula
            assignments[var] = float(eval(formula, params))
        return assignments

    def _check_not_supported(self, model):
        """ check for features in the sbml that are not supported yet """
        if model.getNumConstraints():
            raise Critical_error('Error: Constraints not supported yet')
        if model.getNumEvents():
            raise Critical_error('Error: Events not supported yet')

        p_ids = [p.getId() for p in self._model.getListOfParameters()]
        for ia in self._model.getListOfInitialAssignments():
            if ia.getSymbol() not in p_ids:
                raise Critical_error(
                    'Initial assignments are currently only supported for parameters.')
        for c in model.getListOfCompartments():
            if not c.getConstant():
                raise Critical_error(
                    'Error: Varying compartment sizes not yet supported')

    def _get_enzyme_positions(self):
        """ get the positions of the enzymes in the list of (non constant) species """
        if self._enzyme_positions is not None:
            return self._enzyme_positions
        positions = []
        for pos, species in enumerate([self._model.getSpecies(id) for id in self._species_ids]):
            if species.getSBOTerm() == sbo_type2id['enzyme'] or species.getId().startswith('enzyme'):
                positions.append(pos)
        self._enzyme_positions = positions
        return positions

    def _get_not_enzyme_positions(self):
        """ get the positions of the species without the enzymes """
        if self._not_enzyme_positions is not None:
            return self._not_enzyme_positions
        ret = range(len(self._species_ids))
        for pos in self._get_enzyme_positions():
            ret.remove(pos)
        #[ret.remove(self._species_ids.index(s)) for s in self._get_constant_species()]
        self._not_enzyme_positions = ret
        return ret

    def _get_constant_species(self):
        """ get list of constant species
        @rtype: list of libsbml.species
        @return list of constant species
        """
        # return misc.get_constant_species(self._model)
        is_species_constant = lambda s: not(self._is_species_not_constant(
            s)) and s.getId() not in self._species_changed_by_rule
        return filter(is_species_constant, self._model.getListOfSpecies())

    def _get_not_constant_species(self):
        """ get list of not constant species """
        # return misc.get_not_constant_species(self._model)
        return filter(self._is_species_not_constant, self._model.getListOfSpecies())

    def _get_external_species_conc(self, model):
        """ get external metabolite concentrations as a dict """
        d = {}
        for s in self._get_constant_species():
            if s.isSetInitialConcentration():
                d[s.getId()] = s.getInitialConcentration()
            elif s.isSetInitialAmount():
                d[s.getId()] = s.getInitialAmount()
        return d

    def get_parameter_values(self, parameter_ids=None):
        """
        Get values for the specified parameters
        @type   parameter_ids: list
        @param  parameter_ids: List of strings with parameter ids
        @rtype  numpy.array
        @return array with parameter values
        """
        if parameter_ids is None:
            parameter_ids = self.get_parameter_ids()
        return numpy.array([misc.get_parameter_value(self._model, p) for p in parameter_ids])

    def set_parameter_values(self, parameter_names, parameter_values):
        """
        Set list of parameters to new values
        @type   parameter_names: list
        @param  parameter_names: List of strings with parameter ids
        @type   parameter_values: list
        @param  parameter_values: List of parameter values
        """
        rebuild = False
        for i, p in enumerate(parameter_names):
            if p in self._external_species_ids:
                self._external_species_conc[p] = parameter_values[i]
            else:
                misc.set_parameter_value(self._model, p, parameter_values[i])
                rebuild = True
        if rebuild:
            self._build_math()
        self._reset_calculated_parameters()

    def get_parameter_ids(self):
        """ get list of parameters that are varied in p_elasticities """
        if self._parameter_ids is None:
            const_species_ids = [s.getId()
                                 for s in self._get_constant_species()]
            params_changed_by_rule = [r.getVariable()
                                      for r in self._model.getListOfRules()]
            global_param_ids = []
            for p in self._model.getListOfParameters():
                # if p.getConstant():
                if p.getId() not in params_changed_by_rule:
                    global_param_ids.append(p.getId())
            local_param_ids = []
            for r in self._model.getListOfReactions():
                kl = r.getKineticLaw()
                for p in kl.getListOfParameters():
                    if p.getConstant():
                        local_param_ids.append(p.getId())
            self._parameter_ids = const_species_ids + global_param_ids + local_param_ids
            for p_name in self._parameter_ids:
                if self._parameter_ids.count(p_name) != 1:
                    raise Critical_error(
                        'Parameter ID %s used multiple times. This is valid but not yet supported.' % p_name)
        return self._parameter_ids

    def get_delta_parameters(self, d_param, parameter_names):
        """ enter a an array of parameter deviations and names and get back the corresponding d_param vector for all parameters
        @type    d_param:  numpy.array
        @param   d_param:  vector of parameter changes
        @type    parameter_names:  list of strings
        @param   parameter_names:  list of parameter names
        @rtype:            numpy.array
        @return:           vector of all parameter changes
        """
        all_p_names = self.get_parameter_ids()
        dp = numpy.zeros(len(all_p_names))
        for value, name in zip(d_param, parameter_names):
            dp[all_p_names.index(name)] = value
        return dp

    def get_initial_conc(self, with_rate_rule_params=False, with_initial_assignment=True):
        # TODO: implement initial assignment
        """
        get vector of initial concentrations
        @type    with_rate_rule_params: boolean
        @param   with_rate_rule_params: indicate whehter to include initial values for objects determined by rate rules
        @rtype:                         numpy.array
        @return:                        vector of initial concentrations
        """
        s0 = []
        #species_wo_initial = []
        for s in self._model.getListOfSpecies():
            # if s.getBoundaryCondition() or s.getConstant():
            if s.getConstant() or s.getBoundaryCondition() \
                    or s.getId() in self._assignment_rules.keys():
                continue
            if s.isSetInitialConcentration():
                s0.append(s.getInitialConcentration())
            elif s.isSetInitialAmount():
                s0.append(s.getInitialAmount())
            else:
                #species_wo_initial.append( s.getId() )
                s0.append(0.)
                sys.stderr.write(
                    'No initial value specified for species %s. Setting initial value to 0.\n' % s.getId())
                #raise Critical_error('No initial value specified for species %s' %s.getId())
        if with_rate_rule_params:
            for var in self._rate_rules:
                if not var in self._species_ids:
                    # s0.append(self._replacements[var])
                    s0.append(misc.get_parameter_value(self._model, var))
        return numpy.array(s0)

    def set_initial_conc(self, s0):
        """
        set the initial concentrations
        @type    s0:  numpy.array
        @param   s0:  vector of initial metabolite concentrations
        """
        pos = 0
        for s in self._model.getListOfSpecies():
            if s.getConstant() \
                    or s.getId() in self._assignment_rules.keys():
                continue
            s.setInitialConcentration(s0[pos])
            pos += 1

    def _ast_to_string_libsbml(self, ast):
        """ convert libsbml AST node to string, using libsbml """
        # TODO: replace _ast_to_string method with this one
        formula = libsbml.formulaToString(ast)
        for old, new in [('log', 'math.log')]:
            formula = formula.replace(old, new)
        return formula

    def _ast_to_string(self, ast, mode='python', replace=True):
        """ convert libsbml AST node to string """
        if not ast:
            return
        type = ast.getType()
        l = self._ast_to_string(ast.getLeftChild(), mode, replace)
        r = self._ast_to_string(ast.getRightChild(), mode, replace)
        if type == libsbml.AST_MINUS:
            if r is None:
                return '( - %s )' % l
            return '( %s  - %s )' % (l, r)
        elif type == libsbml.AST_PLUS:
            return '( %s  + %s )' % (l, r)
        elif type == libsbml.AST_TIMES:
            return '( %s  *  %s )' % (l, r)
        elif type == libsbml.AST_DIVIDE:
            return ' ( %s  /  %s ) ' % (l, r)
        elif type == libsbml.AST_FUNCTION_POWER:
            return '( %s  **  %s )' % (l, r)
        elif type == libsbml.AST_INTEGER:
            return str(ast.getInteger())
        elif type == libsbml.AST_REAL:
            return str(ast.getReal())
        elif type == libsbml.AST_REAL_E:
            return str(ast.getReal())
        elif type == libsbml.AST_CONSTANT_PI:
            return str(math.pi)
        elif type == libsbml.AST_CONSTANT_E:
            return str(math.e)
        elif type == libsbml.AST_CONSTANT_FALSE:
            return '0'
        elif type == libsbml.AST_CONSTANT_TRUE:
            return '1'
        elif type == libsbml.AST_FUNCTION_ROOT:
            return '( %s ** (1. / %s) )' % (r, l)
        elif type == libsbml.AST_FUNCTION_EXP:
            return 'math.exp(%s)' % l
        elif type == libsbml.AST_FUNCTION_LN:
            return 'math.log(%s)' % (l)
        elif type == libsbml.AST_FUNCTION_LOG:
            return 'math.log((%s),(%s))' % (l, r)
        elif type == libsbml.AST_FUNCTION_SIN:
            return 'math.sin(%s)' % l
        elif type == libsbml.AST_FUNCTION_SINH:
            return 'math.sinh(%s)' % l
        elif type == libsbml.AST_FUNCTION_TAN:
            return 'math.tan(%s)' % l
        elif type == libsbml.AST_FUNCTION_TANH:
            return 'math.tanh(%s)' % l
        elif type == libsbml.AST_FUNCTION_COS:
            return 'math.cos(%s)' % l
        elif type == libsbml.AST_FUNCTION_COSH:
            return 'math.cosh(%s)' % l
        elif type == libsbml.AST_RELATIONAL_EQ:
            return '( %s == %s )' % (l, r)
        elif type == libsbml.AST_RELATIONAL_GEQ:
            return '( %s >= %s )' % (l, r)
        elif type == libsbml.AST_RELATIONAL_GT:
            return '( %s > %s )' % (l, r)
        elif type == libsbml.AST_RELATIONAL_LEQ:
            return '( %s <= %s )' % (l, r)
        elif type == libsbml.AST_RELATIONAL_LT:
            return '( %s < %s )' % (l, r)
        elif type == libsbml.AST_RELATIONAL_NEQ:
            return '( %s != %s )' % (l, r)
        elif type == libsbml.AST_FUNCTION_PIECEWISE:
            if ast.getNumChildren() != 3:
                raise Critical_error(
                    'AST_FUNCTION_PIECEWISE not yet implemented completely')
            condition = self._ast_to_string(ast.getChild(1), mode, replace)
            return ' ( %s ) * ( %s ) + (1- %s ) *( %s )' % (condition, l, condition, r)
        elif type == libsbml.AST_FUNCTION_CEILING:
            return 'math.ceil( %s )' % l
        elif type == libsbml.AST_FUNCTION_FLOOR:
            return 'math.floor( %s )' % l
        elif type == libsbml.AST_FUNCTION_FACTORIAL:
            return 'math.factorial( %s )' % l
        elif type == libsbml.AST_RATIONAL:
            return str(ast.getReal())
        elif type == libsbml.AST_LOGICAL_AND:
            s = '(%s' % l
            for pos in range(1, ast.getNumChildren()):
                c = self._ast_to_string(ast.getChild(pos), mode, replace)
                s += '  and  ' + c
            s += ')'
            return s
        elif type == libsbml.AST_LOGICAL_OR:
            s = '( %s' % l
            for pos in range(1, ast.getNumChildren()):
                c = self._ast_to_string(ast.getChild(pos), mode, replace)
                s += '  or  ' + c
            s += ')'
            return s
        elif type == libsbml.AST_LOGICAL_XOR:
            s = '(bool( %s )' % l
            for pos in range(1, ast.getNumChildren()):
                c = self._ast_to_string(ast.getChild(pos), mode, replace)
                s += ' ^ bool( %s )' % c
            s += ')'
            return s
        elif type == libsbml.AST_NAME:
            name = ast.getName()
            if name in self._assignment_rules:  # assignment rules are always replaced
                return str(self._assignment_rules[name][replace])
            if replace:
                try:
                    return str(self._replacements[name])
                except:
                    pass
                try:
                    sp_pos = self._species2pos[name]
                    if mode == 'fortran':
                        return 'Y(%s)' % (sp_pos + 1)
                    else:
                        return str(name)
                except:
                    pass
            return name
        elif type == libsbml.AST_FUNCTION:
            children = [self._ast_to_string(ast.getChild(
                x), mode, replace) for x in range(ast.getNumChildren())]
            return self._get_formula(ast.getName(), children)
        elif type == libsbml.AST_NAME_TIME:
            return self._time_variable  # '$TIME$'
        elif type == libsbml.AST_POWER:
            raise Critical_error('AST_POWER not yet implemented')
        else:
            ast_name = misc.ast_code_to_string(type)
            print ast_name
            raise Critical_error(
                'mathematic fucntion %s not yet implemented' % ast_name)

    def _get_formula(self, name, params):
        """ get a string representation for a function definition with parameters already replaced """
        func = self._model.getFunctionDefinition(name)
        p_dict = dict(zip([func.getArgument(x).getName()
                           for x in range(func.getNumArguments())], params))
        # print p_dict
        # TODO: the formulaToString method does not work for functions such as log, exp, etc ...
        # TODO: unify the math conversion (decide whether to use the ast_to_string parser or the libsbml mehtod
        #formula = ' '+ libsbml.formulaToString( func.getBody() ) + ' '
        # for s in ['(',')','*','/','+','-',',']:
        #        formula = formula.replace(s,' %s '%s)
        formula = ' ' + self._ast_to_string(func.getBody()) + ' '
        # print formula
        for param in p_dict:
            formula = formula.replace(' ' + param + ' ', str(p_dict[param]))
        return formula

    def _make_species2conc(self, conc):
        """ construct a dictionary form the given concentration vector"""
        param_ids = []
        for var in self._rate_rules:
            if not var in self._species_ids:
                param_ids.append(var)
        return dict(zip(self._species_ids + param_ids, conc))

    def _dSdt(self, species_conc, t):
        """ compute rate of metabolite change """
        ret = self._species_volume_conversion * \
            numpy.dot(self._N, self._v(species_conc, t))
        # handle rate rules
        for var in self._rate_rules:
            #f = self._rate_rules[var].replace('$TIME$', str(t))
            f = self._rate_rules[var].replace(self._time_variable, str(t))
            species2conc = self._make_species2conc(species_conc)
            species2conc['math'] = globals()['math']
            #rate = eval( f, species2conc, self._external_species_conc )
            rate = eval(f, self._external_species_conc, species2conc)
            if self._species2pos.has_key(var):
                ret[self._species2pos[var]] = rate
            else:
                l = ret.tolist()
                l.append(rate)
                ret = numpy.array(l)
        return ret

    def _v(self, species_conc, t):
        """ get reaction velocities at given point """
        vec = []
        species2conc = self._make_species2conc(species_conc)
        species2conc[self._time_variable] = t
        species2conc['math'] = globals()['math']
        for func in self._kinetic_laws:
            vec.append(eval(func, self._external_species_conc, species2conc))
        return numpy.array(vec)

    def integrate(self, time, steps, s0=None, r_tol=1e-6, a_tol=1e-12):
        """
        integrate the model and return result as matrix
        @type      time:  number
        @param     time:  simulation time
        @type      steps: number
        @param     steps: number of time steps
        @type      s0:    numpy.array
        @param     s0:    vector of initial metabolite concentrations
        @type      r_tol: number
        @param     r_tol: relative tolervance
        @type      a_tol: number
        @param     a_tol: absolute tolerance
        @rtype     numpy.array
        @return    matrix containing a row for each metabolite time course
        """
        if s0 is None:
            s0 = self.get_initial_conc(with_rate_rule_params=True)
        t = numpy.linspace(0, time, steps)
        # return [ t, scipy.integrate.odepack.odeint( self._dSdt, s0, t,
        # rtol=r_tol, atol=a_tol ) ]
        return [t, self._odeint_wrapper(self._dSdt, s0, t, rtol=r_tol, atol=a_tol)]

    def integrate_with_constant_species(self, time, steps, s0=None, r_tol=1e-6, a_tol=1e-12):
        """
        same as integrate function, but also include species
        with boundary condition or constant set to true
        @type      time:  number
        @param     time:  simulation time
        @type      steps: number
        @param     steps: number of time steps
        @type      s0:    numpy.array
        @param     s0:    vector of initial metabolite concentrations
        @type      r_tol: number
        @param     r_tol: relative tolervance
        @type      a_tol: number
        @param     a_tol: absolute tolerance
        @rtype     numpy.array
        @return    matrix containing a row for each metabolite time course
        """
        t, trace = self.integrate(time, steps, s0, r_tol, a_tol)
        for s in self._get_constant_species():
            if s.isSetInitialConcentration():
                conc = s.getInitialConcentration()
            elif s.isSetInitialAmount():
                conc = s.getInitialAmount()
            elif self._assignment_rules.has_key(s.getId()):
                continue
            else:
                raise Noncritical_error(
                    'no initial value given for species %s' % (s.getId()))
            siz = trace.shape
            #new_siz =  ( siz[0], siz[1]+1 )
            #trace = numpy.resize( trace, new_siz )
            #trace[:,-1] = numpy.ones( siz[0] ) * conc
            #trace = numpy.append( trace, numpy.ones(siz[0]).T*conc, axis=0 )
            trace = numpy.insert(trace, siz[1], conc, axis=1)
        return t, trace

    def integrate_return_dict(self, time, steps, s0=None, r_tol=1e-6, a_tol=1e-12, with_constant_species=True, with_assignment_rules=True):
        """ integrates the model and returns a dictionary as result """
        s_ids = self._species_ids + \
            [p for p in self._rate_rules if not p in self._species_ids]
        if with_constant_species:
            t, trace = self.integrate_with_constant_species(
                time, steps, s0, r_tol, a_tol)
            s_ids = s_ids + [s.getId() for s in self._get_constant_species()]
        else:
            t, trace = self.integrate(time, steps, s0, r_tol, a_tol)

        result_d = {s: trace[:, pos] for pos, s in enumerate(s_ids)}
        result_d['time'] = t

        if with_assignment_rules:
            species_assigned = self._compute_assignment_rules(trace, s_ids)
            result_d.update(species_assigned)
        return result_d

    def plot_timecourse(self, time, steps, s0=None, r_tol=1e-6, a_tol=1e-12, legend=True):
        """
        integrate the model and return plot the result
        @type      time:   number
        @param     time:   simulation time
        @type      steps:  number
        @param     steps:  number of time steps
        @type      s0:     numpy.array
        @param     s0:     vector of initial metabolite concentrations
        @type      r_tol:  number
        @param     r_tol:  relative tolervance
        @type      a_tol:  number
        @param     a_tol:  absolute tolerance
        @type      legend: boolean
        @param     legend: indicate whether plot should include a legend
        """
        [t, result] = self.integrate(time, steps, s0, r_tol, a_tol)
        result = result.T
        for p in self._get_not_enzyme_positions() + range(len(self._species_ids), len(result)):
            pylab.plot(t, result[p])
        pylab.xlabel('time')
        pylab.ylabel('concentration')
        if legend:
            pylab.legend([self._species_ids[i] for i in self._get_not_enzyme_positions(
            )] + [p for p in self._rate_rules if not p in self._species_ids])
        pylab.show()

    def get_steady_state(self, s0=None, tol=1e-12, max_steps=100, max_computation_time=180, old_ss=None, d_params=None, old_conc_resp=None, string_output=False):
        """
        finda steady state of the system, if the steady state before a parameter perturbation (old_ss), the
        parameter change (d_params) and the response coefficients of the unperturbed model (old_conc_resp) are 
        given, the initial values are linearly approximated (can be much faster, than with initial values)
        @type   s0:        numpy.array
        @param  s0:        vector of initial concentrations
        @type   tol:       number
        @param  tol:       tolerance
        @type   old_ss:    numpy.array
        @param  old_ss:    previously computed steady stae
        @type   d_params:  numpy.array
        @param  d_params:  vector of parameter deviations
        @type   old_conc_resp: numpy.array
        @param  old_conc_resp: the concentration response coefficients of the model with the old parameter values
        @type   string_output: boolean
        @param  string_output: boolean indicating whether the output should be in string from
        @rtype:      numpy.array
        @return:     vector of steady state metabolite concentrations
        """
        if self._no_steady_state:
            raise Noncritical_error('No steady state could be found')
        if s0 is None:
            s0 = self.get_initial_conc(with_rate_rule_params=True)
        self._ss_s0 = s0
        if self._ss is not None and numpy.linalg.norm(self._ss_s0 - s0) < 0.0001:
            if string_output:
                return misc.matrix2string(numpy.array([self._ss]), self._species_ids)
            return self._ss
        # linear approximation of new ss
        if old_ss is not None and old_conc_resp is not None and d_params is not None:
            s0 = old_ss + numpy.dot(old_conc_resp, d_params)
        result = [s0]

        start_comp_time = time.time()
        sim_time = 0
        for i in range(max_steps):
            stepsize = 10 ** (i / 10 - 1)
            # print stepsize
            t = numpy.linspace(sim_time, sim_time + stepsize, 100)
            s0 = result[-1]
            try:
                result = self._odeint_wrapper(self._dSdt, s0, t)
            except Noncritical_error, e:
                if i > 0:
                    self._no_steady_state = True
                    raise Noncritical_error(
                        'No steady state could be found. This could indicate and oszillating model.')
                else:
                    raise e
            # print numpy.linalg.norm( result[0] - result[-1] )
            if (numpy.linalg.norm(result[0] - result[-1]) / s0.size) < tol:
                break
            sim_time += stepsize
            if i == (max_steps - 1) or (time.time() - start_comp_time > max_computation_time):
                self._no_steady_state = True
                raise Noncritical_error(
                    'No steady state could be found. This could indicate and oszillating model.')

        self._ss = result[-1]
        if string_output:
            return misc.matrix2string(numpy.array([result[1]]).T, ['Conc.'],  self._species_ids, justify='left')
        return result[1]

    def get_steady_state_flux(self, string_output=False):
        """
        get the flux in steady state
        @type       string_output:  boolean
        @param      string_output   boolean indicating whether the output should be a string
        @rtype                      numpy.array
        @return                     vector of reaction rates in steady state
        """
        ss_conc = self.get_steady_state()
        flux = self._v(ss_conc, 1)
        if string_output:
            return misc.matrix2string(numpy.array([flux]).T, ['Flux'], [r.getId() for r in self._model.getListOfReactions()], justify='left')
        return flux

    def solve_steady_state(self):
        """
        solve steady state by a numerical optimization
        @rtype:                 numpy.array
        @return:                vector of steady state concentrations
        """
        # optimization has to be done on the reduced system
        # TODO: implement different comp. sizes
        raise Critical_error(
            "This function does not yet work for differing compartment sizes")
        s0 = self.get_initial_conc()
        si = numpy.dot(self._L_inv, s0)
        t = s0 - numpy.dot(self._L, si)
        f = lambda x: numpy.linalg.norm(
            self._dSdt(numpy.dot(self._L, x) + t, 1))
        ss_i = scipy.optimize.fmin_bfgs(f, si)
        ss = numpy.dot(self._L, ss_i) + t
        return ss

    def _rref(self, mat_in, tol=1e-6):
        """ compute reduced row echelon form of a matrix """
        mat = mat_in.copy().astype(float)
        cols = mat[0].size
        rows = mat.T[0].size
        r = 0
        for c in range(cols):
            # find pivot row
            pos_max = abs(mat.T[c][r:]).argmax() + r
            elem_max = mat.T[c][pos_max]
            # zero row
            if abs(elem_max) < tol:
                #mat.T[c][r:] = zeros(rows-r)
                continue
            # swap
            [mat[r], mat[pos_max]] = [mat[pos_max].copy(), mat[r].copy()]
            # normalize
            mat[r] = mat[r] / elem_max
            # eliminate current column
            for row in range(rows):
                if row == r:
                    continue
                mat[row] = mat[row] - mat[row][c] * mat[r]
            r += 1
            if r == rows:
                break
        return mat

    def _partition_stoich_matrix2(self, mat):
        """ partition the stoichiometric matrix in link matrix and full rank matrix Nr """
        """ this method gave strange results on one model, therefore it was replaced by the mehtod below"""
        cols = mat[0].size
        rows = mat.T[0].size
        NI = numpy.concatenate((mat, numpy.eye(rows, rows)), 1)
        reduced = self._rref(NI)
        Nr_big = reduced.T[:cols].T
        M = reduced.T[cols:].T
        M_inv = numpy.linalg.inv(M)
        # eliminate zero rows
        for i, row in enumerate(Nr_big):
            if numpy.linalg.norm(row) < 1e-6:
                break
        Nr = Nr_big[:i]
        L = M_inv.T[:i].T
        L_inv = M[:i]
        return [L_inv, L, Nr]

    def _partition_stoich_matrix(self, mat):
        """ partition stoichiometric matrix """
        rref, pivot = sympy.Matrix(mat.T).rref(
        )  # compute reduced row echolon form to get the linear indep. rows
        Nr = mat[pivot]  # linear independent rows of N
        # link matrix is L = N*inv(Nr)  [because per definition N = L*Nr]
        L = numpy.dot(mat, numpy.linalg.pinv(Nr))
        try:
            L_inv = numpy.linalg.inv(L)  # inverse of link matrix
        except:
            L_inv = None
        return [L_inv, L, Nr]

    def _get_value_dict(self, species_conc):
        """ get dictionary with standard parameter values """
        value_dict = dict(zip(self._species_ids, species_conc.tolist()))
        value_dict.update(self._replacements)
        value_dict.update(self._external_species_conc)
        return value_dict

    def _get_elasticity(self, formula, value_dict, p_name):
        """ compute elasticity of formula w.r.t p_name """
        value_dict = value_dict.copy()
        x = float(value_dict[p_name])
        if abs(x) > 1e-5:
            dx = abs(x) * 0.001
        else:
            #dx = 0.0005
            dx = 1e-9
        y0 = eval(formula, value_dict, {p_name: x - dx, 'math': math})
        y1 = eval(formula, value_dict, {p_name: x + dx, 'math': math})
        dy = y1 - y0
        return dy / (2 * dx)

    def _get_2nd_elasticity(self, formula, value_dict, p_name1, p_name2):
        """ get the 2nd derivations of the input function with respect to p_name1 and p_name2 """
        value_dict = value_dict.copy()
        default_diff = 1e-5
        x = float(value_dict[p_name1])
        y = float(value_dict[p_name2])
        if abs(x) > 1e-5:
            dx = abs(x) * 0.0001
        else:
            dx = default_diff
        if abs(y) > 1e-5:
            dy = abs(y) * 0.0001
        else:
            dy = default_diff
        identical = (p_name1 == p_name2)
        a = [0.] * 4  # array to store function evaluations
        # list of direction where to move .. 1 step left and down, 1 step right
        # and down ...
        if identical:  # d^2 f/ dxdx
            directions = [-2 * dx, 0, 2 * dx]
            for pos, dx1 in enumerate(directions):
                a[pos] = eval(formula, value_dict, {p_name1: x + dx1, 'math': math})
            der = (a[0] - 2 * a[1] + a[2]) / (4 * dx * dx)
        else:         # d^2 f / dxdy
            directions = [(-dx, -dy), (dx, -dy), (-dx, dy), (dx, dy)]
            for pos, [dx1, dx2] in enumerate(directions):
                a[pos] = eval(formula, value_dict, {p_name1: x + dx1, p_name2: y + dx2, 'math': math})
            der = (a[0] - a[1] - a[2] + a[3]) / (4 * dx * dy)
        return der

    def get_2nd_elasticities(self, parameter_names=None, ss_conc=None):
        """
        get the second order elasticities, if parameter_names is not speciefied for all parameters
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @rtype:                    numpy.array
        @return:                   3d-tensor with second order elasticities, 1st index: reaction, 2nd index: param 1 3rd index: param 2
        """
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if parameter_names is None:
            parameter_names = self.get_parameter_ids()
        value_dict = self._get_value_dict(ss_conc)
        #ee = numpy.zeros( (self._model.getNumReactions(),len(p_names),len(p_names)) )
        ee = numpy.zeros((self._model.getNumReactions(), len(
            parameter_names), len(parameter_names)))
        for r_pos, kl in enumerate([r.getKineticLaw() for r in self._model.getListOfReactions()]):
            formula = self._ast_to_string(
                kl.getMath(), mode='python', replace=False)
            for p1_pos, p1_id in enumerate(parameter_names):
                for p2_pos, p2_id in enumerate(parameter_names[:p1_pos + 1]):
                    ela = self._get_2nd_elasticity(
                        formula, value_dict, p1_id, p2_id)
                    ee[r_pos, p1_pos, p2_pos] = ee[r_pos, p2_pos, p1_pos] = ela
                    #ee[p1_pos, r_pos, p2_pos] = ee[p2_pos, r_pos, p1_pos] = ela
        return ee

    def get_elasticities(self, parameter_names, ss_conc=None, string_output=False, time=None):
        """
        get elasticities for a custom set of parameters
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @rtype:                    numpy.array
        @return:                   matrix with elasticities
        """
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        value_dict = self._get_value_dict(ss_conc)
        if time is not None:
            value_dict[self._time_variable] = time
        ec = numpy.zeros((self._model.getNumReactions(), len(parameter_names)))

        for r_pos, kl in enumerate([r.getKineticLaw() for r in self._model.getListOfReactions()]):
            formula = self._ast_to_string(
                kl.getMath(), mode='python', replace=False)
            for p_pos, p_id in enumerate(parameter_names):
                if not p_id in formula:
                    continue
                e = self._get_elasticity(formula, value_dict, p_id)
                ec[r_pos, p_pos] = e
        if string_output:
            return misc.matrix2string(ec, parameter_names, [r.getId() for r in self._model.getListOfReactions()])
        return ec

    def get_parameter_elasticities(self, ss_conc=None, string_output=False, normalize=None):
        """
        get parameter elasticities ( dv / dp )
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @rtype:                    numpy.array
        @return:                   matrix with elasticities
        """
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None, self.get_parameter_ids()),
                                                       'both': ('flux', self.get_parameter_ids())}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._p_ela is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ep_ss'] - ss_conc) < 1e-8:
                    p_ela = self._p_ela
                    if normalize:
                        p_ela = normalize_wrap(p_ela)
                    if string_output:
                        return misc.matrix2string(p_ela,  self.get_parameter_ids(), [r.getId() for r in self._model.getListOfReactions()])
                    return p_ela
            except:
                pass
        p_names = self.get_parameter_ids()
        p_ela = self.get_elasticities(p_names, ss_conc)
        self._p_ela = p_ela
        if normalize:
            p_ela = normalize_wrap(p_ela, normalize)
        self._reference_params['ep_ss'] = ss_conc
        if string_output:
            return misc.matrix2string(p_ela, p_names, [r.getId() for r in self._model.getListOfReactions()])
        return p_ela

    def get_metabolite_elasticities(self, ss_conc=None, string_output=False, normalize=None):
        """
        get substrate elasticities ( dv / ds )
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @rtype:                    numpy.array
        @return:                   matrix with elasticities
        """
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None,  'conc'),
                                                       'both': ('flux', 'conc')}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._e_ela is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ee_ss'] - ss_conc) < 1e-8:
                    e_ela = self._e_ela
                    if normalize:
                        e_ela = normalize_wrap(e_ela)
                    if string_output:
                        return misc.matrix2string(e_ela, self._species_ids, [r.getId() for r in self._model.getListOfReactions()])
                    return e_ela
            except:
                pass
        e_ela = self.get_elasticities(self._species_ids, ss_conc)
        self._reference_params['ee_ss'] = ss_conc
        self._e_ela = e_ela
        if normalize:
            e_ela = normalize_wrap(e_ela, normalize)
        if string_output:
            return misc.matrix2string(e_ela, self._species_ids, [r.getId() for r in self._model.getListOfReactions()])
        return e_ela

    def get_flux_cc(self, ss_conc=None, string_output=False, normalize=None):
        """
        get flux control coefficients ( dJ / dv )
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing control coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None,  'flux'),
                                                       'both': ('flux', 'flux')}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._fcc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['fcc_ss'] - ss_conc) < 1e-8:
                    fcc = self._fcc
                    if normalize:
                        fcc = normalize_wrap(fcc)
                    if string_output:
                        return misc.matrix2string(fcc, self._species_ids, [r.getId() for r in self._model.getListOfReactions()])
                    return fcc.copy()
            except:
                pass
        e_ela = self.get_metabolite_elasticities(ss_conc)
        ccc = self.get_conc_cc(ss_conc)
        prod = numpy.dot(e_ela, ccc)
        fcc = numpy.eye(prod[0].size) + prod
        self._fcc = fcc
        if normalize:
            fcc = normalize_wrap(fcc, normalize)
        if string_output:
            return misc.matrix2string(fcc, [r.getId() for r in self._model.getListOfReactions()], [r.getId() for r in self._model.getListOfReactions()])
        return fcc

    def get_conc_cc(self, ss_conc=None, string_output=False, normalize=None):
        """
        get concentration control coefficients ( dS / dv )
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing control coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('conc', None),
                                                       'right': (None,  'flux'),
                                                       'both': ('conc', 'flux')}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._ccc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ccc_ss'] - ss_conc) < 1e-8:
                    ccc = self._ccc
                    if normalize:
                        ccc = normalize_wrap(ccc, normalize)
                    if string_output:
                        return misc.matrix2string(ccc, [r.getId() for r in self._model.getListOfReactions()], self._species_ids)
                    return ccc.copy()
            except:
                pass
        e_ela = self.get_metabolite_elasticities(ss_conc)
        M = numpy.dot(numpy.dot(self._Nr, e_ela), self._L)
        try:
            M_inv = numpy.linalg.inv(M)
        except:
            try:
                M_inv = numpy.linalg.inv(
                    M + numpy.eye(len(M)) * self._small_factor)
            except:
                raise Noncritical_error('Error: Singular Jacobian')

        ccc = numpy.dot(-numpy.dot(self._L, M_inv), self._Nr)
        self._ccc = ccc.copy()
        self._reference_params['ccc_ss'] = ss_conc
        if normalize:
            ccc = normalize_wrap(ccc, normalize)
        if string_output:
            return misc.matrix2string(ccc,  [r.getId() for r in self._model.getListOfReactions()], self._species_ids)
        return ccc

    def get_conc_resp(self, ss_conc=None, string_output=False, normalize=None):
        """
        get concentration response coefficients ( dS / dp)
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('conc', None),
                                                       'right': (None, self.get_parameter_ids()),
                                                       'both': ('conc', self.get_parameter_ids())}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._crc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['crc_ss'] - ss_conc) < 1e-8:
                    crc = self._crc
                    if normalize:
                        crc = normalize_wrap(crc, normalize)
                    if string_output:
                        return misc.matrix2string(crc, self.get_parameter_ids(), self._species_ids)
                    return crc
            except:
                pass
        p_ela = self.get_parameter_elasticities(ss_conc)
        ccc = self.get_conc_cc(ss_conc)
        crc = numpy.dot(ccc, p_ela)
        self._crc = crc.copy()
        self._reference_params['crc_ss'] = ss_conc
        if normalize:
            crc = normalize_wrap(crc, normalize)
        if string_output:
            return misc.matrix2string(crc, self.get_parameter_ids(), self._species_ids)
        return crc

    def get_flux_resp(self, ss_conc=None, string_output=False, normalize=None):
        """
        get flux response coefficients ( dJ / dp )
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None, self.get_parameter_ids()),
                                                       'both': ('flux', self.get_parameter_ids())}[normalize])
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        if self._frc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['frc_ss'] - ss_conc) < 1e-8:
                    frc = self._frc
                    if normalize:
                        frc = normalize_wrap(frc, normalize)

                    if string_output:
                        return misc.matrix2string(frc, self.get_parameter_ids(), [r.getId() for r in self._model.getListOfReactions()])
                    return frc
            except:
                pass
        p_ela = self.get_parameter_elasticities(ss_conc)
        fcc = self.get_flux_cc(ss_conc)
        frc = numpy.dot(fcc, p_ela)
        self._frc = frc.copy()
        self._reference_params['frc_ss'] = ss_conc
        if normalize:
            frc = normalize_wrap(frc, normalize)
        if string_output:
            return misc.matrix2string(frc, self.get_parameter_ids(), [r.getId() for r in self._model.getListOfReactions()])
        return frc

    def get_custom_conc_resp(self, parameter_names, ss_conc=None, string_output=False, normalize=None):
        """
        get concentration response coefficients for a custom set of parameters
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        custom_ela = self.get_elasticities(parameter_names, ss_conc)
        ccc = self.get_conc_cc(ss_conc)
        crc = numpy.dot(ccc, custom_ela)
        if string_output:
            return misc.matrix2string(crc, parameter_names, self._species_ids)
        if normalize:
            l, r = {'left': ('conc', None), 'right': (
                None, parameter_names), 'both': ('conc', parameter_names)}[normalize]
            crc = self._normalize_coefficients(crc, left=l, right=r)
        return crc

    def get_custom_flux_resp(self, parameter_names, ss_conc=None, string_output=False, normalize=None):
        """
        get flux response coefficients for a custom set of parameters
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @type    string_output:    boolean
        @param   string_output:    boolean indicating whether the output should be in string form
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   matrix containing response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        custom_ela = self.get_elasticities(parameter_names, ss_conc)
        fcc = self.get_flux_cc(ss_conc)
        frc = numpy.dot(fcc, custom_ela)
        if string_output:
            return misc.matrix2string(frc, parameter_names, self._species_ids)
        if normalize:
            l, r = {'left': ('flux', None), 'right': (
                None, parameter_names), 'both': ('flux', parameter_names)}[normalize]
            frc = self._normalize_coefficients(frc, left=l, right=r)
        return frc

    def _normalize_coefficients(self, coeff, left=None, right=None, first_order=None):
        """ normalize coefficients (2d or 3d) """
        min_value_norm = 1e-18  # minimal value that is divided by in normalization

        if left is not None:
            if left == 'conc':
                d = self.get_steady_state()
            elif left == 'flux':
                d = self._v(self.get_steady_state(), 1)
            elif isinstance(left, list):
                raise Exception('Impelment me...')
            elif isinstance(left, numpy.ndarray):
                d = left
            else:
                raise Exception('Unkown input for normalization')

            # for values very close to zero, we want inf/nan values for the
            # normalization result
            d[abs(d) < min_value_norm] = 0.
            with numpy.errstate(divide='print'):
                L = numpy.diag(1. / d)
                print L
            if coeff.shape.__len__() == 2:
                coeff = numpy.dot(L, coeff)
            if coeff.shape.__len__() == 3:
                coeff = numpy.tensordot(L, coeff, [1, 0])

        if right is not None:
            if isinstance(right, list):
                ss = self.get_steady_state()
                value_d = self._get_value_dict(ss)
                d = numpy.array([value_d[x] for x in right])
            elif isinstance(right, numpy.ndarray):
                d = right
            elif right == 'conc':
                d = self.get_steady_state()
            elif right == 'flux':
                d = self._v(self.get_steady_state(), 1)
            else:
                raise Exception('Unkown input for normalization')
            #d[abs(d)<min_value_norm] = 0.
            R = numpy.diag(d)
            if coeff.shape.__len__() == 2:
                coeff = numpy.dot(coeff, R)
            if coeff.shape.__len__() == 3:
                coeff = numpy.tensordot(
                    numpy.tensordot(coeff, R, [2, 0]), R, [1, 0])

        if first_order is not None:
            # correction terms for normalization of second order terms
            if not isinstance(first_order, numpy.ndarray):
                raise Exception('Not supported yet')

            coeff_correction_r = numpy.zeros(coeff.shape)
            coeff_correction_l = numpy.zeros(coeff.shape)
            if right is not None:
                for i in range(coeff_correction_r.shape[0]):
                    coeff_correction_r[i] = numpy.diag(first_order[i])
            if left is not None:
                coeff_correction_l = misc.matrix2tensor(
                    first_order, first_order)
            coeff = coeff + coeff_correction_r - coeff_correction_l
        if coeff.shape.__len__() == 3 and (left is not None or right is not None) and first_order is None:
            raise Exception('No first order term specified for normalization')
        return coeff

    def _normalize_time_varying_coeff(self, coeff, left=None, right=None):
        coeff_copy = coeff.copy()
        no_coeff = coeff.shape[1]  # number of coeff
        if left is not None:
            # number of right coefficients (here: parameters)
            no_right = no_coeff / left.shape[1]
            if right is not None:
                assert no_right == len(right)
            for i in range(no_coeff):
                try:
                    # divide earch row by left coeff (flux or conc.)
                    coeff_copy[:, i] = coeff_copy[
                        :, i] / left[:, int(i / no_right)]
                except RuntimeWarning:
                    raise Noncritical_error(
                        'Error: Normalization failed. Value to divide by is too small.')
        if right is not None:
            no_right = len(right)
            if left is not None:
                assert no_right == no_coeff / left.shape[1]
            for i in range(no_coeff):
                coeff_copy[:, i] = coeff_copy[:, i] * right[i % no_right]
        return coeff_copy

    def _compute_2nd_resp(self, custom_params=None):
        """ compute the second order responce coefficients """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        parameter_names = custom_params or self.get_parameter_ids()
        p_set = self._species_ids + parameter_names
        ee = self.get_2nd_elasticities(p_set)
        cs = self.get_conc_cc()
        cj = self.get_flux_cc()
        #rs = self.get_conc_resp()
        rs = self.get_custom_conc_resp(parameter_names)
        ep = self.get_elasticities(parameter_names)
        #ep = self.get_parameter_elasticities()
        es = self.get_metabolite_elasticities()
        size_s = len(self._species_ids)
        ess = ee[:, :size_s, :size_s]
        epp = ee[:, size_s:, size_s:]
        esp = ee[:, :size_s, size_s:]
        prod = numpy.tensordot
        gamma = prod(prod(ess, rs, [2, 0]), rs, [1, 0]) \
            + prod(esp, rs, [1, 0]) \
            + prod(esp, rs, [1, 0]).transpose([0, 2, 1]) \
            + epp
        rs2 = prod(cs, gamma, [1, 0])
        rj2 = prod(cj, gamma, [1, 0])
        if not custom_params:
            self._2nd_crc = rs2
            self._2nd_frc = rj2
        return (rs2, rj2)

    def get_2nd_conc_resp(self, normalize=None):
        """
        get the 2nd order concentration response coefficients
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   3d tensor containing second order concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if self._2nd_crc is None:
            self._compute_2nd_resp()
        if normalize:
            fst_order = self.get_conc_resp(normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('conc', None, fst_order),
                                                           'right': (None, self.get_parameter_ids(), fst_order),
                                                           'both': ('conc', self.get_parameter_ids(), fst_order)}[normalize])
            return normalize_wrap(self._2nd_crc, normalize)
        return self._2nd_crc

    def get_2nd_flux_resp(self, normalize=None):
        """
        get the 2nd order flux response coefficients
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   3d tensor containing second order flux response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if self._2nd_frc is None:
            self._compute_2nd_resp()
        if normalize:
            fst_order = self.get_flux_resp(normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('flux', None, fst_order),
                                                           'right': (None, self.get_parameter_ids(), fst_order),
                                                           'both': ('flux', self.get_parameter_ids(), fst_order)}[normalize])
            return normalize_wrap(self._2nd_frc, normalize)
        return self._2nd_frc

    def get_2nd_custom_conc_resp(self, parameter_names, normalize=None):
        """
        get the 2nd order concentration response coefficients for a custom parameter set
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   3d tensor containing second order concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        r2s = self._compute_2nd_resp(parameter_names)[0]
        if normalize:
            fst_order = self.get_custom_conc_resp(
                parameter_names, normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('conc', None, fst_order),
                                                           'right': (None, parameter_names, fst_order),
                                                           'both': ('conc', parameter_names, fst_order)}[normalize])
            return normalize_wrap(r2s, normalize)
        return r2s

    def get_2nd_custom_flux_resp(self, parameter_names, normalize=None):
        """
        get the 2nd order flux response coefficients for a custom parameter set
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   3d tensor containing second order flux response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        rj2 = self._compute_2nd_resp(parameter_names)[1]
        if normalize:
            fst_order = self.get_custom_flux_resp(
                parameter_names, normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('flux', None, fst_order),
                                                           'right': (None, parameter_names, fst_order),
                                                           'both': ('flux', parameter_names, fst_order)}[normalize])
            return normalize_wrap(rj2, normalize)
        return rj2

    def get_2nd_custom_resp(self, parameter_names):
        """
        returns a tuple with ( R2S, R2J ) (much faster than calling both custom methods...)
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @rtype:                    2d tuple of numpy.array
        @return:                   tuple containing second order concentration and flux response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        return self._compute_2nd_resp(parameter_names)

    def get_spectral_conc_cc(self, frequency, ss_conc=None):
        """
        get spectral concentration control coefficients
        @type    frequency:      number
        @param   frequency:      frequency of the perturbation
        @rtype:                  numpy.array
        @return:                 matrix of specetral concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        e_ela = self.get_metabolite_elasticities(ss_conc)
        M = numpy.dot(numpy.dot(self._Nr, e_ela), self._L) - \
            (1j * frequency * numpy.eye(len(self._Nr)))
        try:
            M_inv = numpy.linalg.inv(M)
        except:
            try:
                M_inv = numpy.linalg.inv(
                    M + numpy.eye(len(M)) * self._small_factor)
            except:
                raise Noncritical_error('Error: Singular Jacobian')

        s_ccc = numpy.dot(-numpy.dot(self._L, M_inv), self._Nr)
        return s_ccc

    def get_spectral_flux_cc(self, frequency, ss_conc=None):
        """
        get spectral flux control coefficients 
        @type    frequency:      number
        @param   frequency:      frequency of the perturbation
        @rtype:                  numpy.array
        @return:                 matrix of specetral flux response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        s_ccc = self.get_spectral_conc_cc(frequency, ss_conc)
        e_ela = self.get_metabolite_elasticities(ss_conc)
        prod = numpy.dot(e_ela, s_ccc)
        s_fcc = numpy.eye(prod[0].size) + prod
        return s_fcc

    def get_spectral_conc_resp(self, frequency, ss_conc=None):
        """
        get spectral concentration response coefficients 
        @type    frequency:      number
        @param   frequency:      frequency of the perturbation
        @rtype:                  numpy.array
        @return:                 matrix of specetral concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        s_ccc = self.get_spectral_conc_cc(frequency, ss_conc)
        p_ela = self.get_parameter_elasticities(ss_conc)
        s_crc = numpy.dot(s_ccc, p_ela)
        return s_crc

    def get_spectral_flux_resp(self, frequency, ss_conc=None):
        """
        get spectral flux response coefficients
        @type    frequency:      number
        @param   frequency:      frequency of the perturbation
        @rtype:                  numpy.array
        @return:                 matrix of specetral flux response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.get_steady_state()
        s_fcc = self.get_spectral_flux_cc(frequency, ss_conc)
        p_ela = self.get_parameter_elasticities(ss_conc)
        s_frc = numpy.dot(s_fcc, p_ela)
        return s_frc

    def get_jacobian(self, ss_conc=None):
        """
        get the jacobian (dS/dS)
        @type    ss_conc:        numpy.array
        @param   ss_conc:        vector of steady state metabolite concentrations
        @rtype:                  numpy.array
        @return:                 jacobian matrix
        """
        e_ela = self.get_metabolite_elasticities(ss_conc)
        return numpy.dot(numpy.dot(self._Nr, e_ela), self._L)

    def get_stability(self, ss_conc):
        """
        get the stability (eigenvalues jacobian) of the system at this point
        @type    ss_conc:        numpy.array
        @param   ss_conc:        vector of steady state metabolite concentrations
        @rtype:                  list
        @return:                 list of strings containing messages
        """
        M = self.get_jacobian(ss_conc)
        e, v = numpy.linalg.eig(M)
        msgs = []
        if (e > 0).any():
            msgs.append('Steady state is not stable.')
        else:
            msgs.append('Steady state is stable.')
        if (abs(numpy.imag(e)) > 1e-14).any():
            msgs.append('Steady state might oscillate')
        return msgs

    def _odeint_wrapper(self, d_func, x0, timepoints, rtol=1e-6, atol=1e-12):
        """ wrapper function for scipy.integrate.odepack.odeint adding error handling """
        errors = []

        def dy_error_handling(x, t):
            try:
                return d_func(x, t)
            except Exception, e:
                errors.append(e)
        max_steps = 1000
        max_stepsize = (timepoints[-1] - timepoints[0]) / max_steps
        result, infodict = scipy.integrate.odepack.odeint(dy_error_handling,
                                                          x0,
                                                          timepoints,
                                                          rtol=rtol,
                                                          atol=atol,
                                                          hmax=max_stepsize,
                                                          mxstep=max_steps,
                                                          full_output=True)
        if errors != []:
            # print errors
            #raise errors[0]
            raise Noncritical_error(
                'There was an error integrating the ODEs (' + errors[0].message + ')')
        if infodict['message'] != 'Integration successful.':
            # print infodict
            raise Noncritical_error(
                'There were numerial problems simulating the model.')
        return result

    def __getstate__(self):
        """ flatten object for pickeling """
        d = self.__dict__.copy()
        del d['_model']
        del d['_kinetic_laws']
        d['_doc'] = '<?xml version="1.0" encoding="UTF-8"?>' + self._doc.toSBML()
        return d

    def __setstate__(self, d):
        """ recover the model libsbml element after pickling """
        self.__dict__ = d
        self._doc = libsbml.readSBMLFromString(self._doc)
        self._model = self._doc.getModel()
        misc.make_unique_local_parameters(self._model)
        self._kinetic_laws = self._get_kinetic_laws(self._model)

    def get_elasticities_symbolic(self, parameter_names):
        # TODO: maybe evaluate first parameters / external conc.
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        ec = sympy.zeros((self._model.getNumReactions(), len(parameter_names)))

        for r_pos, kl in enumerate([r.getKineticLaw() for r in self._model.getListOfReactions()]):
            formula = self._ast_to_string(
                kl.getMath(), mode='python', replace=False)
            formula.replace('$TIME$', 'current_time_variable')
            for p_pos, p_id in enumerate(parameter_names):
                if not p_id in formula:
                    continue
                ec[r_pos, p_pos] = sympy.diff(formula, p_id)
        return ec

    def get_time_varying_conc_rc_numerical(self,
                                           end_time,
                                           normalize=None,
                                           initial_cond_as_params=False,
                                           return_multiple=False,
                                           return_species_tc=False):
        """
        get time varying concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        p_id = self.get_parameter_ids()
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self._species_ids)
        l_s = len(self._species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients

        if l_sp > self._max_max_model_size:
            raise Noncritical_error(
                'Model is too large to compute time varying coefficients.')

        def my_ode(x, t):
            # generate ODE system to solve for time varying conc. response
            # coefficients
            # the first l_sp values are the time varying r.c., here they are
            # reshaped to matrix form
            dc = x[:l_sp].reshape((l_s, l_p))
            s = x[l_sp:]  # the species concentrations at current time point
            # species elasticities dv/dS
            es = self.get_elasticities(self._species_ids, s, time=t)
            # parameter elasticities dv/dp
            ep = self.get_elasticities(p_id, s, time=t)
            if initial_cond_as_params:
                # elasticities for initial conditions are zero
                ep = numpy.hstack(
                    (ep, numpy.zeros((self._model.getNumReactions(), l_s))))
            dx = numpy.zeros(l_sp + len(self._species_ids))
            # compuation: d/dt R = N( dv/ds * ds/dp + dv/dp )
            dx[:l_sp] = numpy.reshape(
                numpy.dot(self._N, (numpy.dot(es, dc) + ep)), -1)
            dx[l_sp:] = self._dSdt(s, t)
            return dx

        s0 = self.get_initial_conc(with_rate_rule_params=True)
        c0 = numpy.zeros(l_sp)
        if initial_cond_as_params:
            # initialize derivatives w.r.t. initial conditions to 1.
            c0 = numpy.zeros((l_s, l_p))
            c0[:, len(p_id):] = numpy.eye(l_s)
            c0 = numpy.reshape(c0, -1)
        t = numpy.linspace(0, end_time, 100)
        result = self._odeint_wrapper(my_ode, numpy.concatenate((c0, s0)), t)
        tv_rc = result[:, :l_sp]

        errors = []
        try:
            tv_rc_not_norm = tv_rc.copy()
            if normalize in ['left', 'both']:
                s_t = result[:, l_sp:]
                try:
                    for i in range(l_sp):
                        tv_rc[:, i] = tv_rc[:, i] / s_t[:, int(i / l_p)]
                except RuntimeWarning:
                    #raise Noncritical_error( 'There was an error during normalization.' )
                    raise Noncritical_error(
                        'Error: Normalization failed. Value to divide by is too small.')
            if normalize in ['right', 'both']:
                p_v = self.get_parameter_values(p_id)
                if initial_cond_as_params:
                    p_v = numpy.concatenate((p_v, self.get_initial_conc()))
                for i in range(l_sp):
                    tv_rc[:, i] = tv_rc[:, i] * p_v[i % l_p]
        except Exception, e:
            errors.append(e)
            tv_rc = numpy.zeros(tv_rc_not_norm.shape)

        if return_multiple:
            return [t, tv_rc_not_norm, tv_rc, errors]

        if return_species_tc:
            species_tc = result[:, l_sp:]
            return [t, tv_rc, species_tc]

        return [t, tv_rc]

    def get_time_varying_conc_rc(self,
                                 end_time,
                                 normalize=None,
                                 initial_cond_as_params=False,
                                 return_flat=False,
                                 return_species_tc=False):
        """
        get time varying concentration response coefficients
        """
        if self._rate_rules != {}:
            raise Noncritical_error(
                'MCA methods are not avaiable for explicit ODE systems.')
        p_id = self.get_parameter_ids()
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self._species_ids)
        l_s = len(self._species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients

        if l_sp > self._max_max_model_size:
            raise Noncritical_error(
                'Model is too large to compute time varying coefficients.')

        def my_ode(x, t):
            # generate ODE system to solve for time varying conc. response
            # coefficients
            # the first l_sp values are the time varying r.c., here they are
            # reshaped to matrix form
            dc = x[:l_sp].reshape((l_s, l_p))
            s = x[l_sp:]  # the species concentrations at current time point
            # species elasticities dv/dS
            es = self.get_elasticities(self._species_ids, s, time=t)
            # parameter elasticities dv/dp
            ep = self.get_elasticities(p_id, s, time=t)
            if initial_cond_as_params:
                # elasticities for initial conditions are zero
                ep = numpy.hstack(
                    (ep, numpy.zeros((self._model.getNumReactions(), l_s))))
            dx = numpy.zeros(l_sp + len(self._species_ids))
            # compuation: d/dt R = N( dv/ds * ds/dp + dv/dp )
            dx[:l_sp] = numpy.reshape(
                numpy.dot(self._N, (numpy.dot(es, dc) + ep)), -1)
            dx[l_sp:] = self._dSdt(s, t)
            return dx

        s0 = self.get_initial_conc(with_rate_rule_params=True)
        c0 = numpy.zeros(l_sp)
        if initial_cond_as_params:
            # initialize derivatives w.r.t. initial conditions to 1.
            c0 = numpy.zeros((l_s, l_p))
            c0[:, len(p_id):] = numpy.eye(l_s)
            c0 = numpy.reshape(c0, -1)
        t = numpy.linspace(0, end_time, 100)
        result = self._odeint_wrapper(my_ode, numpy.concatenate((c0, s0)), t)
        tv_rc = result[:, :l_sp]

        if normalize is not None:
            left = right = None
            if normalize in ['left', 'both']:
                left = result[:, l_sp:]
            if normalize in ['right', 'both']:
                right = self.get_parameter_values(p_id)
                if initial_cond_as_params:
                    right = numpy.concatenate((right, self.get_initial_conc()))
            tv_rc = self._normalize_time_varying_coeff(tv_rc, left, right)

        if not return_flat:
            tv_rc = tv_rc.reshape((len(t), l_s, l_p))

        if return_species_tc:
            # return also species timecourse (required for computation of
            # flux-resp-coeff)
            species_tc = result[:, l_sp:]
            return [t, tv_rc, species_tc]

        return [t, tv_rc]

    def get_time_varying_flux_rc(self,
                                 end_time,
                                 normalize=None,
                                 initial_cond_as_params=False,
                                 return_flat=False,
                                 return_species_tc=False):
        """
        get time varying flux response coefficients
        """
        p_id = self.get_parameter_ids()
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self._species_ids)
        l_s = len(self._species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients
        l_v = self._N.shape[1]  # nubmer of fluxes
        timepoints, tv_rcs, species_tc = self.get_time_varying_conc_rc(end_time,
                                                                       normalize=None,
                                                                       initial_cond_as_params=initial_cond_as_params,
                                                                       return_flat=True,
                                                                       return_species_tc=True)
        # initialize flux coeff.
        tv_rc_j = numpy.zeros((len(timepoints), l_v * l_p))

        for i, t in enumerate(timepoints):
            s = species_tc[i, :]
            # species elasticities dv/dS
            es = self.get_elasticities(self._species_ids, s, time=t)
            # parameter elasticities dv/dp
            ep = self.get_elasticities(p_id, s, time=t)
            RS = tv_rcs[i, :].reshape((l_s, l_p))
            RJ = numpy.dot(es, RS) + ep
            tv_rc_j[i, :] = RJ.reshape(-1)

        if normalize is not None:
            left = right = None
            if normalize in ['left', 'both']:
                left = numpy.zeros((len(timepoints), l_v))
                for i, t in enumerate(timepoints):
                    left[i, :] = self._v(species_tc[i, :], t)
            if normalize in ['right', 'both']:
                right = self.get_parameter_values(p_id)
                if initial_cond_as_params:
                    right = numpy.concatenate((right, self.get_initial_conc()))
            tv_rc_j = self._normalize_time_varying_coeff(tv_rc_j, left, right)

        if not return_flat:
            tv_rc_j = tv_rc_j.reshape((len(timepoints), l_v, l_p))

        if return_species_tc:
            # return also species timecourse (required for computation of
            # flux-resp-coeff)
            return [timepoints, tv_rc_j, species_tc]

        return [timepoints, tv_rc_j]

    def get_time_varying_rc_web_if(self,
                                   end_time,
                                   what='conc',
                                   initial_cond_as_params=False):
        errors = []
        if what == 'conc':
            timepoints, tv_rc, species_tc = self.get_time_varying_conc_rc(end_time,
                                                                          normalize=None,
                                                                          initial_cond_as_params=initial_cond_as_params,
                                                                          return_flat=True,
                                                                          return_species_tc=True)
            left = species_tc
        elif what == 'flux':
            timepoints, tv_rc, specie_tc = self.get_time_varying_flux_rc(end_time,
                                                                         normalize=None,
                                                                         initial_cond_as_params=initial_cond_as_params,
                                                                         return_flat=True,
                                                                         return_species_tc=True)
            left = numpy.zeros((len(timepoints), l_v))
            for i, t in enumerate(timepoints):
                left[i, :] = self._v(species_tc[i, :], t)

        right = self.get_parameter_values()
        if initial_cond_as_params:
            right = numpy.concatenate((right, self.get_initial_conc()))

        try:
            tv_rc_norm = self._normalize_time_varying_coeff(tv_rc, left, right)
        except Exception, e:
            errors.append(e)
            tv_rc_norm = numpy.zeros(tv_rc.shape)

        return [timepoints, tv_rc, tv_rc_norm, errors]

    def plot_time_varying_coefficients(self, end_time, what='conc_resp', initial_cond_as_params=False, normalize=None):
        if what == 'conc_resp':
            timepoints, tvrc = self.get_time_varying_conc_rc(end_time,
                                                             normalize=normalize,
                                                             initial_cond_as_params=initial_cond_as_params,
                                                             return_flat=True)
            p_ids = self.get_parameter_ids()
            if initial_cond_as_params:  # if initial conditions are also condidered as parameters
                l_p += len(self._species_ids)
            r_ids = [r.getId() for r in self._model.getListOfReactions()]
            s_ids = self._species_ids
            coeff_names = numpy.reshape(
                [[s_id + '_' + p_id for s_id in s_ids] for p_id in p_ids], -1)
        elif what == 'flux_resp':
            timepoints, tvrc = self.get_time_varying_flux_rc(end_time,
                                                             normalize=normalize,
                                                             initial_cond_as_params=initial_cond_as_params,
                                                             return_flat=True)
            p_ids = self.get_parameter_ids()
            if initial_cond_as_params:  # if initial conditions are also condidered as parameters
                l_p += len(self._species_ids)
            r_ids = [r.getId() for r in self._model.getListOfReactions()]
            coeff_names = numpy.reshape(
                [[r_id + '_' + p_id for r_id in r_ids] for p_id in p_ids], -1)
        pylab.plot(timepoints, tvrc)
        pylab.legend(coeff_names)
        pylab.show()

    # from here on model information methods:

    def get_model_id(self):
        return self._model.getId()

    def get_model_name(self):
        return self._model.getName()

    def get_species_ids(self, include_constant=False):
        """ get species ids of SBML model
            in contrast to self._species_ids this function also returns constant species
        """
        if include_constant:
            return [s.getId() for s in self._model.getListOfSpecies()]
        else:
            return self._species_ids

    def get_species_names(self, include_constant=False, replace_empty_names_with_id=False):
        """ get species names of SBML model
        include_constant specifies whether constant species are also returned
        replace_empty_names_with_id specifies whether empty species names are replaced by the species ID 
        """
        if replace_empty_names_with_id:
            return [self._model.getSpecies(s).getName() or s
                    for s in self.get_species_ids(include_constant)]
        else:
            return [self._model.getSpecies(s).getName()
                    for s in self.get_species_ids(include_constant)]

    def get_species_initial_concentrations(self, include_constant=False):
        if include_constant:
            return [s.getInitialConcentration() for s in self._model.getListOfSpecies()]
        else:
            return [s.getInitialConcentration() for s in misc.get_not_constant_species(self._model)]

    def get_species_initial_amounts(self, include_constant=False):
        if include_constant:
            return [s.getInitialAmount() for s in self._model.getListOfSpecies()]
        else:
            return [s.getInitialAmount() for s in misc.get_not_constant_species(self._model)]

    def get_parameter_names(self, parameter_ids=None, replace_empty_names_with_id=False):
        if parameter_ids is None:
            parameter_ids = self.get_parameter_ids()
        if replace_empty_names_with_id:
            return [misc.get_parameter_name(self._model, p_id) or p_id
                    for p_id in parameter_ids]
        else:
            return [misc.get_parameter_name(self._model, p_id)
                    for p_id in parameter_ids]

    def get_reaction_ids(self, replace_empty_names_with_id=False):
        return [r.getId() for r in self._model.getListOfReactions()]

    def get_reaction_names(self, replace_empty_names_with_id=False):
        if replace_empty_names_with_id:
            return [r.getName() or r.getId()
                    for r in self._model.getListOfReactions()]
        else:
            return [r.getName() for r in self._model.getListOfReactions()]

    def get_reaction_substrates(self):
        return [[s.getSpecies() for s in r.getListOfReactants()] for r in self._model.getListOfReactions()]

    def get_reaction_products(self):
        return [[s.getSpecies() for s in r.getListOfProducts()] for r in self._model.getListOfReactions()]

    def get_reaction_modifiers(self):
        return [[s.getSpecies() for s in r.getListOfProducts()] for r in self._model.getListOfReactions()]

    def get_reaction_kinetics(self):
        return [r.getKineticLaw().getFormula() for r in self._model.getListOfReactions()]

    def get_rule_ids(self):
        return [r.getId() for r in self._model.getListOfRules()]


if __name__ == '__main__':
    import optparse
    import sys

    parser = optparse.OptionParser(usage='usage: %prog [options] model')
    parser.add_option('--string_output', action='store_true',
                      dest='string_output', help='format output as string')
    parser.add_option('-t', '--timecourse', metavar='time',
                      dest='timecourse', help='simulate model for given time')
    parser.add_option('--ss', action='store_true', dest='ss',
                      help='get steady state concentrations and flux')
    parser.add_option('--normalize', dest='normalize',
                      help='normalize coefficients. left | right| both')
    parser.add_option('--p_ela', action='store_true',
                      dest='p_ela', help='compute parameter elasticities')
    parser.add_option('--e_ela', action='store_true',
                      dest='e_ela', help='compute substrate elasticities')
    parser.add_option('--p_ela2', action='store_true', dest='p_ela2',
                      help='compute second order parameter elasticities')
    parser.add_option('--flux_cc', action='store_true',
                      dest='flux_cc', help='compute flux conctrol coefficients')
    parser.add_option('--conc_cc', action='store_true', dest='conc_cc',
                      help='compute concentration conctrol coefficients')
    parser.add_option('--flux_resp', action='store_true',
                      dest='flux_resp', help='compute flux response coefficients')
    parser.add_option('--conc_resp', action='store_true', dest='conc_resp',
                      help='compute concentration response coefficients')
    parser.add_option('--flux_resp2', action='store_true', dest='flux_resp2',
                      help='compute second order flux response coefficients')
    parser.add_option('--conc_resp2', action='store_true', dest='conc_resp2',
                      help='compute second order concentration response coefficients')
    parser.add_option('--flux_resp_tv', metavar='time', dest='flux_resp_tv',
                      help='computes and plots time varying flux response coefficients (parameter: end_time)')
    parser.add_option('--conc_resp_tv', metavar='time', dest='conc_resp_tv',
                      help='computes and plots time varying concentration response coefficients (parameter: end_time)')
    parser.add_option('--spectral_conc_cc', metavar='frequency', dest='spectral_conc_cc',
                      help='compute spectral concentration control coefficients')
    parser.add_option('--spectral_flux_cc', metavar='frequency',
                      dest='spectral_flux_cc', help='compute spectral flux control coefficients')
    parser.add_option('--spectral_conc_resp', metavar='frequency',
                      dest='spectral_conc_resp', help='compute spectral conc response coefficients')
    parser.add_option('--spectral_flux_resp', metavar='frequency',
                      dest='spectral_flux_resp', help='compute spectral flux response coefficients')
    parser.add_option('--stoich_mat', dest='stoich_mat',
                      help='output the stoichiometric matrix')

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    doc = libsbml.readSBML(sys.argv[-1])
    model = doc.getModel()
    j = sbml_mca(model)

    if not options.normalize in [None, 'left', 'right', 'both']:
        raise Exception('Unknown normalization %s' % (options.normalize))

    print 'Reactions:'
    print '\n'.join([r.getId() for r in model.getListOfReactions()])
    print
    print 'Species:'
    print '\n'.join([s.getId() for s in model.getListOfSpecies()])
    print
    print 'Parameters:'
    print '\n'.join(j.get_parameter_ids())
    print
    print

    for param in ['timecourse', 'ss', 'p_ela', 'e_ela', 'p_ela2', 'flux_cc', 'conc_cc', 'flux_resp', 'conc_resp', 'flux_resp2', 'conc_resp2',
                  'flux_resp_tv', 'conc_resp_tv', 'spectral_conc_cc', 'spectral_flux_cc', 'spectral_conc_resp', 'spectral_flux_resp', 'stoich_mat']:
        if getattr(options, param) is not None:
            if param == 'timecourse':
                print 'Simulating the model'
                t = float(getattr(options, param))
                j.plot_timecourse(t, 100)
                print
            if param == 'ss':
                print 'Steady state concentrations:'
                print j.get_steady_state(string_output=options.string_output)
                print 'Steady state flux:'
                print j.get_steady_state_flux(string_output=options.string_output)
                print
            elif param == 'p_ela':
                print 'Parameter elasticities:'
                # if options.normalize!=None:
                #    print 'Warning: No normalization is done for elasticities'
                print j.get_parameter_elasticities(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'e_ela':
                print 'Substrate elasticities'
                # if options.normalize!=None:
                #    print 'Warning: No normalization is done for elasticities'
                print j.get_metabolite_elasticities(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'p_ela2':
                print 'Second order parameter elasticities'
                if options.normalize is not None:
                    print 'Warning: No normalization is done for elasticities'
                ela2 = j.get_2nd_elasticities()
                for pos, r in enumerate(model.getListOfReactions()):
                    print r.getId(), ':'
                    print ela2[pos]
                print
            elif param == 'flux_cc':
                print 'Flux control coefficients:'
                print j.get_flux_cc(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'conc_cc':
                print 'Concentration control coefficients:'
                print j.get_conc_cc(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'flux_resp':
                print 'Flux response coefficients:'
                print j.get_flux_resp(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'conc_resp':
                print 'Concentration response coefficients:'
                print j.get_conc_resp(string_output=options.string_output, normalize=options.normalize)
                print
            elif param == 'flux_resp2':
                print 'Second order flux response coefficients:'
                frc2 = j.get_2nd_flux_resp(normalize=options.normalize)
                for pos, r in enumerate(model.getListOfReactions()):
                    print r.getId(), ':'
                    print frc2[pos]
                print
            elif param == 'conc_resp2':
                print 'Second order concentration response coefficients:'
                crc2 = j.get_2nd_conc_resp(normalize=options.normalize)
                for pos, s in enumerate(j._species_ids):
                    print s, ':'
                    print crc2[pos]
                print
            elif param == 'flux_resp_tv':
                print 'Time varying flux response coefficients:'
                t = float(getattr(options, param))
                j.plot_time_varying_coefficients(
                    t, what='flux_resp', normalize=options.normalize)
                print
            elif param == 'conc_resp_tv':
                print 'Time varying concentration response coefficients:'
                t = float(getattr(options, param))
                j.plot_time_varying_coefficients(
                    t, what='conc_resp', normalize=options.normalize)
                print
            elif param == 'spectral_conc_cc':
                print 'Spectral concentration conctrol coefficients:'
                frequency = float(getattr(options, param))
                print j.get_spectral_conc_cc(frequency)
            elif param == 'spectral_flux_cc':
                print 'Spectral flux conctrol coefficients:'
                frequency = float(getattr(options, param))
                print j.get_spectral_flux_cc(frequency)
            elif param == 'spectral_conc_resp':
                print 'Spectral concentration response coefficients:'
                frequency = float(getattr(options, param))
                print j.get_spectral_conc_resp(frequency)
            elif param == 'spectral_flux_resp':
                print 'Spectral flux response coefficients:'
                frequency = float(getattr(options, param))
                print j.get_spectral_flux_resp(frequency)
            elif param == 'stoich_mat':
                print 'Stoichiometric matrix:'
                if options.string_output:
                    reacs = [r.getId() for r in j._model.getListOfReactions()]
                    print misc.matrix2string(j._N, reacs, j._species_ids)
                else:
                    print misc.matrix2string(j._N)
