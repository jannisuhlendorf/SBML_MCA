#!/usr/bin/env python
from Errors import NoncriticalError, CriticalError
import sys
import numpy
import libsbml
import warnings
import misc
import math
import sympy

sbo_type2id = {'activator': 21,
               'inhibitor': 20,
               'enzyme': 14}

TIME_VARIABLE = 'SBML_MCA_TIME_VARIABLE'
MAX_MODEL_SIZE = 2500


class Model:

    def __init__(self, sbml_model):
        """
        @type model:  libsbml.model or string
        @param model: SBML model, or filename
        """
        if type(sbml_model) == str:
            self._doc = libsbml.readSBML(sbml_model)
            self.sbml_model = self._doc.getModel()
        else:
            self.sbml_model = sbml_model.clone()
            self._doc = libsbml.SBMLDocument(
                sbml_model.getLevel(), sbml_model.getVersion())
            self._doc.setModel(self.sbml_model)

        self._N = None
        self._N_partitioned = None
        self._kinetic_laws = None
        self._external_species_concentrations = None
        self._rate_rules = None
        self._assignment_rules = None
        self._replacements = None
        self._species_2_position = None
        self._species_ids = None
        self._parameter_ids = None
        self._ode_variables = None
        self._enzyme_positions = None
        self._not_enzyme_positions = None
        self._species_volume_prefactor = None
        self._reaction_ids = None

        self._check_not_supported()
        misc.make_unique_local_parameters(self.sbml_model)

    @property
    def N(self):
        if self._N is None:
            """ get the stoichiometric matrix (not including constant and boundary condition species) """
            self._N = numpy.zeros((len(self.species_ids), self.sbml_model.getNumReactions()))
            species_ids_changed_by_rule = self.get_species(species_filter=self.is_species_changed_by_rule)
            for i, r in enumerate(self.sbml_model.getListOfReactions()):
                modes = [(-1, 'Reactants'), (1, 'Products')]
                for sign, direction in modes:
                    for sr in getattr(r, 'getListOf' + direction)():
                        s = self.sbml_model.getSpecies(sr.getSpecies())
                        if s.getBoundaryCondition() \
                                or s.getConstant() \
                                or s.getId() in species_ids_changed_by_rule:  # we dont want no constant species in our stoich matrix
                            continue
                        j = self.species_2_position[sr.getSpecies()]
                        self._N[j, i] += sign * sr.getStoichiometry()
            if len(self._N) == 0:
                raise CriticalError('Empty stoichiometric matrix.')
        return self._N

    @property
    def N_partitioned(self):
        """ get partitioned stoichiometric matrix """
        if self._N_partitioned is None:
            # compute reduced row echolon form to get the linear indep. rows
            rref, pivot = sympy.Matrix(self.N.T).rref()
            Nr = self.N[pivot]  # linear independent rows of N
            # link matrix is L = N*inv(Nr)  [because per definition N = L*Nr]
            L = numpy.dot(self.N, numpy.linalg.pinv(Nr))
            try:
                L_inv = numpy.linalg.inv(L)  # inverse of link matrix
            except:
                L_inv = None
            self._N_partitioned = [L_inv, L, Nr]
        return self._N_partitioned

    @property
    def kinetic_laws(self):
        if self._kinetic_laws is None:
            self._kinetic_laws = self._get_kinetic_laws(compile_formulas=True)
        return self._kinetic_laws

    @property
    def external_species_concentrations(self):
        if self._external_species_concentrations is None:
            self._external_species_concentrations = {}
            for s in self.get_species(species_filter = self.is_species_constant):
                if s.isSetInitialConcentration():
                    self._external_species_concentrations[s.getId()] = s.getInitialConcentration()
                elif s.isSetInitialAmount():
                    self._external_species_concentrations[s.getId()] = s.getInitialAmount()
        return self._external_species_concentrations

    @property
    def rate_rules(self):
        if self._rate_rules is None:
            self._rate_rules = {}
            for rule in self.sbml_model.getListOfRules():
                var = rule.getVariable()
                formula = misc.ast_to_string(rule.getMath(), self.assignment_rules, self.replacements, mode='')
                if rule.isRate():
                    self._rate_rules[var] = formula
                elif rule.isAlgebraic():
                    raise CriticalError('Algebraic rules not supported')
            self._rate_rules = self._replace_flux_symbols(self._rate_rules)
        return self._rate_rules

    @property
    def assignment_rules(self):
        if self._assignment_rules is None:
            is_loop = True
            self._assignment_rules = {}
            while is_loop:  # loop until no assignment rule is dependent on another assignment
                for rule in self.sbml_model.getListOfRules():
                    if rule.isAssignment():
                        var = rule.getVariable()
                        formula = misc.ast_to_string(rule.getMath(),
                                                     self.assignment_rules,
                                                     self.replacements,
                                                     mode='',
                                                     replace=True)
                        formula_wo_replace = misc.ast_to_string(rule.getMath(),
                                                                self.assignment_rules,
                                                                self.replacements,
                                                                mode='',
                                                                replace=False)
                        self._assignment_rules[var] = {True: formula, False: formula_wo_replace}
                # check dependencies
                is_loop = False
                for var1 in self._assignment_rules:
                    for var2 in self._assignment_rules:
                        if var2 in self._assignment_rules[var1][True]:
                            is_loop = True
        return self._assignment_rules

    @property
    def replacements(self):
        """ get dictionary of parameter values and compartmets"""
        if self._replacements is None:
            # do not take include parameters that are modified by a rule
            self._replacements = {}
            params_changed_by_rule = [r.getVariable() for r in self.sbml_model.getListOfRules()]
            for (pos, base) in enumerate([self.sbml_model] +
                                                 [r.getKineticLaw() for r in self.sbml_model.getListOfReactions()]):
                for p in base.getListOfParameters():
                    if not p.getId() in params_changed_by_rule:
                        self._replacements[p.getId()] = p.getValue()
            for comp in self.sbml_model.getListOfCompartments():
                s = 1.
                if comp.isSetSize():
                    s = comp.getSize()
                elif comp.isSetVolume():
                    s = comp.getVolume()
                self._replacements[comp.getId()] = s
            # handle initial assignments (they might be dependent on each other,
            # therefore try 50 evaluations)
            max_iter = 50
            # this has to be done manually, because get_initial_conc method depends
            # on parameters which are not yet known
            params_with_species = {}
            for s in self.sbml_model.getListOfSpecies():
                if s.isSetInitialConcentration():
                    params_with_species[s.getId()] = s.getInitialConcentration()
                elif s.isSetInitialAmount():
                    params_with_species[s.getId()] = s.getInitialAmount()

            while any([math.isnan(self._replacements[x]) for x in self._replacements]) and max_iter > 0:
                params_with_species.update(self._replacements)
                ass = self._evaluate_initial_assignments(params_with_species)
                for p in ass:
                    if self._replacements.has_key(p) and math.isnan(self._replacements[p]):
                        self._replacements[p] = float(ass[p])
                    max_iter -= 1
        return self._replacements

    @property
    def species_2_position(self):
        if self._species_2_position is None:
            self._species_2_position = dict(zip(self.species_ids, range(self.species_ids.__len__())))
        return self._species_2_position

    @property
    def species_ids(self):
        if self._species_ids is None:
            self._species_ids = [s.getId() for s in filter(self.is_species_not_constant,
                                                           self.sbml_model.getListOfSpecies())]
        return self._species_ids

    @property
    def parameter_ids(self):
        """ get list of parameters that are varied in p_elasticities """
        if self._parameter_ids is None:
            const_species_ids = [s.getId() for s in self.get_species(species_filter=self.is_species_constant)]
            params_changed_by_rule = [r.getVariable() for r in self.sbml_model.getListOfRules()]
            global_param_ids = []
            for p in self.sbml_model.getListOfParameters():
                if p.getId() not in params_changed_by_rule:
                    global_param_ids.append(p.getId())
            local_param_ids = []
            for r in self.sbml_model.getListOfReactions():
                kl = r.getKineticLaw()
                for p in kl.getListOfParameters():
                    if p.getConstant():
                        local_param_ids.append(p.getId())
            self._parameter_ids = const_species_ids + global_param_ids + local_param_ids
            for p_name in self._parameter_ids:
                if self._parameter_ids.count(p_name) != 1:
                    raise CriticalError(
                        'Parameter ID %s used multiple times. This is valid but not yet supported.' % p_name)
        return self._parameter_ids

    @property
    def ode_variables(self):
        """ get list of species IDs and parameter IDs which are modified by an ODE (ether take part in reaction
         or are changed by rate rule """
        if self._ode_variables is None:
            self._ode_variables = self.species_ids \
                                  + [p_id for p_id in self._rate_rules if not p_id in self.species_ids]
        return self._ode_variables

    @property
    def enzyme_positions(self):
        if self._enzyme_positions is None:
            self._enzyme_positions = []
            for pos, species in enumerate([self.sbml_model.getSpecies(id) for id in self.species_ids]):
                if species.getSBOTerm() == sbo_type2id['enzyme'] or species.getId().startswith('enzyme'):
                    self._enzyme_positions.append(pos)
        return self._enzyme_positions

    @property
    def not_enzyme_positions(self):
        if self._not_enzyme_positions is None:
            self._not_enzyme_positions = [i for i in range(len(self.species_ids)) if i not in self.enzyme_positions]
        return self._not_enzyme_positions

    def get_species(self, species_filter=None):
        if species_filter is None:
            return self.sbml_model.getListOfSpecies()
        else:
            return filter(species_filter, self.sbml_model.getListOfSpecies())

    @property
    def species_volume_prefactor(self):
        """ get a vector of conversion factors (particle number to concentration) for each species  """
        if self._species_volume_prefactor is None:
            comp_size = {comp.getId(): comp.getSize() for comp in self.sbml_model.getListOfCompartments()}
            for c in comp_size:
                if math.isnan(comp_size[c]):
                    comp_size[c] = 1.
            factors = []
            for s in [self.sbml_model.getSpecies(id) for id in self.species_ids]:
                if s.getHasOnlySubstanceUnits():
                    factors.append(1.)
                    continue
                factors.append(1. / comp_size[s.getCompartment()])
            self._species_volume_prefactor = numpy.array(factors)
        return self._species_volume_prefactor

    @property
    def reaction_ids(self):
        if self._reaction_ids is None:
            self._reaction_ids = [r.getId() for r in self.sbml_model.getListOfReactions()]
        return self._reaction_ids

    def get_parameter_values(self, parameter_ids=None):
        """
        Get values for the specified parameters
        @type   parameter_ids: list
        @param  parameter_ids: List of strings with parameter ids
        @rtype  numpy.array
        @return array with parameter values
        """
        if parameter_ids is None:
            parameter_ids = self.parameter_ids
        return numpy.array([misc.get_parameter_value(self.sbml_model, p) for p in parameter_ids])

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
            if p in self.external_species_concentrations:
                self.external_species_concentrations[p] = parameter_values[i]
            else:
                misc.set_parameter_value(self.sbml_model, p, parameter_values[i])
                rebuild = True

    def get_delta_parameters(self, d_param, parameter_names):
        """ enter a an array of parameter deviations and names and get back the corresponding d_param vector for all parameters
        @type    d_param:  numpy.array
        @param   d_param:  vector of parameter changes
        @type    parameter_names:  list of strings
        @param   parameter_names:  list of parameter names
        @rtype:            numpy.array
        @return:           vector of all parameter changes
        """
        all_p_names = self.parameter_ids
        dp = numpy.zeros(len(all_p_names))
        for value, name in zip(d_param, parameter_names):
            dp[all_p_names.index(name)] = value
        return dp

    def get_initial_conc(self, with_rate_rule_params=False):
        """
        get vector of initial concentrations
        @type    with_rate_rule_params: boolean
        @param   with_rate_rule_params: indicate whehter to include initial values for objects determined by rate rules
        @rtype:                         numpy.array
        @return:                        vector of initial concentrations
        """
        s0 = []
        for s in self.sbml_model.getListOfSpecies():
            if s.getConstant() or s.getBoundaryCondition() \
                    or s.getId() in self.assignment_rules.keys():
                continue
            if s.isSetInitialConcentration():
                s0.append(s.getInitialConcentration())
            elif s.isSetInitialAmount():
                s0.append(s.getInitialAmount())
            else:
                s0.append(0.)
                sys.stderr.write('No initial value specified for species %s. Setting initial value to 0.\n' % s.getId())
                # raise CriticalError('No initial value specified for species %s' %s.getId())
        if with_rate_rule_params:
            for var in self.rate_rules:
                if not var in self.species_ids:
                    s0.append(misc.get_parameter_value(self.sbml_model, var))
        return numpy.array(s0)

    def set_initial_conc(self, s0):
        """
        set the initial concentrations
        @type    s0:  numpy.array
        @param   s0:  vector of initial metabolite concentrations
        """
        pos = 0
        for s in self.sbml_model.getListOfSpecies():
            if s.getConstant() \
                    or s.getId() in self.assignment_rules.keys():
                continue
            s.setInitialConcentration(s0[pos])
            pos += 1

    def is_species_constant(self, s):
        return s.getConstant() or s.getBoundaryCondition() or self.is_species_changed_by_rule(s)

    def is_species_not_constant(self, s):
        return not self.is_species_constant(s)

    def is_species_changed_by_rule(self, s):
        rule_variables = [r.getVariable() for r in self.sbml_model.getListOfRules()]
        return s.getId() in rule_variables

    def get_model_id(self):
        return self.sbml_model.getId()

    def get_model_name(self):
        return self.sbml_model.getName()

    def get_species_ids(self, include_constant=False):
        """ get species ids of SBML model
            in contrast to self.model.species_ids this function also returns constant species
        """
        if include_constant:
            return [s.getId() for s in self.sbml_model.getListOfSpecies()]
        else:
            return self.species_ids

    def get_species_names(self, include_constant=False, replace_empty_names_with_id=False):
        """ get species names of SBML model
        include_constant specifies whether constant species are also returned
        replace_empty_names_with_id specifies whether empty species names are replaced by the species ID
        """
        if replace_empty_names_with_id:
            return [self.sbml_model.getSpecies(s).getName() or s
                    for s in self.get_species_ids(include_constant)]
        else:
            return [self.sbml_model.getSpecies(s).getName()
                    for s in self.get_species_ids(include_constant)]

    def get_species_initial_concentrations(self, include_constant=False):
        if include_constant:
            return [s.getInitialConcentration() for s in self.sbml_model.getListOfSpecies()]
        else:
            return [s.getInitialConcentration() for s in misc.get_not_constant_species(self.sbml_model)]

    def get_species_initial_amounts(self, include_constant=False):
        if include_constant:
            return [s.getInitialAmount() for s in self.sbml_model.getListOfSpecies()]
        else:
            return [s.getInitialAmount() for s in misc.get_not_constant_species(self.sbml_model)]

    def get_parameter_names(self, parameter_ids=None, replace_empty_names_with_id=False):
        if parameter_ids is None:
            parameter_ids = self.parameter_ids
        if replace_empty_names_with_id:
            return [misc.get_parameter_name(self.sbml_model, p_id) or p_id
                    for p_id in parameter_ids]
        else:
            return [misc.get_parameter_name(self.sbml_model, p_id)
                    for p_id in parameter_ids]

    def get_reaction_ids(self, replace_empty_names_with_id=False):
        return [r.getId() for r in self.sbml_model.getListOfReactions()]

    def get_reaction_names(self, replace_empty_names_with_id=False):
        if replace_empty_names_with_id:
            return [r.getName() or r.getId()
                    for r in self.sbml_model.getListOfReactions()]
        else:
            return [r.getName() for r in self.sbml_model.getListOfReactions()]

    def get_reaction_substrates(self):
        return [[s.getSpecies() for s in r.getListOfReactants()] for r in self.sbml_model.getListOfReactions()]

    def get_reaction_products(self):
        return [[s.getSpecies() for s in r.getListOfProducts()] for r in self.sbml_model.getListOfReactions()]

    def get_reaction_modifiers(self):
        return [[s.getSpecies() for s in r.getListOfProducts()] for r in self.sbml_model.getListOfReactions()]

    def get_reaction_kinetics(self):
        return [r.getKineticLaw().getFormula() for r in self.sbml_model.getListOfReactions()]

    def get_rule_ids(self):
        return [r.getId() for r in self.sbml_model.getListOfRules()]

    """ private functions """

    def _get_kinetic_laws(self, compile_formulas=True):
        """ gather string representation of the kinetic laws """
        formulas = []
        for pos, kl in enumerate([r.getKineticLaw() for r in self.sbml_model.getListOfReactions()]):
            formula = misc.ast_to_string(kl.getMath(),
                                         self.assignment_rules,
                                         self.replacements,
                                         mode='python')
            formula = formula.replace(' ', '')
            if compile_formulas:
                formula = compile(formula, 'kl_compile', 'eval')
            formulas.append(formula)
        return formulas

    def _replace_flux_symbols(self, formula_dict):
        """ replace references to reaction fluxes in rate rules """
        r_ids = [r.getId() for r in self.sbml_model.getListOfReactions()]
        kinetic_laws = self._get_kinetic_laws(compile_formulas=False)
        for id in formula_dict:
            for (r_id, kl) in zip(r_ids, kinetic_laws):
                formula_dict[id] = formula_dict[id].replace(r_id, kl)
        return formula_dict

    def _evaluate_initial_assignments(self, params):
        """ evaluate the intial assignment rules of the model with the given parameters """
        params = params.copy()
        assignments = {}
        for ass in self.sbml_model.getListOfInitialAssignments():
            var = ass.getSymbol()
            formula = misc.ast_to_string_libsbml(ass.getMath())
            # print formula
            assignments[var] = float(eval(formula, params))
        return assignments

    def _check_not_supported(self):
        """ check for features in the sbml that are not supported yet """
        if self.sbml_model.getNumConstraints():
            raise CriticalError('Error: Constraints not supported yet')
        if self.sbml_model.getNumEvents():
            raise CriticalError('Error: Events not supported yet')

        p_ids = [p.getId() for p in self.sbml_model.getListOfParameters()]
        for ia in self.sbml_model.getListOfInitialAssignments():
            if ia.getSymbol() not in p_ids:
                raise CriticalError(
                    'Initial assignments are currently only supported for parameters.')
        for c in self.sbml_model.getListOfCompartments():
            if not c.getConstant():
                raise CriticalError(
                    'Error: Varying compartment sizes not yet supported')
