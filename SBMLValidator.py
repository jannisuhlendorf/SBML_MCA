import libsbml

class SBMLValidator:
    """ class for validating an SBML model """

    _maximal_model_size = 4000

    def __init__(self, model):
        """
        @type model:  libsbml.model or string
        @param model: SBML model, or filename
        """
        if type(model)==str:
            self._doc = libsbml.readSBML(model)
            m = self._doc.getModel()
            self._model = m
        else:
            self._model=model.clone()
            self._doc = libsbml.SBMLDocument(model.getLevel(), model.getVersion())
            self._doc.setModel(self._model)

    def validate(self):
        """
        validate the model (this function should be called)
        @return: [warnings,errors] (each of them is a list of strings)
        """
        warnings = []
        errors = []
        errors += self._check_model_size()
        errors += self._check_not_supported_sbml_mca()
        warnings += self._check_initial_conditions()
        w,e = self._get_libsbml_errors()
        warnings += w
        errors += e
        return [warnings,errors]

    def _check_model_size(self):
        """ check the size of the model """
        errors = []
        no_sp = self._model.getNumSpecies()
        no_par = self._model.getNumParameters() # global parameters
        for r in self._model.getListOfReactions(): # local parameters
            no_par += r.getKineticLaw().getNumParameters()
        if no_sp*no_par > self._maximal_model_size:
            errors.append('Model is larger than maximal allowed size.')
        return errors

    def _check_not_supported_sbml_mca(self):
        """ check for features in the sbml that are not supported yet """
        errors = []
        if self._model.getNumConstraints():
            errors.append('Constraints not supported')
        if self._model.getNumEvents():
            errors.append('Events not supported')
        p_ids = [p.getId() for p in self._model.getListOfParameters()]
        for ia in self._model.getListOfInitialAssignments():
            if ia.getSymbol() not in p_ids:
                errors.append('Initial assignments are currently only supported for parameters.')
        for c in self._model.getListOfCompartments():
            if not c.getConstant():
                errors.append('Varying compartment sizes not supported')
        for rule in self._model.getListOfRules():
            if rule.isAlgebraic():
                errors.append('Algebraic rules not supported')

        parameter_ids = [p.getId() for p in self._model.getListOfParameters()]
        for rule in self._model.getListOfRules():
            if not rule.isRate():
                continue
            var = rule.getVariable()
            if var in parameter_ids:
                errors.append('Parameter %s is changed by a rate rule. This is not supported.' %(var))
        return errors

    def _check_initial_conditions(self):
        """ check if for each species initial conditions are given """
        warnings = []
        for s in self._model.getListOfSpecies():
            if s.getConstant():
                continue
            if not(s.isSetInitialConcentration() or s.isSetInitialAmount()):
                warnings.append('No initial value specified for species %s. Setting initial value to 0.' %s.getId())
        return warnings

    def _get_libsbml_errors(self):
        """ get errors from libsbml """
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, True)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, False)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, False)
        self._doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, False)
        # TODO: find out why this gives zero errors for cloned models
        errors = []
        warnings = []
        self._doc.checkConsistency()
        for i in range(self._doc.getNumErrors()):
            e = self._doc.getError(i)
            msg = e.getMessage()
            line = e.getLine()
            fatal = e.isFatal()
            e_type = 'warning'
            if e.isError():
                e_type = 'error'
            elif e.isInfo():
                e_type = 'info'
            error_msg = 'Libsbml %s in line %d: %s' %(e_type, line, msg)

            if fatal:
                errors.append(error_msg)
            else:
                warnings.append(error_msg)
        return [warnings,errors]
