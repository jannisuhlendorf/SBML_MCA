import Model
import Simulator
from Errors import NoncriticalError, CriticalError
import misc
import numpy
import sympy
import math
import pylab


class MCA:
    """ Class implementing various metabolic control analysis methods """
    
    _SMALL_FACTOR = 1e-12  # used to add to the trace of a singular jacobian
    MIN_VALUE_NORM = 1e-18  # minimal value that is divided by in normalization

    def __init__(self, model):
        """
        constructor for MCA class
        @type model: Model.Model
        @param model: model object
        """
        if not isinstance(model, Model.Model):
            raise CriticalError("MCA class requires Model object")
        self.model = model
        self.simulator = Simulator.Simulator(model)

        # fields to save already computed coeffiecients
        self._p_ela = None
        self._e_ela = None
        self._ccc = None
        self._fcc = None
        self._crc = None
        self._frc = None
        self._2nd_crc = None
        self._2nd_frc = None
        # save parameters for which the coefficients have been computed
        self._reference_params = {}

    def get_2nd_elasticities(self, parameter_names=None, ss_conc=None):
        """
        get the second order elasticities, if parameter_names is not speciefied for all parameters
        @type    parameter_names:  list
        @param   parameter_names:  list of parameter names
        @type    ss_conc:          numpy.array
        @param   ss_conc:          vector of steady state concentrations
        @rtype:                    numpy.array
        @return:                   3d-tensor with second order elasticities,
                                   1st index: reaction, 2nd index: param 1 3rd index: param 2
        """
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if parameter_names is None:
            parameter_names = self.model.parameter_ids
        value_dict = self._get_value_dict(ss_conc)
        ee = numpy.zeros((self.model.sbml_model.getNumReactions(), len(
            parameter_names), len(parameter_names)))
        for r_pos, kl in enumerate([r.getKineticLaw() for r in self.model.sbml_model.getListOfReactions()]):
            formula = misc.ast_to_string(kl.getMath(),
                                         self.model.assignment_rules,
                                         self.model.replacements,
                                         mode='python',
                                         replace=False)
            for p1_pos, p1_id in enumerate(parameter_names):
                for p2_pos, p2_id in enumerate(parameter_names[:p1_pos + 1]):
                    ela = self._get_2nd_elasticity(formula, value_dict, p1_id, p2_id)
                    ee[r_pos, p1_pos, p2_pos] = ee[r_pos, p2_pos, p1_pos] = ela
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
            ss_conc = self.simulator.get_steady_state()
        value_dict = self._get_value_dict(ss_conc)
        if time is not None:
            value_dict[Model.TIME_VARIABLE] = time
        ec = numpy.zeros((self.model.sbml_model.getNumReactions(), len(parameter_names)))

        for r_pos, kl in enumerate([r.getKineticLaw() for r in self.model.sbml_model.getListOfReactions()]):
            formula = misc.ast_to_string(kl.getMath(),
                                         self.model.assignment_rules,
                                         self.model.replacements,
                                         mode='python',
                                         replace=False)
            for p_pos, p_id in enumerate(parameter_names):
                if p_id not in formula:
                    continue
                e = self._get_elasticity(formula, value_dict, p_id)
                ec[r_pos, p_pos] = e
        if string_output:
            return misc.matrix2string(ec, parameter_names,
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
                                                       'right': (None, self.model.parameter_ids),
                                                       'both': ('flux', self.model.parameter_ids)}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._p_ela is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ep_ss'] - ss_conc) < 1e-8:
                    p_ela = self._p_ela
                    if normalize:
                        p_ela = normalize_wrap(p_ela)
                    if string_output:
                        return misc.matrix2string(p_ela, self.model.parameter_ids,
                                                  [r.getId() for r in self.model.sbml_model.getListOfReactions()])
                    return p_ela
            except:
                pass
        p_names = self.model.parameter_ids
        p_ela = self.get_elasticities(p_names, ss_conc)
        self._p_ela = p_ela
        if normalize:
            p_ela = normalize_wrap(p_ela, normalize)
        self._reference_params['ep_ss'] = ss_conc
        if string_output:
            return misc.matrix2string(p_ela, p_names,
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
                                                       'right': (None, 'conc'),
                                                       'both': ('flux', 'conc')}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._e_ela is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ee_ss'] - ss_conc) < 1e-8:
                    e_ela = self._e_ela
                    if normalize:
                        e_ela = normalize_wrap(e_ela)
                    if string_output:
                        return misc.matrix2string(e_ela, self.model.species_ids,
                                                  [r.getId() for r in self.model.sbml_model.getListOfReactions()])
                    return e_ela
            except:
                pass
        e_ela = self.get_elasticities(self.model.species_ids, ss_conc)
        self._reference_params['ee_ss'] = ss_conc
        self._e_ela = e_ela
        if normalize:
            e_ela = normalize_wrap(e_ela, normalize)
        if string_output:
            return misc.matrix2string(e_ela, self.model.species_ids,
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None, 'flux'),
                                                       'both': ('flux', 'flux')}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._fcc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['fcc_ss'] - ss_conc) < 1e-8:
                    fcc = self._fcc
                    if normalize:
                        fcc = normalize_wrap(fcc)
                    if string_output:
                        return misc.matrix2string(fcc, self.model.species_ids,
                                                  [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
            return misc.matrix2string(fcc, [r.getId() for r in self.model.sbml_model.getListOfReactions()],
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('conc', None),
                                                       'right': (None, 'flux'),
                                                       'both': ('conc', 'flux')}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._ccc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['ccc_ss'] - ss_conc) < 1e-8:
                    ccc = self._ccc
                    if normalize:
                        ccc = normalize_wrap(ccc, normalize)
                    if string_output:
                        return misc.matrix2string(ccc,
                                                  [r.getId() for r in self.model.sbml_model.getListOfReactions()],
                                                  self.model.species_ids)
                    return ccc.copy()
            except:
                pass
        e_ela = self.get_metabolite_elasticities(ss_conc)
        [L_inv, L, Nr] = self.model.N_partitioned
        M = numpy.dot(numpy.dot(Nr, e_ela), L)
        try:
            M_inv = numpy.linalg.inv(M)
        except:
            try:
                M_inv = numpy.linalg.inv(M + numpy.eye(len(M)) * self._SMALL_FACTOR)
            except:
                raise NoncriticalError('Error: Singular Jacobian')

        ccc = numpy.dot(-numpy.dot(L, M_inv), Nr)
        self._ccc = ccc.copy()
        self._reference_params['ccc_ss'] = ss_conc
        if normalize:
            ccc = normalize_wrap(ccc, normalize)
        if string_output:
            return misc.matrix2string(ccc,
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()],
                                      self.model.species_ids)
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('conc', None),
                                                       'right': (None, self.model.parameter_ids),
                                                       'both': ('conc', self.model.parameter_ids)}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._crc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['crc_ss'] - ss_conc) < 1e-8:
                    crc = self._crc
                    if normalize:
                        crc = normalize_wrap(crc, normalize)
                    if string_output:
                        return misc.matrix2string(crc, self.model.parameter_ids, self.model.species_ids)
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
            return misc.matrix2string(crc, self.model.parameter_ids, self.model.species_ids)
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                      {'left': ('flux', None),
                                                       'right': (None, self.model.parameter_ids),
                                                       'both': ('flux', self.model.parameter_ids)}[normalize])
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        if self._frc is not None:
            try:
                if numpy.linalg.norm(self._reference_params['frc_ss'] - ss_conc) < 1e-8:
                    frc = self._frc
                    if normalize:
                        frc = normalize_wrap(frc, normalize)

                    if string_output:
                        return misc.matrix2string(frc, self.model.parameter_ids,
                                                  [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
            return misc.matrix2string(frc,
                                      self.model.parameter_ids,
                                      [r.getId() for r in self.model.sbml_model.getListOfReactions()])
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        custom_ela = self.get_elasticities(parameter_names, ss_conc)
        ccc = self.get_conc_cc(ss_conc)
        crc = numpy.dot(ccc, custom_ela)
        if string_output:
            return misc.matrix2string(crc, parameter_names, self.model.species_ids)
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
        if self.model.rate_rules != {}:
            raise NoncriticalError('MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        custom_ela = self.get_elasticities(parameter_names, ss_conc)
        fcc = self.get_flux_cc(ss_conc)
        frc = numpy.dot(fcc, custom_ela)
        if string_output:
            return misc.matrix2string(frc, parameter_names, self.model.species_ids)
        if normalize:
            l, r = {'left': ('flux', None), 'right': (
                None, parameter_names), 'both': ('flux', parameter_names)}[normalize]
            frc = self._normalize_coefficients(frc, left=l, right=r)
        return frc

    def get_2nd_conc_resp(self, normalize=None):
        """
        get the 2nd order concentration response coefficients
        @type    normalize:        string or None
        @param   normalize:        one of (None, left, right, both) indicating how the coefficients should be normalized
        @rtype:                    numpy.array
        @return:                   3d tensor containing second order concentration response coefficients
        """
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if self._2nd_crc is None:
            self._compute_2nd_resp()
        if normalize:
            fst_order = self.get_conc_resp(normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('conc', None, fst_order),
                                                           'right': (None, self.model.parameter_ids, fst_order),
                                                           'both': ('conc', self.model.parameter_ids, fst_order)}[
                                                              normalize])
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not available for explicit ODE systems.')
        if self._2nd_frc is None:
            self._compute_2nd_resp()
        if normalize:
            fst_order = self.get_flux_resp(normalize=normalize)
            normalize_wrap = lambda mat, normalize: apply(self._normalize_coefficients, (mat,) +
                                                          {'left': ('flux', None, fst_order),
                                                           'right': (None, self.model.parameter_ids, fst_order),
                                                           'both': ('flux', self.model.parameter_ids, fst_order)}[
                                                              normalize])
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
        e_ela = self.get_metabolite_elasticities(ss_conc)
        L_inv, L, Nr = self.model.N_partitioned
        M = numpy.dot(numpy.dot(Nr, e_ela), L) - \
            (1j * frequency * numpy.eye(len(Nr)))
        try:
            M_inv = numpy.linalg.inv(M)
        except:
            try:
                M_inv = numpy.linalg.inv(
                    M + numpy.eye(len(M)) * self._SMALL_FACTOR)
            except:
                raise NoncriticalError('Error: Singular Jacobian')

        s_ccc = numpy.dot(-numpy.dot(L, M_inv), Nr)
        return s_ccc

    def get_spectral_flux_cc(self, frequency, ss_conc=None):
        """
        get spectral flux control coefficients
        @type    frequency:      number
        @param   frequency:      frequency of the perturbation
        @rtype:                  numpy.array
        @return:                 matrix of specetral flux response coefficients
        """
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
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
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        if ss_conc is None:
            ss_conc = self.simulator.get_steady_state()
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
        L_inv, L, Nr = self.model.N_partitioned
        return numpy.dot(numpy.dot(Nr, e_ela), L)

    def get_stability(self, ss_conc=None):
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

    def get_elasticities_symbolic(self, parameter_names):
        # TODO: maybe evaluate first parameters / external conc.
        if self.model.rate_rules != {}:
            raise NoncriticalError('MCA methods are not avaiable for explicit ODE systems.')
        ec = sympy.zeros(self.model.sbml_model.getNumReactions(), len(parameter_names))
        for r_pos, kl in enumerate([r.getKineticLaw()
                                    for r in self.model.sbml_model.getListOfReactions()]):
            formula = misc.ast_to_string(kl.getMath(),
                                         self.model.assignment_rules,
                                         self.model.replacements,
                                         mode='python',
                                         replace=False)
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
        @param end_time: time to compute time varying coeffs.
        @type end_time: float
        @param normalize: one of (None, left, right, both) indicating how the coefficients should be normalized
        @type normalize: string or None
        @param initial_cond_as_params: include the initial concentration as influencing parameters
        @type initial_cond_as_params: bool
        @param return_multiple: if True, return normalized, non-normalized coefficients and errors
        @type return_multiple: bool
        @param return_species_tc: return also time courses of species
        @type return_species_tc: bool
        @return: [time, time_varying_concentration_response_coefficients]
        @rtype: list
        """
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        p_id = self.model.parameter_ids
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self.model.species_ids)
        l_s = len(self.model.species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients

        if l_sp > Model.MAX_MODEL_SIZE:
            raise NoncriticalError(
                'Model is too large to compute time varying coefficients.')

        def my_ode(x, t):
            # generate ODE system to solve for time varying conc. response
            # coefficients
            # the first l_sp values are the time varying r.c., here they are
            # reshaped to matrix form
            dc = x[:l_sp].reshape((l_s, l_p))
            s = x[l_sp:]  # the species concentrations at current time point
            # species elasticities dv/dS
            es = self.get_elasticities(self.model.species_ids, s, time=t)
            # parameter elasticities dv/dp
            ep = self.get_elasticities(p_id, s, time=t)
            if initial_cond_as_params:
                # elasticities for initial conditions are zero
                ep = numpy.hstack(
                    (ep, numpy.zeros((self.model.sbml_model.getNumReactions(), l_s))))
            dx = numpy.zeros(l_sp + len(self.model.species_ids))
            # compuation: d/dt R = N( dv/ds * ds/dp + dv/dp )
            dx[:l_sp] = numpy.reshape(
                numpy.dot(self.model.N, (numpy.dot(es, dc) + ep)), -1)
            dx[l_sp:] = self.simulator.dS_dt(s, t)
            return dx

        s0 = self.model.get_initial_conc(with_rate_rule_params=True)
        c0 = numpy.zeros(l_sp)
        if initial_cond_as_params:
            # initialize derivatives w.r.t. initial conditions to 1.
            c0 = numpy.zeros((l_s, l_p))
            c0[:, len(p_id):] = numpy.eye(l_s)
            c0 = numpy.reshape(c0, -1)
        t = numpy.linspace(0, end_time, 100)
        result = self.simulator.odeint_wrapper(my_ode, numpy.concatenate((c0, s0)), t)
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
                    # raise NoncriticalError( 'There was an error during normalization.' )
                    raise NoncriticalError(
                        'Error: Normalization failed. Value to divide by is too small.')
            if normalize in ['right', 'both']:
                p_v = self.model.get_parameter_values(p_id)
                if initial_cond_as_params:
                    p_v = numpy.concatenate((p_v, self.model.get_initial_conc()))
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
        @param end_time: time to compute time varying coeffs.
        @type end_time: float
        @param normalize: one of (None, left, right, both) indicating how the coefficients should be normalized
        @type normalize: string or None
        @param initial_cond_as_params: include the initial concentration as influencing parameters
        @type initial_cond_as_params: bool
        @param return_flat: if True, return flattened (2d) version of coefficients
        @type return_flat: bool
        @param return_species_tc: return also time courses of species
        @type return_species_tc: bool
        @return: [time, time_varying_concentration_response_coefficients]
        @rtype: list
        """
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        p_id = self.model.parameter_ids
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self.model.species_ids)
        l_s = len(self.model.species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients

        if l_sp > Model.MAX_MODEL_SIZE:
            raise NoncriticalError(
                'Model is too large to compute time varying coefficients.')

        def my_ode(x, t):
            # generate ODE system to solve for time varying conc. response
            # coefficients
            # the first l_sp values are the time varying r.c., here they are
            # reshaped to matrix form
            dc = x[:l_sp].reshape((l_s, l_p))
            s = x[l_sp:]  # the species concentrations at current time point
            # species elasticities dv/dS
            es = self.get_elasticities(self.model.species_ids, s, time=t)
            # parameter elasticities dv/dp
            ep = self.get_elasticities(p_id, s, time=t)
            if initial_cond_as_params:
                # elasticities for initial conditions are zero
                ep = numpy.hstack((ep,
                                   numpy.zeros((self.model.sbml_model.getNumReactions(), l_s))))
            dx = numpy.zeros(l_sp + len(self.model.species_ids))
            # compuation: d/dt R = N( dv/ds * ds/dp + dv/dp )
            dx[:l_sp] = numpy.reshape(
                numpy.dot(self.model.N, (numpy.dot(es, dc) + ep)), -1)
            dx[l_sp:] = self.simulator.dS_dt(s, t)
            return dx

        s0 = self.model.get_initial_conc(with_rate_rule_params=True)
        c0 = numpy.zeros(l_sp)
        if initial_cond_as_params:
            # initialize derivatives w.r.t. initial conditions to 1.
            c0 = numpy.zeros((l_s, l_p))
            c0[:, len(p_id):] = numpy.eye(l_s)
            c0 = numpy.reshape(c0, -1)
        t = numpy.linspace(0, end_time, 100)
        result = self.simulator.odeint_wrapper(my_ode, numpy.concatenate((c0, s0)), t)
        tv_rc = result[:, :l_sp]

        if normalize is not None:
            left = right = None
            if normalize in ['left', 'both']:
                left = result[:, l_sp:]
            if normalize in ['right', 'both']:
                right = self.model.get_parameter_values(p_id)
                if initial_cond_as_params:
                    right = numpy.concatenate((right, self.model.get_initial_conc()))
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
        @param end_time: time to compute time varying coeffs.
        @type end_time: float
        @param normalize: one of (None, left, right, both) indicating how the coefficients should be normalized
        @type normalize: string or None
        @param initial_cond_as_params: include the initial concentration as influencing parameters
        @type initial_cond_as_params: bool
        @param return_flat: if True, return flattened (2d) version of coefficients
        @type return_flat: bool
        @param return_species_tc: return also time courses of species
        @type return_species_tc: bool
        @return: [timepoints, time_varying_flux_response_coefficients]
        @rtype: list
        """
        p_id = self.model.parameter_ids
        l_p = len(p_id)
        if initial_cond_as_params:  # if initial conditions are also condidered as parameters
            l_p += len(self.model.species_ids)
        l_s = len(self.model.species_ids)  # number of ODE species
        l_sp = l_p * l_s  # number of concentration response coefficients
        l_v = self.model.N.shape[1]  # nubmer of fluxes
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
            es = self.get_elasticities(self.model.species_ids, s, time=t)
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
                    left[i, :] = self.simulator.flux(species_tc[i, :], t)
            if normalize in ['right', 'both']:
                right = self.model.get_parameter_values(p_id)
                if initial_cond_as_params:
                    right = numpy.concatenate((right, self.model.get_initial_conc()))
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
        """
        get time varying response coefficients
        method for web-interface
        @param end_time: time to compute time varying coeffs.
        @type end_time: float
        @param what: conc or flux for concentration or flux response coefficients
        @type what: str
        @param initial_cond_as_params: include the initial concentration as influencing parameters
        @type initial_cond_as_params: bool
        @return:  [timepoints, tv_rc, tv_rc_norm, errors]
        @rtype: list
        """
        errors = []
        l_v = self.model.N.shape[1]  # nubmer of fluxes
        if what == 'conc':
            timepoints, tv_rc, species_tc = self.get_time_varying_conc_rc(end_time,
                                                                          normalize=None,
                                                                          initial_cond_as_params=initial_cond_as_params,
                                                                          return_flat=True,
                                                                          return_species_tc=True)
            left = species_tc
        elif what == 'flux':
            timepoints, tv_rc, species_tc = self.get_time_varying_flux_rc(end_time,
                                                                         normalize=None,
                                                                         initial_cond_as_params=initial_cond_as_params,
                                                                         return_flat=True,
                                                                         return_species_tc=True)
            left = numpy.zeros((len(timepoints), l_v))
            for i, t in enumerate(timepoints):
                left[i, :] = self.simulator.flux(species_tc[i, :], t)
        else:
            raise NoncriticalError('Unknown method %s' %what)
        right = self.model.get_parameter_values()
        if initial_cond_as_params:
            right = numpy.concatenate((right, self.model.get_initial_conc()))

        try:
            tv_rc_norm = self._normalize_time_varying_coeff(tv_rc, left, right)
        except Exception, e:
            errors.append(e)
            tv_rc_norm = numpy.zeros(tv_rc.shape)
        return [timepoints, tv_rc, tv_rc_norm, errors]

    def plot_time_varying_coefficients(self,
                                       end_time,
                                       what='conc_resp',
                                       initial_cond_as_params=False,
                                       normalize=None):
        """
        plot time varying coefficients using pylab
        @param end_time: time to compute time varying coeffs.
        @type end_time: float
        @param what: conc or flux for concentration or flux response coefficients
        @type what: str
        @param initial_cond_as_params: include the initial concentration as influencing parameters
        @type initial_cond_as_params: bool
        @param normalize: one of (None, left, right, both) indicating how the coefficients should be normalized
        @type normalize: string or none
        @return: None
        @rtype: None
        """
        if what == 'conc_resp':
            timepoints, tvrc = self.get_time_varying_conc_rc(end_time,
                                                             normalize=normalize,
                                                             initial_cond_as_params=initial_cond_as_params,
                                                             return_flat=True)
            p_ids = self.model.parameter_ids
            if initial_cond_as_params:  # if initial conditions are also condidered as parameters
                l_p += len(self.model.species_ids)
            r_ids = [r.getId() for r in self.model.sbml_model.getListOfReactions()]
            s_ids = self.model.species_ids
            coeff_names = numpy.reshape(
                [[s_id + '_' + p_id for s_id in s_ids] for p_id in p_ids], -1)
        elif what == 'flux_resp':
            timepoints, tvrc = self.get_time_varying_flux_rc(end_time,
                                                             normalize=normalize,
                                                             initial_cond_as_params=initial_cond_as_params,
                                                             return_flat=True)
            p_ids = self.model.parameter_ids
            if initial_cond_as_params:  # if initial conditions are also condidered as parameters
                l_p += len(self.model.species_ids)
            r_ids = [r.getId() for r in self.model.sbml_model.getListOfReactions()]
            coeff_names = numpy.reshape(
                [[r_id + '_' + p_id for r_id in r_ids] for p_id in p_ids], -1)
        pylab.plot(timepoints, tvrc)
        pylab.legend(coeff_names)
        pylab.show()

    def _normalize_coefficients(self, coeff, left=None, right=None, first_order=None):
        """
        normalize coefficients (2d or 3d)
        @param coeff: coefficients to be normalized
        @type coeff: numpy.array
        @param left: conc for steady state concentrations, flux for steady state fluxes
                     or numpy.ndarray for custom normalization
        @type left: str or numpy.ndarray
        @param right: conc for steady state concentrations, flux for steady state fluxes
                      or numpy.ndarray for custom normalization
        @type right: str or numpy.ndarray
        @param first_order: correction terms for normalization of second order terms
        @type first_order: numpy.ndarray
        @return: normalized coefficients
        @rtype: numpy.ndarray
        """
        if left is not None:
            if left == 'conc':
                d = self.simulator.get_steady_state()
            elif left == 'flux':
                d = self.simulator.flux(self.simulator.get_steady_state(), 1)
            elif isinstance(left, list):
                raise Exception('Impelment me...')
            elif isinstance(left, numpy.ndarray):
                d = left
            else:
                raise NoncriticalError('Unkown input for normalization')
            # for values very close to zero, we want inf/nan values for the
            # normalization result
            d[abs(d) < self.MIN_VALUE_NORM] = 0.
            with numpy.errstate(divide='print'):
                L = numpy.diag(1. / d)
            if coeff.shape.__len__() == 2:
                coeff = numpy.dot(L, coeff)
            if coeff.shape.__len__() == 3:
                coeff = numpy.tensordot(L, coeff, [1, 0])

        if right is not None:
            if isinstance(right, list):
                ss = self.simulator.get_steady_state()
                value_d = self._get_value_dict(ss)
                d = numpy.array([value_d[x] for x in right])
            elif isinstance(right, numpy.ndarray):
                d = right
            elif right == 'conc':
                d = self.simulator.get_steady_state()
            elif right == 'flux':
                d = self.simulator.flux(self.simulator.get_steady_state(), 1)
            else:
                raise NoncriticalError('Unkown input for normalization')
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
        """
        normalize time varying response coefficients
        @param coeff: time varying coefficient
        @type coeff: numpy.ndarray
        @param left: left side normalization (concentrations or fluxes)
        @type left: None or numpy.ndarray
        @param right: right side normalization
        @type right: None or numpy.ndarray
        @return: normalized coefficients
        @rtype: numpy.ndarray
        """
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
                    raise NoncriticalError(
                        'Error: Normalization failed. Value to divide by is too small.')
        if right is not None:
            no_right = len(right)
            if left is not None:
                assert no_right == no_coeff / left.shape[1]
            for i in range(no_coeff):
                coeff_copy[:, i] = coeff_copy[:, i] * right[i % no_right]
        return coeff_copy

    def _compute_2nd_resp(self, custom_params=None):
        """
        compute the second order responce coefficients
        @param custom_params: list of custom parameter to compute 2nd order response coefficients for
        @type custom_params: list
        @return: 2nd concentration response coefficients, 2nd flux response coefficients
        @rtype: list
        """
        if self.model.rate_rules != {}:
            raise NoncriticalError(
                'MCA methods are not avaiable for explicit ODE systems.')
        parameter_names = custom_params or self.model.parameter_ids
        p_set = self.model.species_ids + parameter_names
        ee = self.get_2nd_elasticities(p_set)
        cs = self.get_conc_cc()
        cj = self.get_flux_cc()
        rs = self.get_custom_conc_resp(parameter_names)
        ep = self.get_elasticities(parameter_names)
        es = self.get_metabolite_elasticities()
        size_s = len(self.model.species_ids)
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
        return rs2, rj2

    def _get_value_dict(self, species_conc):
        """
        get dictionary with standard parameter values
        @param species_conc: list of species concentrations
        @type species_conc: list
        @return: dictionary {species_id: concentration}
        @rtype: dict
        """
        value_dict = dict(zip(self.model.species_ids, species_conc.tolist()))
        value_dict.update(self.model.replacements)
        value_dict.update(self.model.external_species_concentrations)
        return value_dict

    @staticmethod
    def _get_elasticity(formula, value_dict, p_name):
        """
        compute elasticity (partial derivative) of formula w.r.t p_name
        @param formula: formula to derive
        @type formula: str
        @param value_dict: dictionary with values for parameters
        @type value_dict: dict
        @param p_name: parameter name to derive
        @type p_name: str
        @return: d(formula)/d(p_name)
        @rtype: float
        """
        value_dict = value_dict.copy()
        x = float(value_dict[p_name])
        if abs(x) > 1e-5:
            dx = abs(x) * 0.001
        else:
            # dx = 0.0005
            dx = 1e-9
        y0 = eval(formula, value_dict, {p_name: x - dx, 'math': math})
        y1 = eval(formula, value_dict, {p_name: x + dx, 'math': math})
        dy = y1 - y0
        return dy / (2 * dx)

    @staticmethod
    def _get_2nd_elasticity(formula, value_dict, p_name1, p_name2):
        """
        get the 2nd derivations of the input function with respect to p_name1 and p_name2
        @param formula: formula to derive
        @type formula: str
        @param value_dict: dictionary with values for parameters
        @type value_dict: dict
        @param p_name1: first parameter name to derive
        @type p_name1: str
        @param p_name2: second parameter name to derive
        @type p_name2: str
        @return:  d^2 f / dxdy
        @rtype: float
        """
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
        else:  # d^2 f / dxdy
            directions = [(-dx, -dy), (dx, -dy), (-dx, dy), (dx, dy)]
            for pos, [dx1, dx2] in enumerate(directions):
                a[pos] = eval(formula, value_dict, {p_name1: x + dx1, p_name2: y + dx2, 'math': math})
            der = (a[0] - a[1] - a[2] + a[3]) / (4 * dx * dy)
        return der

