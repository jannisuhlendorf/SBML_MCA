from Errors import NoncriticalError, CriticalError
import Model
import misc
import numpy
import scipy
import scipy.integrate
import scipy.optimize
import time
import pylab
import math


class Simulator:
    """ Class for the simulation of SBML models """

    def __init__(self, model):
        """
        @param model: SBML_MCA model object
        @type model: Model.Model
        """
        if not isinstance(model, Model.Model):
            raise CriticalError("Simulator requires Model object")
        self.model = model

        # save computed steady state
        self._steady_state = None
        self._no_steady_state = False
        self._steady_state_s0 = None

    def get_steady_state(self, s0=None, tol=1e-12, max_steps=100, max_computation_time=180, old_ss=None,
                         d_params=None, old_conc_resp=None, string_output=False):
        """
        find a steady state of the system, if the steady state before a parameter perturbation (old_ss), the
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
            raise NoncriticalError('No steady state could be found')
        if s0 is None:
            s0 = self.model.get_initial_conc(with_rate_rule_params=True)
        self._steady_state_s0 = s0
        if self._steady_state is not None and numpy.linalg.norm(self._steady_state_s0 - s0) < 0.0001:
            if string_output:
                return misc.matrix2string(numpy.array([self._steady_state]), self.model.species_ids)
            return self._steady_state
        # linear approximation of new ss
        if old_ss is not None and old_conc_resp is not None and d_params is not None:
            s0 = old_ss + numpy.dot(old_conc_resp, d_params)
        result = [s0]

        start_comp_time = time.time()
        sim_time = 0
        for i in range(max_steps):
            stepsize = 10 ** (i / 10 - 1)
            t = numpy.linspace(sim_time, sim_time + stepsize, 100)
            s0 = result[-1]
            try:
                result = self.odeint_wrapper(self.dS_dt, s0, t)
            except NoncriticalError, e:
                if i > 0:
                    self._no_steady_state = True
                    raise NoncriticalError(
                        'No steady state could be found. This could indicate and oszillating model.')
                else:
                    raise e
            if (numpy.linalg.norm(result[0] - result[-1]) / s0.size) < tol:
                break
            sim_time += stepsize
            if i == (max_steps - 1) or (time.time() - start_comp_time > max_computation_time):
                self._no_steady_state = True
                raise NoncriticalError(
                    'No steady state could be found. This could indicate and oszillating model.')
        self._steady_state = result[-1]
        if string_output:
            return misc.matrix2string(numpy.array([result[1]]).T,
                                      ['Conc.'],
                                      self.model.species_ids,
                                      justify='left')
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
        flux = self.flux(ss_conc, 1)
        if string_output:
            return misc.matrix2string(numpy.array([flux]).T,
                                      ['Flux'],
                                      self.model.reaction_ids,
                                      justify='left')
        return flux

    def solve_steady_state(self):
        """
        solve steady state by a numerical optimization
        @rtype:                 numpy.array
        @return:                vector of steady state concentrations
        """
        # optimization has to be done on the reduced system
        # TODO: implement different comp. sizes
        s0 = self.model.get_initial_conc()
        [L_inv, L, _] = self.model.N_partitioned
        si = numpy.dot(L_inv, s0)
        t = s0 - numpy.dot(L, si)
        f = lambda x: numpy.linalg.norm(
            self.dS_dt(numpy.dot(L, x) + t, 1))
        ss_i = scipy.optimize.fmin_bfgs(f, si)
        ss = numpy.dot(L, ss_i) + t
        return ss

    def dS_dt(self, species_conc, t):
        """
        compute rate of metabolite change
        @param species_conc: list of species concentrations
        @type species_conc: list[float]
        @param t: current time
        @type t: float
        @return: array of rate of changes for each substrate
        @rtype: numpy.array
        """
        ret = self.model.species_volume_prefactor * numpy.dot(self.model.N, self.flux(species_conc, t))
        # handle rate rules
        for var in self.model.rate_rules:
            f = self.model.rate_rules[var].replace(Model.TIME_VARIABLE, str(t))
            f = self.model.rate_rules[var].replace(Model.TIME_VARIABLE, str(t))
            species2conc = dict(zip(self.model.ode_variables, species_conc))
            species2conc['math'] = globals()['math']
            # rate = eval( f, species2conc, self._external_species_conc )
            rate = eval(f, self.model.external_species_concentrations, species2conc)
            if self.model.species_2_position.has_key(var):
                ret[self.model.species_2_position[var]] = rate
            else:
                l = ret.tolist()
                l.append(rate)
                ret = numpy.array(l)
        return ret

    def flux(self, species_conc, t):
        """
        get reaction velocities at given point
        @param species_conc: list of species concentrations
        @type species_conc: list[float]
        @param t: current time
        @type t: float
        @return: array of reaction rate velocities
        @rtype: numpy.array
        """
        vec = []
        species2conc = dict(zip(self.model.ode_variables, species_conc))
        species2conc[Model.TIME_VARIABLE] = t
        species2conc['math'] = globals()['math']

        for func in self.model.kinetic_laws:
            vec.append(eval(func, self.model.external_species_concentrations, species2conc))
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
            s0 = self.model.get_initial_conc(with_rate_rule_params=True)
        t = numpy.linspace(0, time, steps)
        return [t, self.odeint_wrapper(self.dS_dt, s0, t, rtol=r_tol, atol=a_tol)]

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
        for s in self.model.get_species(species_filter=self.model.is_species_constant):
            if s.isSetInitialConcentration():
                conc = s.getInitialConcentration()
            elif s.isSetInitialAmount():
                conc = s.getInitialAmount()
            elif self.model.assignment_rules.has_key(s.getId()):
                continue
            else:
                raise NoncriticalError(
                    'no initial value given for species %s' % (s.getId()))
            siz = trace.shape
            trace = numpy.insert(trace, siz[1], conc, axis=1)
        return t, trace

    def integrate_return_dict(self, time, steps, s0=None, r_tol=1e-6, a_tol=1e-12, with_constant_species=True,
                              with_assignment_rules=True):
        """
        integrates the model and returns a dictionary as result
        @param time: current time
        @type time: float
        @param steps: number of steps
        @type steps: int
        @param s0: steady state concentration (optional)
        @type s0: numpy.array
        @param r_tol: relative solver tolerance
        @type r_tol: float
        @param a_tol: absolute solver tolernace
        @type a_tol: float
        @param with_constant_species: include constant species
        @type with_constant_species: bool
        @param with_assignment_rules: include assignment rules
        @type with_assignment_rules: bool
        @return: dictionary with simulation results {species: result}
        @rtype: dict
        """
        s_ids = self.model.species_ids + \
                [p for p in self.model.rate_rules if not p in self.model.species_ids]
        if with_constant_species:
            t, trace = self.integrate_with_constant_species(
                time, steps, s0, r_tol, a_tol)
            s_ids = s_ids + [s.getId() for s in self.model.get_species(species_filter=self.model.is_species_constant)]

        else:
            t, trace = self.integrate(time, steps, s0, r_tol, a_tol)

        result_d = {s: trace[:, pos] for pos, s in enumerate(s_ids)}
        result_d['time'] = t

        if with_assignment_rules:
            species_assigned = self.compute_assignment_rules(trace, s_ids)
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
        for p in self.model.not_enzyme_positions + range(len(self.model.species_ids), len(result)):
            pylab.plot(t, result[p])
        pylab.xlabel('time')
        pylab.ylabel('concentration')
        if legend:
            pylab.legend([self.model.species_ids[i] for i in self.model.not_enzyme_positions]
                         + [p for p in self.model.rate_rules if not p in self.model.species_ids])
        pylab.show()

    def compute_assignment_rules(self, simulation_result, species_ids):
        """
        compute assignment rules for species after a simulation
        @param simulation_result: matrix of simulation result
        @type simulation_result: numpy.ndarray
        @param species_ids: list of species IDs
        @type species_ids: list[str]
        @return: dictionary of assignment rule results {assignment_rule: result}
        @rtype: dict
        """
        simulation_result_dict = {}
        for pos, id in enumerate(species_ids):
            simulation_result_dict[id] = simulation_result[:, pos]
        result = {}
        for ar in self.model.assignment_rules:
            result[ar] = eval(self.model.assignment_rules[ar][True],
                              self.model.external_species_concentrations,
                              simulation_result_dict)
        return result

    @staticmethod
    def odeint_wrapper(d_func, x0, timepoints, rtol=1e-6, atol=1e-12):
        """
        wrapper function for scipy.integrate.odepack.odeint adding error handling
        @param d_func: derivative function
        @type d_func: function
        @param x0: list of intial values
        @type x0: numpy.array
        @param timepoints: list of time points
        @type timepoints: numpy.array
        @type r_tol:  number
        @param r_tol:  relative tolervance
        @type a_tol:  number
        @param a_tol:  absolute tolerance
        @return: result matrix
        @rtype: numpy.ndarray
        """
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
            raise NoncriticalError(
                'There was an error integrating the ODEs (' + errors[0].message + ')')
        if infodict['message'] != 'Integration successful.':
            raise NoncriticalError('There were numerial problems simulating the model.')
        return result
