#!/usr/bin/env python

import misc
import sys
import numpy
import libsbml
import scipy.optimize


class Thermodynamics:

    """ Class to compute perform thermodynamic checks on SBML models """

    def __init__(self, model):
        """
        @param model: libsbml.Model
        @type model: libsbml.Model
        """
        self._model = model
        self._N = self._get_stoich_mat(model)

    @staticmethod
    def _get_stoich_mat(model):
        """
         get stoichiometric matrix
        @param model: libsbml.model
        @type model: libsbml.model
        @return: stoich. matrix
        @rtype: numpy.array
        """
        species2pos = dict(zip([s.getId() for s in model.getListOfSpecies()], range(model.getNumSpecies())))
        N = numpy.zeros((model.getNumSpecies(), model.getNumReactions()))
        for i, r in enumerate(model.getListOfReactions()):
            modes = [(-1,'Reactants'), (1,'Products')]
            for sign, direction in modes:
                for sr in getattr(r, 'getListOf' + direction)():
                    s = model.getSpecies(sr.getSpecies())
                    j = species2pos[sr.getSpecies()]
                    N[j, i] = sign*sr.getStoichiometry()
        not_enzymes = misc.get_species_wo_enzymes(model)
        positions = [species2pos[sid] for sid in [s.getId() for s in not_enzymes]]
        return N[positions]

    def check_flux_signs(self, ss, fluxes, keq, enzymes=None, tol=1e-9):
        """
        check the sign of the fluxes
        @param ss: steady state values
        @type ss: numpy.array
        @param fluxes: flux values
        @type fluxes: numpy.array
        @param keq: equilibrium constants
        @type keq: numpy.array
        @param enzymes: list of enzyme concentrations (optional)
        @type enzymes: numpy.array
        @param tol: tolercance (default: 1e-9)
        @type tol: float
        @return: True (flux signs are OK) or Flase (not)
        @rtype: bool
        """
        if not(0. in ss or 0. in keq):
            ln_ss = numpy.log(ss)
            allowed_directions = numpy.dot(self._N.T, ln_ss) - numpy.log(keq) < tol    # sum n_i * ln(s_i) < ln(q) ?
        else:
            # the above is numerically not feasible if we have zero entries
            prod = numpy.ones(len(keq))
            for i, row in enumerate(self._N):
                if ss[i] == 0.:
                    continue
                prod *= ss[i] ** row
            allowed_directions = prod - keq < tol # mathematically equivalent to the above (hopefylly:)
        real_directions = fluxes >= tol
        diff = allowed_directions - real_directions

        if enzymes:
            zeroflux = flux < tol
            enzymezero = enzymes < tol
            bothzero = zeroflux * enzymezero
            diff = diff * 1 - bothzero # if an enzyme and a flux is zero, the eq const doesnt matter

        if sum(diff) > 0:
            sys.stderr.write("Flux for reaction %s has wrong sign.\n" %diff.tolist().index(True))
            return False
        return True

    def check_equilibrium_constants(self, keqs, tol=1e-14):
        """
        check if the equilibrium constants are thermodynamically feasible
        @param keqs: list of equilibrium constants
        @type keqs: numpy.array
        @param tol: tolerance (default: 1e-14)
        @type tol: float
        @return: True (OK) Flase (no)
        @rtype: bool
        """
        k = misc.nullspace (self._N)
        if k.shape[1] > 0:
            v = numpy.dot(k.T, numpy.log(keqs))
            if numpy.linalg.norm(v) > tol * len(v):
                return False
        return True

    def find_steady_state_flux(self, flux, ss, eq):
        """
        find a steady state flux which is consistent with the
        steady state concentrations and equilibrium constants
        @param flux: steady state flux
        @type flux: numpy.array
        @param ss: steady state concentrations
        @type ss: numpy.array
        @param eq: list of equilibrium constants
        @type eq: numnpy.array
        @return: new steady state flux
        @rtype: numpy.array
        """
        indices = [i for i, s in enumerate(misc.get_species_wo_enzymes(self._model))
                   if not(s.getConstant() or s.getBoundaryCondition())]
        f = lambda v: numpy.linalg.norm( numpy.dot(self._N[indices], v))
        nep = misc.get_not_enzyme_positions(self._model)
        ss_l = numpy.log(ss[nep])
        eq_l = numpy.log(eq)

        constraints=[]
        for i, row in enumerate(self._N.T):
            constraints.append(lambda x: -numpy.sign(numpy.dot(row, ss_l) - eq_l[i]) * numpy.sign(x[i]))
            
        new_flux = scipy.optimize.fmin_cobyla(f, flux, constraints, rhobeg=0.01, maxfun=1e7)
        return new_flux

    def check_steady_state(self, flux, tol=1e-3):
        """
        check whether the given flux vector satisfies the steady state condition ds/dt=0
        @param flux: steady state flux
        @type flux: numpy.array
        @param tol: tolerance (default: 1e-3)
        @type tol: float
        @return: None (raises Exception if steady state assumption is violated)
        @rtype: None
        """
        dsdt = numpy.dot(self._N, flux)
        flag=False
        for i,s in enumerate(misc.get_species_wo_enzymes(self._model)):
            if s.getBoundaryCondition() or s.getConstant():
                continue
            if abs(dsdt[i])>tol:
                sys.stderr.write('Steady state assumption for species %s violated\n' %s.getId())
                print s.getName()
                print dsdt[i]
                flag=True
        if flag:
            raise Exception('This is not a steady state')
        
if __name__ == '__main__':
    doc = libsbml.readSBML(sys.argv[1])
    model = doc.getModel()
    t = Thermodynamics(model)
    keq = numpy.ones(model.getNumReactions())
    keq[0]=2
    keq[1]=.501
    fluxes = numpy.ones(model.getNumReactions())
    ss = numpy.ones(model.getNumSpecies()) * 2
    t.check_flux_signs(ss, fluxes, keq)