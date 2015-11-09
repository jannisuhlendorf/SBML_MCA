import libsbml, sys, numpy
import misc, sbml_mca

def add_regulation( model , reaction, species, mode, k_mod, steady_state, test=False):
    #### TODO
    #### add ki/ka parameter to reaction
    if not mode in ['inh', 'act']:
        raise Exception('Unknown regulation type')
    
    m = model.clone()
    spec = m.getSpecies(species)
    reactions = [r.getId() for r in m.getListOfReactions()]
    
    #def const(s): return not(s.getBoundaryCondition() or s.getConstant())
    def const(s): return not(s.getConstant())
    species_ids             =    [ s.getId() for s in filter(const, m.getListOfSpecies()) ]
    try:
        species_ss_conc = steady_state[ species_ids.index(species) ]
    except:
        # constant species is regulator .. oculd make sense ....
        species_ss_conc = float(model.getSpecies(species).getInitialConcentration())
    
    # add a prefactor to the kinetic
    pos = reactions.index(reaction)

    if not k_mod:
        k_mod = species_ss_conc

    enzyme = misc.get_enzyme_for_reaction(m,m.getReaction(reaction))
    e_conc_old = float(enzyme.getInitialConcentration())

    if mode=='inh':
        k_mod_name='ki_'+reaction
        c_mod = '(%s/%s)' %(spec.getId(),k_mod_name)
        prefactor = '(1 / ( 1 + %s )) * ' %c_mod
        e_conc_new = (1 + (species_ss_conc/k_mod)) * e_conc_old
    else:
        if species_ss_conc<1e-5:
            raise Exception('Activator has a too small concentration')
        k_mod_name='ka_'+reaction
        c_mod = '(%s/%s)' %(spec.getId(),k_mod_name)
        prefactor = '(%s / (1 + %s)) * ' %(c_mod,c_mod)
        e_conc_new = (1+ (k_mod/species_ss_conc)) * e_conc_old
        if spec.getInitialConcentration()==0:
            sys.stderr.write('Initial concentration for essential activator 0. Setting to 0.001\n')
            spec.setInitialConcentration(0.001)
        
        
    kl = m.getReaction(reaction).getKineticLaw()
    param = libsbml.Parameter( k_mod_name, k_mod )
    kl.addParameter(param)
    kl.setFormula( prefactor + kl.getFormula() )

    ## adjust enzyme concentration
    enzyme.setInitialConcentration( e_conc_new )
    
    if test:
        # check wether we did everything right
        pos = species_ids.index( enzyme.getId() )
        ss2 = steady_state.copy()
        ss2[pos]=e_conc_new
        j1 = sbml_mca.sbml_mca(model)  # original model
        j2 = sbml_mca.sbml_mca(m)      # regulated model
        nep = j1._get_not_enzyme_positions()
        if numpy.linalg.norm(j1._v(steady_state,1.) - j2._v(ss2,1)) > 1e-6:
            raise Exception('Error: regulated steady state differst from original')
    return m

if __name__=='__main__':
    pass
    #d = libsbml
