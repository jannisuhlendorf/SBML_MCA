import cStringIO, copy, libsbml, numpy, pylab
import sbml_mca

def tensor2string( ten, x=None, y=None, z=None, justify='center', delim=' | '):
    ret=''
    for pos,slice in enumerate(ten):
        if x:
            ret += '%s:\n' %x[pos]
        ret += matrix2string( slice, y, z, justify, delim) + '\n'
    return ret

def matrix2string( mat, head=None, left=None, justify='center', delim=' | '):
    max_length =  [max( (len(str(elem)) for elem in col) ) for col in mat.T]
    mat = mat.tolist()
    if head!=None:
        if len(mat[0])!=len(head):
            print len(mat[0]), len(head)
            raise Exception('Head dos not fix matrix')
        max_head = [len(str(x)) for x in head]
        max_length = [max(x) for x in zip(max_length,max_head)]
    if left!=None:
        if len(mat)!=len(left):
            print len(mat), len(left)
            raise Exception('Head dos not fix matrix')
        max_length = [max(len(str(x)) for x in left)]  + max_length
        head = [''] + copy.deepcopy(head)
        mat = [ [l]+rest for l,rest in zip(left,mat)  ]
    justify = {'left': str.ljust,'right': str.rjust, 'center': str.center}[justify.lower()]
    formatline = lambda line: delim.join( justify(x[0],x[1]) for x in zip( (str(y) for y in line), max_length ) )
    ret=''
    for line in [head or '']+mat:
        ret += formatline( map( str, line) ) + '\n'
    return ret

def matrix2python( mat ):
    ret='['
    for row in mat:
        ret+='['
        for e in row:
            ret+= str(e) + ', '
        ret+='], ' 
    ret+=']'
    return ret



def tensor2matlab( a, name ):
    ret='%s = zeros(%s,%s,%s)\n' %(name, a.shape[0], a.shape[1], a.shape[2])
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                ret += '%s(%s,%s,%s) = %s\n' %(name,i+1,j+1,k+1,a[i,j,k])
    return ret

def matrix2matlab( mat, name=None ):
    if not name:
        name='matrix'
    ret='%s = sparse( zeros(%s,%s) );\n' %(name,mat.shape[0],mat.shape[1])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i,j]==0:
                continue
            ret += '%s(%s,%s) = %s;\n' %(name, i+1, j+1, mat[i,j] )
    return ret

def make_unique_local_parameters( model ):
    """ Function to give all local parameters unique ids """
    #for pos,kl in enumerate([ r.getKineticLaw() for r in model.getListOfReactions() ]):
    for pos,r in enumerate( model.getListOfReactions() ):
        kl = r.getKineticLaw()
        for p in kl.getListOfParameters():
            new_p_id = '%s_local_%s' %(p.getId(),pos)
            f=' %s '%(kl.getFormula())
            for s in ['(',')','*','/','+','-',',']:
                f = f.replace(s,' %s '%s)
            old_id = p.getId()
            new_f = f.replace( ' '+old_id+' ', new_p_id )
            p.setId( new_p_id )
            reaction_base = r.getName() or r.getId()
            p.setName( reaction_base + '_' + old_id )
            kl.setFormula( new_f )
    return model

def get_constant_species( model ):
    def const(s): return ( s.getConstant() or s.getBoundaryCondition() )
    return filter(const, model.getListOfSpecies())
        
def get_not_constant_species( model ):
    def not_const(s): return not( s.getConstant() or s.getBoundaryCondition() )
    return filter(not_const, model.getListOfSpecies())


def is_enzyme( species ):
    return species.getSBOTerm()==14 or species.getId().startswith('enzyme')

def get_enzymes( model ):
    enzymes=[]
    for r in model.getListOfReactions(): 
        for s in [model.getSpecies(m.getSpecies()) for m in r.getListOfModifiers()]:
            if is_enzyme(s):
                enzymes.append(s)
                break
    if len(enzymes)!=model.getNumReactions():
        raise Exception('no enzyme found for reaction %s' %r.getId())
    return enzymes

def get_enzyme_positions(model):
    return [i for i,s in enumerate(model.getListOfSpecies()) if is_enzyme(s)]

def get_not_enzyme_positions(model):
    return [i for i,s in enumerate(model.getListOfSpecies()) if not is_enzyme(s)]

def get_species_wo_enzymes(model):
    enzymes = get_enzymes(model)
    return [ s for s in model.getListOfSpecies() if not s in enzymes]


def get_enzyme_for_reaction(model,reaction):
    is_enzyme = lambda s: s.getSBOTerm()==14 or s.getId().startswith('enzyme')
    for m in reaction.getListOfModifiers():
        s=model.getSpecies( m.getSpecies() )
        if is_enzyme(s):
            return s
    raise Exception('No enzyme was found for reaction %s' %reaction.getId() )


def get_enzyme_catalysed_reactions( model ):
    """ Find the reactions that are being catalysed by an enzyme (enzymatic reactions) """
    l=[]
    for r in model.getListOfReactions():
        try:
            get_enzyme_for_reaction(model,r)
            l.append(r)
        except:
            pass            
    return l

def get_positions_enzyme_catalysed_reactions( model ):
    ecr = get_enzyme_catalysed_reactions( model )
    lor = [ r for r in model.getListOfReactions() ]
    return [ lor.index(r) for r in ecr ]


def add_enzymes_to_reactions( model ):
    """ Create an enzyme for each reaction, if none exists (and modify kinetic law) """
    for r in model.getListOfReactions():
        try:
            get_enzyme_for_reaction(model,r)
        except:
            e = model.createSpecies()
            e.setId('enzyme_'+r.getId())
            e.setInitialConcentration(1)
            e.setSBOTerm(14)
            m=r.createModifier()
            m.setSpecies(e.getId())
            kl=r.getKineticLaw()
            kl.setFormula(e.getId()+' * '+kl.getFormula())
            
            
def get_parameter_value( model, _id ):
    for p in model.getListOfSpecies():
        if p.getId()==_id:
            return p.getInitialConcentration()
    for p in model.getListOfParameters():
        if p.getId()==_id:
            return p.getValue()
    for p in model.getListOfCompartments():
        if p.getId()==_id:
            return p.getVolume()
    for kl in [ r.getKineticLaw() for r in model.getListOfReactions() ]:
        for p in kl.getListOfParameters():
            if p.getId()==_id:
                return p.getValue()

def set_parameter_value( model, _id, value ):
    for p in model.getListOfSpecies():
        if p.getId()==_id:
            p.setInitialConcentration(value)
            return
    for p in model.getListOfParameters():
        if p.getId()==_id:
            p.setValue(value)
            return
    for p in model.getListOfCompartments():
        if p.getId()==_id:
            p.setVolume(value)
            return
    for kl in [ r.getKineticLaw() for r in model.getListOfReactions() ]:
        for p in kl.getListOfParameters():
            if p.getId()==_id:
                p.setValue(value)
                return
    raise Exception('Parameter %s not found' %_id)

def get_parameter_name( model, _id ):
    for p in model.getListOfParameters():
        if p.getId()==_id:
            return p.getName()
    for kl in [ r.getKineticLaw() for r in model.getListOfReactions() ]:
        for p in kl.getListOfParameters():
            if p.getId()==_id:
                return p.getName()
    return ''

def plot_ss_quantity( model, param, quantity, region=[0,10], steps=5 ):
    """ plot the steady state value of some quantity in dependence of some parameter """
    # find quantity
    sim = sbml_mca.sbml_mca(model)
    sp_pos = r_pos = None
    try:
        sp_pos = sim._species_ids.index(quantity)
    except:  pass
    try:
        r_pos  = [r.getId() for r in sim._model.getListOfReactions()].index(quantity)
    except: pass
    if sp_pos==None and r_pos==None:
        raise Exception('Quantity not found')
    
    p_values = numpy.linspace(region[0],region[1],steps)
    result=[]
    for p in p_values:
        set_parameter_value( model, param, p )
        sim = sbml_mca.sbml_mca(model)
        ss = sim.get_steady_state()
        if sp_pos!=None:
            x = ss[sp_pos]
        else:
            v = sim._v(ss,1)
            x = v[r_pos]
        result.append(x)
    pylab.plot( p_values, result )
    pylab.show()

        
def nullspace( matrix, tol=1e-14):
    """ Get the nullspace (kernel) of a matrix"""
    u,s,vh = numpy.linalg.svd(matrix)
    return vh[ s<tol ].T


def matrix2tensor( A, B ):
    # make a tensor from 2 matrices Tikl = Aik * Bil
    ret = numpy.zeros((A.shape[0], A.shape[1], B.shape[1] ))
    for i in range(ret.shape[0]):
        ret[i] = numpy.dot( numpy.matrix( A[i] ).T, numpy.matrix( B[i] ) )       
    return ret



def load_matrix_from_file( f_name ):
    f = open(f_name,'r')
    return numpy.array( [[float(x) for x in line.split()] for line in f.readlines()] )



def cov2cor( mat ):
    # convert covariance to correlation matrix
    mat = mat + numpy.diag( numpy.ones( mat.shape[0] )*10e-16 )
    d = numpy.diag( 1./numpy.sqrt(mat.diagonal()) )
    return numpy.dot( numpy.dot( d, mat ), d )


def print_reactions( model ):
    for i,r in enumerate(model.getListOfReactions()):
        d={'Reactants':[], 'Products':[]}
        for tp in d.keys():
            for sr in getattr( r, 'getListOf'+tp )():
                name = sr.getSpecies()
                if sr.getStoichiometry()!=1:
                    name = str(sr.getStoichiometry()) + ' ' + name
                d[tp].append(name)
        
        print i+1, '\t', r.getId()#, ': '
        #print ' + '.join(d['Reactants']) + ' -> ' + ' + '.join(d['Products'])


def plot_spectrum( model, quantity, param, region=[0,1], steps=100 ):
    sim = sbml_mca.sbml_mca(model)
    
    r_pos  = [r.getId() for r in sim._model.getListOfReactions()].index(quantity)
    p_pos = sim.get_parameter_ids().index(param)


    step=float(region[1]-region[0]) / (steps+1)

    x=[]
    result=[]
    for i in range(steps):
        w = region[0] + i*step
        x.append(w)
        result.append( abs( sim.get_spectral_flux_resp(w)[r_pos][p_pos] ) )

    pylab.plot( x, result )
    pylab.show()



def ast_code_to_string( ast_code ):
    ast_names = [ 'AST_CONSTANT_E', 'AST_CONSTANT_FALSE', 'AST_CONSTANT_PI', 'AST_CONSTANT_TRUE', 'AST_DIVIDE', 'AST_FUNCTION', 'AST_FUNCTION_ABS',
                  'AST_FUNCTION_ARCCOS', 'AST_FUNCTION_ARCCOSH', 'AST_FUNCTION_ARCCOT', 'AST_FUNCTION_ARCCOTH', 'AST_FUNCTION_ARCCSC', 'AST_FUNCTION_ARCCSCH',
                  'AST_FUNCTION_ARCSEC', 'AST_FUNCTION_ARCSECH', 'AST_FUNCTION_ARCSIN', 'AST_FUNCTION_ARCSINH', 'AST_FUNCTION_ARCTAN', 'AST_FUNCTION_ARCTANH',
                  'AST_FUNCTION_CEILING', 'AST_FUNCTION_COS', 'AST_FUNCTION_COSH', 'AST_FUNCTION_COT', 'AST_FUNCTION_COTH', 'AST_FUNCTION_CSC',
                  'AST_FUNCTION_CSCH', 'AST_FUNCTION_DELAY', 'AST_FUNCTION_EXP', 'AST_FUNCTION_FACTORIAL', 'AST_FUNCTION_FLOOR', 'AST_FUNCTION_LN',
                  'AST_FUNCTION_LOG', 'AST_FUNCTION_PIECEWISE', 'AST_FUNCTION_POWER', 'AST_FUNCTION_ROOT', 'AST_FUNCTION_SEC', 'AST_FUNCTION_SECH',
                  'AST_FUNCTION_SIN', 'AST_FUNCTION_SINH', 'AST_FUNCTION_TAN', 'AST_FUNCTION_TANH', 'AST_INTEGER', 'AST_LAMBDA', 'AST_LOGICAL_AND',
                  'AST_LOGICAL_NOT', 'AST_LOGICAL_OR', 'AST_LOGICAL_XOR', 'AST_MINUS', 'AST_NAME', 'AST_NAME_AVOGADRO', 'AST_NAME_TIME',
                  'AST_PLUS', 'AST_POWER', 'AST_RATIONAL', 'AST_REAL', 'AST_REAL_E', 'AST_RELATIONAL_EQ', 'AST_RELATIONAL_GEQ', 'AST_RELATIONAL_GT',
                  'AST_RELATIONAL_LEQ', 'AST_RELATIONAL_LT', 'AST_RELATIONAL_NEQ', 'AST_TIMES', 'AST_UNKNOWN' ]

    for name in ast_names:
        if ast_code == getattr( libsbml, name ): # they don't seem to be within a certain range
            return name
    raise Exception('AST code %s unknown' %ast_code)














if __name__=='__main__':
    import sys
    #print load_matrix_from_file(sys.argv[1])

    d = libsbml.readSBML(sys.argv[-1])
    m = d.getModel()
    plot_spectrum( m, 'reaction_1', 'A_GO_0005623' )












