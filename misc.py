from Errors import NoncriticalError, CriticalError
import Model
import copy
import libsbml
import numpy
import math

"""
helper functions
"""


def matrix2string(mat, head=None, left=None, justify='center', delim=' | '):
    """
    convert matrix to a string representation
    @param mat: matrix
    @type mat: numpy.ndarray
    @param head: list of column names
    @type head: list[str]
    @param left: list of row names
    @type left: list[str]
    @param justify: justify (left|center|right)
    @type justify: str
    @param delim: delimiter (default |)
    @type delim: str
    @return: string representation of matrix
    @rtype: str
    """
    max_length = [max((len(str(elem)) for elem in col)) for col in mat.T]
    mat = mat.tolist()
    if head is not None:
        if len(mat[0]) != len(head):
            raise Exception('Head dos not fix matrix')
        max_head = [len(str(x)) for x in head]
        max_length = [max(x) for x in zip(max_length,max_head)]
    if left is not None:
        if len(mat) != len(left):
            raise Exception('Head dos not fix matrix')
        max_length = [max(len(str(x)) for x in left)]  + max_length
        head = [''] + copy.deepcopy(head)
        mat = [[l] + rest for l, rest in zip(left,mat)]
    justify = {'left': str.ljust,'right': str.rjust, 'center': str.center}[justify.lower()]
    formatline = lambda line: delim.join(justify(x[0],x[1]) for x in zip((str(y) for y in line), max_length))
    ret=''
    for line in [head or ''] + mat:
        ret += formatline(map( str, line)) + '\n'
    return ret


def make_unique_local_parameters(model):
    """
    Function to give all local parameters unique ids
    @param model: libsbml.model
    @type model: libsbml.model
    @return: libsbml.model
    @rtype: libsbml.model
    """
    for pos,r in enumerate(model.getListOfReactions()):
        kl = r.getKineticLaw()
        for p in kl.getListOfParameters():
            new_p_id = '%s_local_%s' %(p.getId(),pos)
            f=' %s '%(kl.getFormula())
            for s in ['(', ')', '*', '/', '+', '-', ',']:
                f = f.replace(s,' %s '%s)
            old_id = p.getId()
            new_f = f.replace(' '+old_id+' ', new_p_id)
            p.setId(new_p_id)
            reaction_base = r.getName() or r.getId()
            p.setName(reaction_base + '_' + old_id)
            kl.setFormula(new_f)
    return model


def get_not_constant_species(model):
    """
    get species of the model that are not constant
    @param model: libsbml.model
    @type model: libsbml.model
    @return: list of species
    @rtype: list[libsbml.species]
    """
    def not_const(s): return not( s.getConstant() or s.getBoundaryCondition() )
    return filter(not_const, model.getListOfSpecies())


def is_enzyme(species):
    """
    is the given species an enzyme?
    @param species: libsbml.species
    @type species: libsbml.species
    @return: bool
    @rtype: bool
    """
    return species.getSBOTerm()==14 or species.getId().startswith('enzyme')


def get_enzymes(model):
    """
    get the enzymes of the model
    @param model: libsbml.model
    @type model: libsmbl.model
    @return: list of enzymes
    @rtype: list[libsbml.species]
    """
    enzymes=[]
    for r in model.getListOfReactions(): 
        for s in [model.getSpecies(m.getSpecies()) for m in r.getListOfModifiers()]:
            if is_enzyme(s):
                enzymes.append(s)
                break
    if len(enzymes)!=model.getNumReactions():
        raise Exception('no enzyme found for reaction %s' %r.getId())
    return enzymes


def get_not_enzyme_positions(model):
    """
    get positions of species that are not enzymes
    @param model: libsbml.model
    @type model: libsbml.model
    @return: list of positions
    @rtype: list[int]
    """
    return [i for i,s in enumerate(model.getListOfSpecies()) if not is_enzyme(s)]


def get_species_wo_enzymes(model):
    """
    gett species that are not enzymes
    @param model: libsbml.model
    @type model: libsbml.model
    @return: list of species
    @rtype: list[libsbml.species]
    """
    enzymes = get_enzymes(model)
    return [s for s in model.getListOfSpecies() if s not in enzymes]


def get_enzyme_for_reaction(model,reaction):
    """
    get the enzyme for a given reaction
    @param model: libsbml.model
    @type model: libsbml.model
    @param reaction: libsbml.reaction
    @type reaction: libsbml.reaction
    @return: enzme for the reaction
    @rtype: libsbml.species
    """
    is_enzyme = lambda s: s.getSBOTerm()==14 or s.getId().startswith('enzyme')
    for m in reaction.getListOfModifiers():
        s=model.getSpecies(m.getSpecies())
        if is_enzyme(s):
            return s
    raise Exception('No enzyme was found for reaction %s' %reaction.getId() )


def get_enzyme_catalysed_reactions(model):
    """
    Find the reactions that are being catalysed by an enzyme (enzymatic reactions)
    @param model: libsbml.model
    @type model: libsbml.model
    @return: list of reactions
    @rtype: list[libsbml.reaction]
    """
    l = []
    for r in model.getListOfReactions():
        try:
            get_enzyme_for_reaction(model, r)
            l.append(r)
        except:
            pass            
    return l


def add_enzymes_to_reactions(model):
    """
    Create an enzyme for each reaction, if none exists (and modify kinetic law)
    @param model: libsbml.model
    @type model: libsbml.model
    """
    for r in model.getListOfReactions():
        try:
            get_enzyme_for_reaction(model, r)
        except:
            e = model.createSpecies()
            e.setId('enzyme_'+r.getId())
            e.setInitialConcentration(1)
            e.setSBOTerm(14)
            m = r.createModifier()
            m.setSpecies(e.getId())
            kl = r.getKineticLaw()
            kl.setFormula(e.getId() + ' * ' + kl.getFormula())
            
            
def get_parameter_value(model, _id):
    """
    get a certain paraemter value
    @param model: libsbml.model
    @type model: libsbml.model
    @param _id: id of paraemter
    @type _id: str
    @return: parameter value
    @rtype: float
    """
    for p in model.getListOfSpecies():
        if p.getId() == _id:
            return p.getInitialConcentration()
    for p in model.getListOfParameters():
        if p.getId() == _id:
            return p.getValue()
    for p in model.getListOfCompartments():
        if p.getId() == _id:
            return p.getVolume()
    for kl in [r.getKineticLaw() for r in model.getListOfReactions()]:
        for p in kl.getListOfParameters():
            if p.getId() == _id:
                return p.getValue()


def set_parameter_value(model, _id, value):
    """
    set a parameter value
    @param model: libsbml.model
    @type model: libsbml.model
    @param _id: id of paraemter
    @type _id: str
    @param value: parameter value
    @type value: float
    """
    for p in model.getListOfSpecies():
        if p.getId() == _id:
            p.setInitialConcentration(value)
            return
    for p in model.getListOfParameters():
        if p.getId() == _id:
            p.setValue(value)
            return
    for p in model.getListOfCompartments():
        if p.getId() == _id:
            p.setVolume(value)
            return
    for kl in [ r.getKineticLaw() for r in model.getListOfReactions() ]:
        for p in kl.getListOfParameters():
            if p.getId() == _id:
                p.setValue(value)
                return
    raise Exception('Parameter %s not found' %_id)


def get_parameter_name(model, _id):
    """
    get the name of a parameter
    @param model: libsbml.model
    @type model: libsbml.model
    @param _id: id of paraemter
    @type _id: str
    @return: parmaeter name or empty string if no nameis given
    @rtype: str
    """
    for p in model.getListOfParameters():
        if p.getId()==_id:
            return p.getName()
    for kl in [ r.getKineticLaw() for r in model.getListOfReactions() ]:
        for p in kl.getListOfParameters():
            if p.getId()==_id:
                return p.getName()
    return ''


def nullspace(matrix, tol=1e-14):
    """
    Get the nullspace (kernel) of a matrix
    @param matrix: marix
    @type matrix: numpy.ndarray
    @param tol: tolerance (default 1e-14)
    @type tol: float
    @return: nullspace
    @rtype: numpy.ndarray
    """
    u,s,vh = numpy.linalg.svd(matrix)
    return vh[s<tol].T


def matrix2tensor(A, B):
    """
    make a tensor from 2 matrices Tikl = Aik * Bil
    @param A: matrix A
    @type A: numpy.ndarray
    @param B: matrix B
    @type B: numpy.ndarray
    @return: tensor
    @rtype: numpy.ndarray
    """
    ret = numpy.zeros((A.shape[0], A.shape[1], B.shape[1] ))
    for i in range(ret.shape[0]):
        ret[i] = numpy.dot(numpy.matrix(A[i]).T, numpy.matrix(B[i]))
    return ret


def ast_code_to_string(ast_code):
    """
    get string representation of libsbml AST code
    @param ast_code: ast code
    @type ast_code: int
    @return: string
    @rtype: str
    """
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


def get_formula(name, sbml_model, params, assignment_rules, replacements):
    """
    get a string representation for a function definition with parameters already replaced
    @param name: function name
    @type name: str
    @param sbml_model: libsbml model
    @type sbml_model: libsbml.model
    @param params: parameter values
    @type params: list
    @param assignment_rules: dictionary of assignment rules
    @type assignment_rules: dict
    @param replacements: dictionary of replacements
    @type replacements: dict
    @return: formula
    @rtype: str
    """
    func = sbml_model.getFunctionDefinition(name)
    p_dict = dict(zip([func.getArgument(x).getName() for x in range(func.getNumArguments())], params))
    # print p_dict
    # TODO: the formulaToString method does not work for functions such as log, exp, etc ...
    # TODO: unify the math conversion (decide whether to use the ast_to_string parser or the libsbml mehtod
    formula = ' ' + ast_to_string(func.getBody(), sbml_model, assignment_rules, replacements) + ' '
    for param in p_dict:
        formula = formula.replace(' ' + param + ' ', str(p_dict[param]))
    return formula


def ast_to_string(ast, sbml_model, assignment_rules, replacements, mode='python', replace=True):
    """
    convert libsbml AST node to string
    @param ast: ast node
    @type ast: libsbml.ast
    @param sbml_model: sbml model
    @type sbml_model: libsbml.model
    @param assignment_rules: dictionary of assignment rules
    @type assignment_rules: dict
    @param replacements: dictionary of replacements
    @type replacements: dict
    @param mode: conversion mode (python or fortran)
    @type mode: str
    @param replace: whether to replace parameters with their values
    @type replace: bool
    @return: string representation of the formula
    @rtype: str
    """
    if ast is None:
        return
    type = ast.getType()
    l = ast_to_string(ast.getLeftChild(), sbml_model, assignment_rules, replacements, mode, replace)
    r = ast_to_string(ast.getRightChild(), sbml_model, assignment_rules, replacements, mode, replace)
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
            raise CriticalError(
                'AST_FUNCTION_PIECEWISE not yet implemented completely')
        condition = ast_to_string(ast.getChild(1), sbml_model, assignment_rules, replacements, mode, replace)
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
            c = ast_to_string(ast.getChild(pos), sbml_model, assignment_rules, replacements, mode, replace)
            s += '  and  ' + c
        s += ')'
        return s
    elif type == libsbml.AST_LOGICAL_OR:
        s = '( %s' % l
        for pos in range(1, ast.getNumChildren()):
            c = ast_to_string(ast.getChild(pos), sbml_model, assignment_rules, replacements, mode, replace)
            s += '  or  ' + c
        s += ')'
        return s
    elif type == libsbml.AST_LOGICAL_XOR:
        s = '(bool( %s )' % l
        for pos in range(1, ast.getNumChildren()):
            c = ast_to_string(ast.getChild(pos), sbml_model, assignment_rules, replacements, mode, replace)
            s += ' ^ bool( %s )' % c
        s += ')'
        return s
    elif type == libsbml.AST_NAME:
        name = ast.getName()
        if name in assignment_rules:  # assignment rules are always replaced
            return str(assignment_rules[name][replace])
        if replace:
            try:
                return str(replacements[name])
            except:
                pass
        return name
    elif type == libsbml.AST_FUNCTION:
        children = [ast_to_string(ast.getChild(x), sbml_model, assignment_rules, replacements, mode, replace)
                    for x in range(ast.getNumChildren())]
        return get_formula(ast.getName(), sbml_model, children, assignment_rules, replacements)
    elif type == libsbml.AST_NAME_TIME:
        return Model.TIME_VARIABLE
    elif type == libsbml.AST_POWER:
        raise CriticalError('AST_POWER not yet implemented')
    else:
        ast_name = ast_code_to_string(type)
        print ast_name
        raise CriticalError('mathematic fucntion %s not yet implemented' % ast_name)


def ast_to_string_libsbml(ast):
    """
    convert libsbml AST node to string, using libsbml
    @param ast: libsbml.ast
    @type ast: libsbml.ast
    @return: string
    @rtype: str
    """
    formula = libsbml.formulaToString(ast)
    for old, new in [('log', 'math.log')]:
        formula = formula.replace(old, new)
    return formula
