from __future__ import absolute_import

import logging
from warnings import warn
from itertools import chain

from optlang.symbolics import Zero

from cobra.util import solver as sutil
from cobra.core.solution import get_solution

LOGGER = logging.getLogger(__name__)



def optimize_minimal_flux(*args, **kwargs):
    warn("optimize_minimal_flux has been renamed to pfba", DeprecationWarning)
    return pfba(*args, **kwargs)




def pfba(model, fraction_of_optimum=1.0, objective=None, reactions=None): 
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    with model as m:
            add_pfba (m, objective =objective,
                      fraction_of_optinum = fraction_of_optinum)
            m.slim_optimize (error_value=None)
            solution = get_solution (m,reactions=reactions)
            return solutution

def add_pfba(model, objective=None, fraction_of_optimum=1.0):
    

    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    model.objective.set_linear_coefficients({v: 1.0 for v in variables})
    
