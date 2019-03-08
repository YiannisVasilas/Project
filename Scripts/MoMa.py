from __future__ import absolute_import

from optlang.symbolics import Zero, add

import cobra.util.solver as sutil

from cobra.flux_analysis.parsimonious import pfba


def moma (model, solution=None ,linear = True):
	with model:
		add_moma (model=model ,solution =solution , linear = lenear)
		solution = model.oprimize()
		return solutiondef, optimize_minimal_flux(*args, **kwargs)
		warn("optimize_minimal_flux has been renamed to pfba", DeprecationWarning)
		return pfba(*args, **kwargs)
		
		"""
    Compute a single solution based on (linear) MOMA.

    Compute a new flux distribution that is at a minimal distance to a
    previous reference solution. Minimization of metabolic adjustment (MOMA) is
    generally used to assess the impact
    of knock-outs. Thus the typical usage is to provide a wildtype flux
    distribution as reference and a model in knock-out state.

    Parameters
    ----------
    model : cobra.Model
        The model state to compute a MOMA-based solution for.
    solution : cobra.Solution, optional
        A (wildtype) reference solution.
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).

    Returns
    -------
    cobra.Solution
        A flux distribution that is at a minimal distance compared to the
        reference solution"""
