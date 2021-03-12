# -*- coding: utf-8 -*-
# Taken from https://github.com/pescap/cald
from __future__ import division
import bempp.api
import numpy as np
import scipy.sparse.linalg
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator
from bempp.api.grid import grid_from_element_data
from scipy.sparse.linalg.interface import LinearOperator as _LinearOperator
from bempp.api.utils.logging import timeit as _timeit

def get_h(grid):
    #Routine that returns the Mesh Quality Information
    elements1 = list(grid.leaf_view.entity_iterator(1))
    vol = []
    for el1 in elements1:
        vol.append(el1.geometry.volume)
    vol = np.array(vol)
    return [vol.min(), vol.max(), vol.mean()]

class _it_counter(object):
    def __init__(self, store_residuals):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []

    def __call__(self, x):
        self._count += 1
        if self._store_residuals:
            self._residuals.append(np.linalg.norm(x))
            #print ('iteration -', self._count,"|| residual -", np.linalg.norm(x))
    @property
    def count(self):
        return self._count

    @property
    def residuals(self):
        return self._residuals

"""
Modification of linalg/iterative_solvers/gmres function: Allows to plug numpy matrices, LinearOperator or weak_forms() as an input
"""

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, use_strong_form=False, return_residuals=False):
    """Interface to the scipy.sparse.linalg.gmres function.

    This function behaves like the scipy.sparse.linalg.gmres function. But
    instead of a linear operator and a vector b it can take a boundary operator
    and a grid function. In this case, the result is returned as a grid function in the
    correct space.
    
    """
    import bempp.api
    import time

    if not isinstance(A, BoundaryOperator) and use_strong_form==True:
        raise ValueError("Strong form only with BoundaryOperator")
    if isinstance(A, BoundaryOperator) and not isinstance(b, GridFunction):
        raise ValueError("Instance Error")
    if not isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        raise ValueError("Instance Error")

    # Assemble weak form before the logging messages

    if isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        if use_strong_form:
            if not A.range.is_compatible(b.space):
                raise ValueError(
                    "The range of A and the space of A must have the same number of unknowns if the strong form is used.")
            A_op = A.strong_form()
            b_vec = b.coefficients
        else:
            A_op = A.weak_form()
            b_vec = b.projections(A.dual_to_range)
    else:
        A_op = A
        b_vec = b

    callback = _it_counter(return_residuals)

    
    start_time = time.time()
    x, info = scipy.sparse.linalg.gmres(A_op, b_vec,
                                        tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    
    if isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        res_fun = GridFunction(A.domain, coefficients=x.ravel())

    else:
        res_fun = x
    if return_residuals:
        return res_fun, info, callback.residuals
    else:
        return res_fun, info
