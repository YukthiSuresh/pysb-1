from __future__ import absolute_import
from .base import SimulationResult, Simulator
from .scipyode import SerialExecutor
from pydas.dassl import DASSL
from concurrent.futures import ProcessPoolExecutor, Executor, Future
from functools import partial
import numpy as np
import sympy
import itertools
import pysb.bng
import re
try:
    import cython
except ImportError:
    cython = None


class DaeSimulator(Simulator):
    """
    Differential Algebraic Equation simulator using PyDAS and DASSL

    .. warning::
        The interface for this class is considered experimental and may
        change without warning as PySB is updated.

    This Simulator integrates PySB with the Fortran-based differential
    algebraic system solver [DASSL]_. It requires the [PyDAS]_ library to be
    installed - the easiest way to do this is with [Anaconda]_::

        conda install -c rmg pydas

    See the [PyDAS]_ README for more detailed installation instructions.

    .. [DASSL] https://www.osti.gov/scitech/servlets/purl/5882821 (PDF)
    .. [PyDAS] https://github.com/ReactionMechanismGenerator/PyDAS
    .. [Anaconda] https://www.anaconda.com/download/

    Parameters
    ----------
    model : pysb.Model
        Model to simulate.
    tspan : vector-like, optional
        Time values over which to simulate. The first and last values define
        the time range. Returned trajectories are sampled at every value unless
        the simulation is interrupted for some reason, e.g., due to
        satisfaction of a logical stopping criterion (see 'tout' below).
    initials : vector-like or dict, optional
        Values to use for the initial condition of all species. Ordering is
        determined by the order of model.species. If not specified, initial
        conditions will be taken from model.initial_conditions (with
        initial condition parameter values taken from `param_values` if
        specified).
    param_values : vector-like or dict, optional
        Values to use for every parameter in the model. Ordering is
        determined by the order of model.parameters.
        If passed as a dictionary, keys must be parameter names.
        If not specified, parameter values will be taken directly from
        model.parameters.
    verbose : bool or int, optional (default: False)
        Sets the verbosity level of the logger. See the logging levels and
        constants from Python's logging module for interpretation of integer
        values. False is equal to the PySB default level (currently WARNING),
        True is equal to DEBUG.
    **kwargs : dict
        Extra keyword arguments, including:

        * ``eqn_mode``: Equation evaluation mode - one of `cython`
        or `python`. They should give identical results, differing only in
        speed. `cython` does likewise with the `cython`
        library. `python` uses `sympy`'s `lambdify` function to execute the
        system of equations in pure Python, which is the most compatible but
        slowest option. The `python` option might be necessary if your model
        contains non-standard `sympy` functions for example. In general,
        we recommend `cython` where possible.

        * ``atol``: Absolute tolerance (default: 1e-16)

        * ``rtol``: Relative tolerance (default: 1e-8)

    Notes
    -----
    If ``tspan`` is not defined, it may be defined in the call to the
    ``run`` method.

    Examples
    --------
    Simulate a model and display the results for an observable:

    >>> from pysb.examples.robertson import model
    >>> import numpy as np
    >>> np.set_printoptions(precision=4)
    >>> sim = DaeSimulator(model, tspan=np.linspace(0, 40, 10), \
                           eqn_mode='python')
    >>> simulation_result = sim.run()
    >>> print(simulation_result.observables['A_total']) \
        #doctest: +NORMALIZE_WHITESPACE
    [ 1.      0.899   0.8506  0.8179  0.793   0.7728  0.7557  0.7408  0.7277
    0.7158]

    For further information on retrieving trajectories (species,
    observables, expressions over time) from the ``simulation_result``
    object returned by :func:`run`, see the examples under the
    :class:`SimulationResult` class.

    """
    _supports = {'multi_initials': True,
                 'multi_param_values': True}

    def __init__(self, model, tspan=None, initials=None,
                 param_values=None, verbose=False, compiler='cython', **kwargs):
        super(DaeSimulator, self).__init__(
            model=model, tspan=tspan, initials=initials,
            param_values=param_values, verbose=verbose, **kwargs)

        pysb.bng.generate_equations(self._model, verbose=self.verbose)

        self._eqn_subs = {e: e.expand_expr(expand_observables=True) for
                          e in self._model.expressions}
        n_species = len(self.model.species)
        dydt = sympy.symbols(','.join('__dydt{}'.format(i)
                                      for i in range(n_species)))
        dae_mat = sympy.Matrix([(self.model.odes[i] - dydt[i]) for i in
                                range(n_species)]).subs(self._eqn_subs)

        self._atol = kwargs.pop('atol', 1e-16)
        self._rtol = kwargs.pop('rtol', 1e-8)

        self._compiler = compiler
        if compiler == 'cython':
            if cython is None:
                raise ImportError('"cython" library not installed. Install '
                                  'cython or switch eqn_mode to "python"')
            code_eqs = '\n'.join(['delta[%d] = %s' %
                                  (i, sympy.ccode(o))
                                  for i, o in enumerate(dae_mat)])
            code_eqs = self._eqn_substitutions(code_eqs)
            code_eqs += '\nreturn delta'
            code_eqs = 'cdef double delta[{}]\n'.format(len(dae_mat)) + \
                       code_eqs

            # self._pdmodel = PyDasModelCython(code_eqs)
        elif compiler == 'python':
            self._symbols = sympy.symbols(','.join('__s%d' % sp_id for sp_id in
                                                   range(n_species))
                                          + ',') + tuple(model.parameters) + dydt

            code_eqs = (self._symbols, sympy.flatten(dae_mat))
            # self._pdmodel = PyDasModel(code_eqs_py)
        else:
            raise ValueError('Unknown compiler: {}'.format(compiler))
        self._code_eqs = code_eqs

        if kwargs:
            raise ValueError('Unknown keyword argument(s): {}'.format(
                ', '.join(kwargs.keys())
            ))

    def _eqn_substitutions(self, eqns):
        """String substitutions on the sympy C code for the ODE RHS and
        Jacobian functions to use appropriate terms for variables and
        parameters."""
        # Substitute 'y[i]' for '__si'
        eqns = re.sub(r'\b__s(\d+)\b',
                      lambda m: 'y[%s]' % (int(m.group(1))),
                      eqns)

        # Substitute 'dydt[i]' for '__dydti'
        eqns = re.sub(r'\b__dydt(\d+)\b',
                      lambda m: 'dydt[%s]' % (int(m.group(1))),
                      eqns)

        # Substitute 'p[i]' for any named parameters
        for i, p in enumerate(self._model.parameters):
            eqns = re.sub(r'\b(%s)\b' % p.name, 'p[%d]' % i, eqns)
        return eqns

    def run(self, tspan=None, initials=None, param_values=None,
            num_processors=1):
        """
        Run a simulation and returns the result (trajectories)

        .. note::
            In early versions of the Simulator class, ``tspan``, ``initials``
            and ``param_values`` supplied to this method persisted to future
            :func:`run` calls. This is no longer the case.

        Parameters
        ----------
        tspan
        initials
        param_values
            See parameter definitions in :class:`ScipyOdeSimulator`.
        num_processors : int
            Number of processes to use (default: 1). Set to a larger number
            (e.g. the number of CPU cores available) for parallel execution of
            simulations. This is only useful when simulating with more than one
            set of initial conditions and/or parameters.

        Returns
        -------
        A :class:`SimulationResult` object
        """
        super(DaeSimulator, self).run(tspan=tspan,
                                      initials=initials,
                                      param_values=param_values)

        if num_processors == 1:
            self._logger.debug('Single processor (serial) mode')
        else:
            self._logger.debug('Multi-processor (parallel) mode using {} '
                               'processes'.format(num_processors))

        num_species = len(self._model.species)

        with SerialExecutor() if num_processors == 1 else \
                ProcessPoolExecutor(max_workers=num_processors) as executor:
            sim_partial = partial(_run_dassl,
                                  code_eqs=self._code_eqs,
                                  n_species=num_species,
                                  tspan=self.tspan,
                                  compiler=self._compiler,
                                  atol=self._atol,
                                  rtol=self._rtol)

            results = [executor.submit(sim_partial, *args)
                       for args in zip(self.initials, self.param_values)]

            try:
                trajectories, tout_all = zip(*[r.result() for r in results])
            finally:
                for r in results:
                    r.cancel()

        self._logger.info('All simulation(s) complete')
        return SimulationResult(self, tout_all, trajectories)


def _run_dassl(initials, param_values, code_eqs, tspan, n_species,
               compiler, atol, rtol):
    if compiler == 'python':
        pdmodel = PyDasModel(sympy.lambdify(*code_eqs))
    else:
        pdmodel = PyDasModelCython(code_eqs)

    pdmodel.param_vec = param_values
    t0 = tspan[0]
    tout = [t0]
    y0 = initials
    trajectories = np.ndarray((len(tspan), n_species))
    trajectories[0] = y0
    dydt0 = -pdmodel.residual(
        t0, y0, np.zeros(len(y0), np.float64)
    )[0]
    pdmodel.initialize(
        t0, y0, dydt0, atol=atol, rtol=rtol
    )
    for t_ind in range(len(tspan) - 1):
        pdmodel.advance(tspan[t_ind + 1])
        tout.append(pdmodel.t)
        trajectories[t_ind + 1] = pdmodel.y

    return trajectories, tout


class PyDasModelCython(DASSL):
    def __init__(self, code_eqs):
        super(PyDasModelCython, self).__init__()
        self._code_eqs = code_eqs
        self.param_vec = None

    def residual(self, t, y, dydt):
        return np.array(cython.inline(self._code_eqs, p=self.param_vec)), 0


class PyDasModel(DASSL):
    def __init__(self, code_eqs):
        super(PyDasModel, self).__init__()
        self._code_eqs = code_eqs
        self.param_vec = None

    def residual(self, t, y, dydt):
        return np.array(self._code_eqs(*itertools.chain(
            y, self.param_vec, dydt))), 0
