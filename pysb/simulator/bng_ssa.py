from pysb.simulator.base import Simulator, SimulationResult, SimulationResultNF
from pysb.bng import generate_equations, BngFileInterface
import numpy as np
import os


class BngSimulator(Simulator):
    _supports = {
        'multi_initials':     True,
        'multi_param_values': True
    }

    def __init__(self, model, tspan=None, cleanup=True, verbose=False):
        super(BngSimulator, self).__init__(model, tspan=tspan,
                                           verbose=verbose)
        self.cleanup = cleanup
        self._outdir = None
        generate_equations(self._model,
                           cleanup=self.cleanup,
                           verbose=self.verbose)

    def run(self, tspan=None, initials=None, param_values=None, n_sim=1,
            method='ssa', output_dir=None, output_file_basename=None,
            cleanup=True, verbose=False, **additional_args):
        """
        Simulate a model with BNG's SSA simulator and return the trajectories.

        Parameters
        ----------
        tspan: vector-like
            time span of simulation
        initials: vector-like, optional
            initial conditions of model
        param_values : vector-like or dictionary, optional
                Values to use for every parameter in the model. Ordering is
                determined by the order of model.parameters.
                If not specified, parameter values will be taken directly from
                model.parameters.
        n_sim: int, optional
            number of simulations to run
        method : str
            Type of simulation to run. Must be one of
            ['ssa', 'nf', 'pla', 'ode']
        output_dir : string, optional
            Location for temporary files generated by BNG. If None (the
            default), uses a temporary directory provided by the system. A
            temporary directory with a random name is created within the
            supplied location.
        output_file_basename : string, optional
            This argument is used as a prefix for the temporary BNG
            output directory, rather than the individual files.
        cleanup : bool, optional
            If True (default), delete the temporary files after the simulation is
            finished. If False, leave them in place. Useful for debugging.
        verbose: bool, optional
            If True, print BNG screen output.
        additional_args: kwargs, optional
            Additional arguments to pass to BioNetGen

        """
        super(BngSimulator, self).run(tspan=tspan,
                                      initials=initials,
                                      param_values=param_values,
                                      n_sim=n_sim,
                                      )

        if method not in ['ssa', 'nf', 'pla', 'ode']:
            print("Method must be one of ['ssa', 'nf', 'pla', 'ode']")
            quit()
        additional_args['method'] = method
        additional_args['t_end'] = np.max(self.tspan)
        additional_args['n_steps'] = len(self.tspan)
        additional_args['verbose'] = verbose
        params_names = [g.name for g in self._model.parameters]

        with BngFileInterface(self._model, verbose=verbose,
                              output_dir=output_dir,
                              output_prefix=output_file_basename,
                              cleanup=cleanup) as bngfile:
            bngfile.action('generate_network', overwrite=True, verbose=verbose)
            bngfile.action('saveConcentrations')
            if output_file_basename is None:
                prefix = 'pysb'
            else:
                prefix = output_file_basename
            for i in range(n_sim):

                tmp = additional_args.copy()
                tmp['prefix'] = prefix + str(i)
                for n, p_val in enumerate(self.param_values[i]):
                    bngfile.set_parameter('setParameter', params_names[n],
                                          self.param_values[i][n])
                for n, p_val in enumerate(self.initials[i]):
                    bngfile.set_concentration('setConcentration',
                                              self._model.species[n],
                                              self.initials[i][n])
                bngfile.action('simulate', **tmp)
            bngfile.execute()
            tout, yout = read_multi_simulation_results(n_sim,
                                                       bngfile.base_filename,
                                                       method)

        if method == 'nf':
            return SimulationResultNF(self, tout=tout, trajectories=yout)

        return SimulationResult(self, tout=tout, trajectories=yout)


def read_multi_simulation_results(n_sims, base_filename, method):
    """
    Reads the results of a BNG simulation and parses them into a numpy
    array
    """
    # Read concentrations data

    trajectories = [None] * n_sims
    tout = []
    # load the data

    for n in range(n_sims):
        filename = base_filename + str(n) + '.cdat'

        if method == 'nf':
            filename = base_filename + str(n) + '.gdat'
            # species = base_filename+str(n)+'.species'
            # with open(species, 'r') as f:
            #     out = f.read()

            if not os.path.isfile(filename):
                raise Exception("Cannot find input file " + filename)
            data = np.loadtxt(filename)
        else:
            data = np.loadtxt(filename, skiprows=1)

        # store data
        tout.append(data[:, 0])
        trajectories[n] = data[:, 1:]
    return np.array(tout), np.array(trajectories)
