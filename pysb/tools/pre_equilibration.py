import numpy as np
from pysb.simulator import ScipyOdeSimulator
from itertools import compress


def pre_equilibration(model, time_search, ligand, ligand_value=None, parameters=None, tolerance=1e-6):
    """

    :param model: PySB model
    :param ligand_idx: Species whose value want to be changed.
    :param time_search: time span array to be used to find the equilibrium
    :param tolerance: (tolerance, -tolerance) Range within equilibrium is considered as reached
    :param ligand_value: Initial condition of ligand (usually zero)
    :param parameters: Model parameters, must have same order as model.parameters
    :return:
    """
    if parameters is not None:
        # accept vector of parameter values as an argument
        if len(parameters) != len(model.parameters):
            raise Exception("param_values must be the same length as model.parameters")
        if not isinstance(parameters, np.ndarray):
            parameters = np.array(parameters)
    else:
        # create parameter vector from the values in the model
        parameters = np.array([p.value for p in model.parameters])

    param_dict = dict((p.name, parameters[i]) for i, p in enumerate(model.parameters))

    # Check if ligand name to be used for pre equilibration is provided
    if not isinstance(ligand, str):
        raise Exception('ligand must be a string with the parameter name')

    if ligand_value is not None:
        param_dict[ligand] = ligand_value
    else:
        param_dict[ligand] = 0

    # Solve system for the time span provided
    solver = ScipyOdeSimulator(model, tspan=time_search, param_values=param_dict).run()
    y = solver.species.T
    print (y)
    dt = time_search[1] - time_search[0]

    time_to_equilibration = [0, 0]
    for idx, sp in enumerate(y):
        sp_eq = False
        derivative = np.diff(sp) / dt
        derivative_range = ((derivative < tolerance) & (derivative > -tolerance))
        # Indexes of values less than tolerance and greater than -tolerance
        derivative_range_idxs = list(compress(range(len(derivative_range)), derivative_range))
        for i in derivative_range_idxs:
            # Check if derivative is close to zero in the time points ahead
            if (derivative[i + 3] < tolerance) | (derivative[i + 3] > -tolerance):
                sp_eq = True
                if time_search[i] > time_to_equilibration[0]:
                    time_to_equilibration[0] = time_search[i]
                    time_to_equilibration[1] = i
            if not sp_eq:
                raise Exception('Equilibrium can not be reached within the time_search input')
            if sp_eq:
                break
        else:
            raise Exception('Species s{0} has not reached equilibrium'.format(idx))

    conc_eq = y[:, time_to_equilibration[1]]
    return time_to_equilibration, conc_eq
