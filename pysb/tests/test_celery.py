from pysb.simulator.tasks import simulate
from pysb.examples import earm_1_0
import numpy as np
import tempfile
import os


def test_celery_simulate():
    init_kwargs = {'model': earm_1_0.model,
                   'tspan': np.linspace(0, 100, 101)}
    run_kwargs = [{'param_values': {'L_0': 10000}}]

    simulate('ScipyOdeSimulator', init_kwargs, run_kwargs)
    with tempfile.NamedTemporaryFile() as tf:
        if os.name == 'nt':
            tf.close()
        simulate('ScipyOdeSimulator', init_kwargs, run_kwargs,
                 save_location=tf.name, save_kwargs={'append': True})
