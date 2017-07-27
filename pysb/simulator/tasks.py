import celery
try:
    import raven
    import raven.contrib.celery
except ImportError:
    raven = None
import os
import pysb.simulator

RESULT_TO_BACKEND = '__RESULT_TO_BACKEND__'


class Celery(celery.Celery):
    def on_configure(self):
        if 'PYSB_SENTRY_DSN' in os.environ:
            if raven is None:
                raise ImportError('This option requires the "raven" library')

            client = raven.Client(os.environ['PYSB_SENTRY_DSN'])

            # Register a custom filter to filter out duplicate logs
            raven.contrib.celery.register_logger_signal(client)

            # Hook into the Celery error handler
            raven.contrib.celery.register_signal(client)


class CeleryConfig(object):
    task_serializer = 'pickle'
    result_serializer = 'pickle'
    accept_content = ['pickle']
    broker_url = os.environ.get('PYSB_CELERY_URL', None)

app = Celery()
app.config_from_object(CeleryConfig)
if 'PYSB_CELERY_CONFIG' is os.environ:
    app.config_from_envvar('PYSB_CELERY_CONFIG')


@app.task
def simulate(sim_class_name, init_kwargs, run_kwargs_list,
             save_location=RESULT_TO_BACKEND, save_kwargs=None):
    """
    Task for running simulations with Celery_

    Prerequisites:

      * You will need to install the "celery" library, and the "raven"
        library if you wish to use Sentry_ for logging
      * A task broker is required, which is a centralized queue for storing
        tasks. Consult the Celery_ documentation_ for examples.

    Configuration (using environment variables):

      * Set the environment variable PYSB_CELERY_URL with the URL of a task
        broker, e.g. amqp://example.com:5672/myqueue
      * (Optional) Set the environment variable PYSB_SENTRY_DSN with a
        Sentry_ DSN, which allows centralized, aggregated logging of task
        messages and errors (requires the "raven" library)
      * (Optional) Customize the Celery_ configuration by specifying a
        configuration module (Python file) using the environment variable
        PYSB_CELERY_CONFIG

    Starting workers:

      * Start one or more Celery workers, which process tasks::

        celery -A pysb.simulator.tasks worker

    Create simulation tasks:

      * Example of a simulation task::

        >>> from pysb.simulator.tasks import simulate
        >>> from pysb.examples import earm_1_0
        >>> import numpy as np
        >>> simulation_init_kwargs = {                    \
                'tspan': np.linspace(0, 100, 101),        \
                'model': earm_1_0.model                   \
            }

        Empty run() arguments for a single simulation
        >>> run_kwargs = [{}]
        >>> simulate.delay('ScipyOdeSimulator',                 \
                           init_kwargs=simulation_init_kwargs,  \
                           run_kwargs_list=run_kwargs,          \
                           save_location='/tmp/example.h5'      \
            )

    .. _Sentry: http://www.sentry.io/
    .. _Celery: http://www.celeryproject.org/
    .. _documentation: http://docs.celeryproject.org/en/latest/index.html

    Parameters
    ----------
    sim_class_name: str
        Name of a Simulator class, e.g. 'ScipyOdeSimulator'
    init_kwargs: dict
        Dictionary of keyword arguments to the Simulator constructor
    run_kwargs_list: list of dict
        List of dictionaries of keyword arguments to the Simulator run()
        method, one for each set of simulations
    save_location: str
        Save the results to an HDF5 file using
        :py:func:`pysb.simulator.SimulationResult.save()` if this is set to a
        string, which is interpreted as a filename. Alternatively, set to the
        constant pysb.simulator.tasks.RESULT_TO_BACKEND to save the list of
        :py:func:`pysb.simulator.SimulationResult` objects to the user's chosen
        Celery result backend
    save_kwargs: dict
        Dictionary of keyword arguments to the
        :py:func:`pysb.simulator.SimulationResult.save()` method

    Returns
    -------
    list of SimulationResult or list of bool
        List of SimulationResult objects if save_location is set to
        RESULT_TO_BACKEND. Otherwise, a list of bool objects which are True
        if the simulation succeeded (the results are saved to disk at the
        location specified by save_location).
    """
    try:
        sim_class = getattr(pysb.simulator, sim_class_name)
    except AttributeError:
        raise ValueError('{} is not a recognized simulator class'.format(
            sim_class_name))

    sim = sim_class(**init_kwargs)
    simres_list = []
    for simset_run_kwargs in run_kwargs_list:
        simres = sim.run(**simset_run_kwargs)
        if save_location == RESULT_TO_BACKEND:
            simres_list.append(simres)
        else:
            if save_kwargs is None:
                save_kwargs = {}
            simres.save(save_location, **save_kwargs)

    return simres_list if simres_list else True
