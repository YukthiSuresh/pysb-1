#!/usr/bin/env python
"""Simulate the bax_pore model and plot the results."""

from __future__ import print_function
from pylab import *
from pysb.simulator.scipy import ScipyOdeSimulator

from bax_pore import model


t = linspace(0, 100)
print("Simulating...")
x = ScipyOdeSimulator.execute(model, tspan=t)

plt.plot(t, x['BAX4'])
plt.plot(t, x['BAX4_inh'])
plt.legend(['BAX4', 'BAX4_inh'], loc='upper left')
plt.show()
