from pysb.importers.boolean import model_from_boolean
from pysb.bng import run_ssa
import numpy as np
import matplotlib.pyplot as plt

model = model_from_boolean('ABC_example.txt')

tspan = np.linspace(0, 10, 101)

n_sims = 1000
traj = []
for n in range(n_sims):
    print n
    traj.append(run_ssa(model, t_end=tspan[-1], n_steps=len(tspan)-1, verbose=False))

# transpose to aid in plotting
a_traj = np.array([tr['A_True_obs'] for tr in traj]).T
b_traj = np.array([tr['B_True_obs'] for tr in traj]).T
c_traj = np.array([tr['C_True_obs'] for tr in traj]).T

# plot the mean concentrations at each time point
plt.plot(tspan, a_traj.mean(axis=1), lw=2, label="A")
plt.plot(tspan, b_traj.mean(axis=1), lw=2, label="B")
plt.plot(tspan, c_traj.mean(axis=1), lw=2, label="C")

plt.xlabel('Time (au)')
plt.ylabel('Frequency')
plt.legend(loc=0)
 
plt.show()
