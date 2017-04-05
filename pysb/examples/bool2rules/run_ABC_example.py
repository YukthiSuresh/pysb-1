from pysb.importers.boolean import model_from_boolean
from pysb.bng import run_ssa
import numpy as np
import matplotlib.pyplot as plt

model_gsp = model_from_boolean('ABC_example.txt', mode='GSP')
model_roa = model_from_boolean('ABC_example.txt', mode='ROA')

n_sims = 1000

###### Gillespie updating #####
tspan = np.linspace(0, 10, 101)
traj = []
for n in range(n_sims):
    print n
    x = run_ssa(model_gsp, t_end=tspan[-1], n_steps=len(tspan)-1, verbose=False)
    traj.append(x)

# plot the mean concentrations at each time point
plt.figure('Gillespie updating')

# transpose to aid in plotting
a_traj = np.array([tr['A_True_obs'] for tr in traj]).T
b_traj = np.array([tr['B_True_obs'] for tr in traj]).T
c_traj = np.array([tr['C_True_obs'] for tr in traj]).T

plt.plot(tspan, a_traj.mean(axis=1), lw=2, label="A")
plt.plot(tspan, b_traj.mean(axis=1), lw=2, label="B")
plt.plot(tspan, c_traj.mean(axis=1), lw=2, label="C")

plt.xlabel('Time (au)')
plt.ylabel('Frequency')
plt.legend(loc=0)

##### Random order asynchronous updating #####
rounds = 10
nodes = len(model_roa.monomers)-1 # exclude RESET Monomer
traj = []
for n in range(n_sims):
    print n
    x = run_ssa(model_roa, t_end=nodes*rounds*100, n_steps=1, verbose=False, 
                output_step_interval=(nodes*2+2), 
                max_sim_steps=((nodes*2+2)*rounds))
    traj.append(x)

# plot the mean concentrations at each update round
plt.figure('Random order asynchronous updating')

# transpose to aid in plotting
a_traj = np.array([tr['A_True_obs'] for tr in traj]).T
b_traj = np.array([tr['B_True_obs'] for tr in traj]).T
c_traj = np.array([tr['C_True_obs'] for tr in traj]).T

plt.plot(range(rounds+1), a_traj.mean(axis=1), '-o', lw=2, ms=12, label='A')
plt.plot(range(rounds+1), b_traj.mean(axis=1), '-o', lw=2, ms=12, label='B')
plt.plot(range(rounds+1), c_traj.mean(axis=1), '-o', lw=2, ms=12, label='C')

plt.xlabel('Round')
plt.ylabel('Frequency')
plt.legend(loc=0)

plt.show()
