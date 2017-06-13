from earm.lopez_embedded import model
from pysb.tools import pre_equilibration as equil
import numpy as np

time = np.linspace(0, 20000, 100)

a = equil.pre_equilibration(model, time, ligand='L_0', tolerance=1e-3)
