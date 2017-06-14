from earm.lopez_embedded import model
from pysb.tools.helper_functions import pre_equilibration as equil
import numpy as np

time = np.linspace(0, 20000, 100)

a = equil(model, time, ligand='L_0', tolerance=1e-3)
