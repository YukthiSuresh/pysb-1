from __future__ import print_function
# from pysb.examples.earm_1_0 import model
from earm.lopez_embedded import model
from pysb.tools.model_visualization import run_visualization
import numpy as np
import csv
# type Bax cluster: 5400
# type Bak cluster: 4052

# type1 Bid cluster: 3905
# type2 Bid cluster: 2415

f = open('/home/oscar/Documents/tropical_earm_different_rates/parameters_5000/pars_embedded_2415.txt')
data = csv.reader(f)
parames = [float(i[1]) for i in data]
tspan = np.linspace(0, 20160, 100)
run_visualization(model, tspan, parameters=parames, render_type='flux', save_video=False, verbose=True)

