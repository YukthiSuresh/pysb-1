from __future__ import print_function
from pysb.examples.earm_1_0 import model
from pysb.tools.model_visualization import run_visualization
import numpy as np

tspan = np.linspace(0, 20160, 100)
run_visualization(model, tspan, render_type='flux', save_video=True)

