from pysb.examples.simple_reaction_pyurdme import model
from pysb.tools.pysb_pyurdme import run_pyurdme
import numpy as np
import matplotlib.pyplot as plt
import pyurdme
from pysb.integrate import odesolve



model.diffusivities = [('E(b=None)',0.005), ('S(b=None)',0.003), ('E(b=1) % S(b=1)',0.001), ('P()',0.002)]



mesh = pyurdme.URDMEMesh.generate_unit_square_mesh(40,40)

tspan = np.linspace(0, 5, 500)
result = run_pyurdme(model, tspan, mesh)
result.export_to_vtk("__s0", "simple_output_s0")
result.export_to_vtk("__s1", "simple_output_s1")
result.export_to_vtk("__s2", "simple_output_s2")
result.export_to_vtk("__s3", "simple_output_s3")