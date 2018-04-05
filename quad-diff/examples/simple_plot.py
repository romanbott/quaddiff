import os
import sys
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import quaddiff as qd
import numpy as np

quad = qd.QD()

# Build quadratic differential
ZEROS = [
    1.83498611552-0.62068496262j,
    -1.76274513263-4.74555749048j,
    3.41166252546-3.15694630057j
]

SMPL_POLES = [
    -1.03421848056-0.569727324956j,
]

DBL_POLES = [
    -0.189776680622+4.14349669114j,
    -1.75040778904-4.92983287716j,
    1.90158216804+4.30355758953j
]

for zero in ZEROS:
    quad.add_zero(zero)

for pole in SMPL_POLES:
    quad.add_simple_pole(pole)

for pole in DBL_POLES:
    quad.add_double_pole(pole)

# Create points mesh
N = 6
quad.add_points_mesh(N=N)

# Make phase plot
quad.phase_plot(1, save='/tmp/phase_plot_qd.png', show=False)
