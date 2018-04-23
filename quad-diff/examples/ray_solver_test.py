import os
import sys
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import cProfile
import pstats
import StringIO

from quaddiff.core.quaddiff import QuadraticDifferential
from quaddiff.core.trajectory import (
    calculate_ray, calculate_ray_2, calculate_ray_3)


quad = QuadraticDifferential()

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

quad.zeros = ZEROS
quad.smplpoles = SMPL_POLES
quad.dblpoles = DBL_POLES

init_point = 0

pr1 = cProfile.Profile()
pr1.enable()
t1 = calculate_ray(init_point, quad)
pr1.disable()
with open('ray_solver_1.txt', 'w') as stream:
    ps = pstats.Stats(pr1, stream=stream)
    ps.sort_stats('tottime').print_stats()

pr2 = cProfile.Profile()
pr2.enable()
t2 = calculate_ray_2(init_point, quad)
pr2.disable()
with open('ray_solver_2.txt', 'w') as stream:
    ps = pstats.Stats(pr2, stream=stream)
    ps.sort_stats('tottime').print_stats()

pr3 = cProfile.Profile()
pr3.enable()
t3 = calculate_ray_3(init_point, quad)
pr3.disable()
with open('ray_solver_3.txt', 'w') as stream:
    ps = pstats.Stats(pr3, stream=stream)
    ps.sort_stats('tottime').print_stats()


