import os
import sys
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from time import time
import numpy as np

from quaddiff.core.trajectory import calculate_ray
from quaddiff.core.trajectory import calculate_ray_2
import quaddiff as qd

quad = qd.QuadraticDifferential()

# Number of loops
N_LOOPS = 100

optimized_time = 0
legacy_time = 0
for _ in range(N_LOOPS):
    print("Building new Quadratic Differential")
    # Build quadratic differential
    N_ZEROS = 2
    ZEROS = np.random.uniform(-5, 5, size=[N_ZEROS, 2])
    
    N_SMPL_POLES = 2
    SMPL_POLES  = np.random.uniform(-5, 5, size=[N_SMPL_POLES, 2])
    
    N_DBL_POLES = 2
    DBL_POLES = np.random.uniform(-5, 5, size=[N_DBL_POLES, 2])
    
    for zero in ZEROS:
        quad.add_zero(complex(*zero))
    
    for pole in SMPL_POLES:
        quad.add_smplpole(complex(*pole))
    
    for pole in DBL_POLES:
        quad.add_dblpole(complex(*pole))
    
    plotpoint = complex(*np.random.uniform(-5, 5, size=[2]))

    print("Calculating with new code")
    t = time()
    calculate_ray(plotpoint, quad)
    optimized_time += time() - t

    print("Calculating with old code")
    t = time()
    plotpoint = complex(*np.random.uniform(-5, 5, size=[2]))
    calculate_ray_2(plotpoint, quad)
    legacy_time += time() - t

mesage = ''' ======== RESULTS =======
|{:^10}|{:^10}|{:^10}|{:^10}|
----------------------------------------------
|{:^10}|{:^10}|{:^10}|{:^10}|
----------------------------------------------
'''.format("new-total", "new-per", "old-total", "old-per", round(optimized_time, 5), round(optimized_time / N_LOOPS, 5), round(legacy_time, 5), round(legacy_time / N_LOOPS, 5))
print(mesage)
