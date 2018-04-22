"""Testing refactored trayectory solver"""
# pylint: disable=no-member
# pylint: disable=wrong-import-position
# pylint: disable=import-error

from __future__ import print_function

import os
import sys
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from time import time
import numpy as np
from tqdm import tqdm

from quaddiff.core.trajectory import calculate_ray
from quaddiff.core.trajectory import calculate_ray_2
import quaddiff as qd

N_ZEROS = 4
N_SMPL_POLES = 4
N_DBL_POLES = 6
N_LOOPS = 300


def main():
    """Calculate mean processing times for solver methods"""
    quad = qd.QuadraticDifferential()

    opt_time = 0
    old_time = 0
    opt_errs = 0
    old_errs = 0
    for _ in tqdm(range(N_LOOPS)):
        # build quadratic differential
        zeros = np.random.uniform(-5, 5, size=[N_ZEROS, 2])
        smpl_poles = np.random.uniform(-5, 5, size=[N_SMPL_POLES, 2])
        dbl_poles = np.random.uniform(-5, 5, size=[N_DBL_POLES, 2])

        for zero in zeros:
            quad.add_zero(complex(*zero))

        for pole in smpl_poles:
            quad.add_smplpole(complex(*pole))

        for pole in dbl_poles:
            quad.add_dblpole(complex(*pole)) #bla bal

        plotpoint = complex(*np.random.uniform(-5, 5, size=[2]))

        start = time()
        try:
            calculate_ray(plotpoint, quad)
            opt_time += time() - start
        except:
            opt_errs += 1

        start = time()
        try:
            calculate_ray_2(plotpoint, quad)
            old_time += time() - start
        except:
            old_errs += 1

    # Message header
    message = ("=" * 19) + "RESULTS" + ("=" * 19) + "\n"
    message += "|{:^10}|{:^10}|{:^10}|{:^10}|\n".format("version", "total", "per", "errors")
    message += "|" + "-" * 43 + "|" + "\n"

    # New version information
    opt_per = opt_time / (N_LOOPS - opt_errs)
    message += "|{:^10}|{:^10.4f}|{:^10.4f}|{:^10}|\n".format("new", opt_time, opt_per, opt_errs)
    message += "|" + "-" * 43 + "|" + "\n"

    # Old version information
    old_per = old_time / (N_LOOPS - old_errs)
    message += "|{:^10}|{:^10.4f}|{:^10.4f}|{:^10}|\n".format("old", old_time, old_per, old_errs)
    message += "-" * 45 + "\n"
    message += "N_LOOPS = {}".format(N_LOOPS)

    print(message)

if __name__ == "__main__":
    main()
