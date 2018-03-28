"""Quadratic Differential unit tester."""
import unittest  # pylint: disable=wrong-import-position

import os
import sys
import time
import cmath as cm
import numpy as np
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import quaddiff as qd  # pylint: disable=wrong-import-position
import quaddiff.core.constants as constants


class QuadraticDifferentialTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()  # pylint: disable=invalid-name

    def test_non_empty_quaddiff(self):
        self.qd.add_zero(1)
        self.assertEqual(-1, self.qd(0))

    def test_quaddiff_zeros(self):
        self.qd.zeros += [1, 2, 3, 4, 5, 6, 7]
        for zero in self.qd.zeros:
            self.assertEqual(0, self.qd(zero))

    def test_poles(self):
        self.qd.add_smplpole(1)
        self.assertEqual(self.qd(1), qd.INF)

        self.qd.add_dblpole(2)
        self.assertEqual(self.qd(2), qd.INF)

    def test_proximity(self):
        self.qd.add_smplpole(1)
        self.qd.add_dblpole(-1)
        close = constants.CLOSE_2POLE
        self.assertFalse(self.qd.close_2pole(2, close))
        self.assertTrue(self.qd.close_2pole(1.001, close))
        self.assertTrue(self.qd.close_2pole(-1.001, close))
        self.qd.sensitivity = 1
        self.assertTrue(self.qd.close_2pole(1, close))

    def test_phase_change(self):
        self.qd.add_zero(1)
        first = self.qd(1 + 1j)

        self.qd.phase = 1j
        second = self.qd(1 + 1j)

        self.assertNotEqual(first, second)
        self.assertEqual(1j*first, second)

    def test_save_and_load(self):
        self.qd.add_zero(1)
        self.qd.add_zero(0)
        self.qd.add_zero(-1j)

        self.qd.add_smplpole(1+1j)
        self.qd.add_smplpole(-2 + 1j)
        self.qd.add_smplpole(3j)
        self.qd.add_smplpole(2 - 1j)

        self.qd.add_dblpole(2)
        self.qd.add_dblpole(3 + 0.5j)

        self.qd.phase = cm.exp(.45 * cm.pi * 2j)

        qd_repr = str(self.qd)

        name = 'test_quad'
        self.qd.save('/tmp', name=name)

        reconstructed = qd.QuadraticDifferential.from_file('/tmp', name=name)
        reconstructed_repr = str(reconstructed)

        self.assertEqual(qd_repr, reconstructed_repr)
        os.remove(os.path.join('/tmp', name + '.json'))


class MonodromyTests(unittest.TestCase):
    def setUp(self):
        self.mono = qd.Monodromy(1+0*1j)

    def test_trivial_dist(self):
        self.assertEqual(self.mono(1), 1)

    def test_monodromy_change(self):
        point = 1 + 0*1j
        dist1 = self.mono(point)

        trajectory = [cm.exp(t*cm.pi*2j) for t in np.linspace(0, 1, 10)]
        for point in trajectory:
            self.mono.update(point)

        dist2 = self.mono(point)

        self.assertNotEqual(dist1, dist2)

    def test_monodromy_return(self):
        point = 1 + 0 * 1j
        dist1 = self.mono(point)

        trajectory = [cm.exp(t*cm.pi*1j) for t in np.linspace(0, 4, 100)]
        for tpoint in trajectory:
            self.mono.update(tpoint)

        dist2 = self.mono(point)

        self.assertEqual(dist1, dist2)

    def test_monodromy_continuous(self):
        point = 1 + 0 * 1j
        dist_old = self.mono(point)

        trajectory = [cm.exp(t*cm.pi*1j) for t in np.linspace(0, 8, 500)]

        is_close = True
        for point in trajectory:
            dist = self.mono(point)
            if abs(dist-dist_old) > 0.2:
                is_close = False
                break
            dist_old = dist

        self.assertTrue(is_close)

    def test_monodromy_continuous_brownian(self):
        point = 1 + 0 * 1j
        dist_old = self.mono(point)
        phase = 0.2

        is_close = True
        for k in range(50000):
            phase += np.random.normal(scale=0.08)
            point = cm.rect(1, phase)
            dist = self.mono(point)
            if abs(dist-dist_old) > 0.2:
                is_close = False
                break
            dist_old = dist

        self.assertTrue(is_close)

    def test_monodromy_continuous_2dbrownian(self):
        point = 1 + 1j
        dist_old = self.mono(point)
        real = 1
        imag = 0.2

        is_close = True
        for n in range(20):
            for k in range(5000):
                real = np.random.normal(scale=0.05)
                real = min(real, 3.1)
                imag += np.random.normal(scale=0.15)
                point = cm.exp(complex(real, imag))
                dist = self.mono(point)
                is_close = is_close and abs(dist-dist_old) < 0.5
                if abs(dist-dist_old) > 0.5:
                    print(abs(dist-dist_old), point)
                dist_old = dist

        self.assertTrue(is_close)


if __name__ == '__main__':
    unittest.main(verbosity=2)
