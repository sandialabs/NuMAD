import unittest
import numpy as np
import matplotlib.pyplot as plt

from pynumad.utils.interpolation import interpolator_wrap

class TestInterpolator(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.x = np.linspace(0,10,5)
        v1 = np.sin(self.x)
        v2 = self.x**2
        self.v = np.array((v1,v2)).transpose()
        self.xq = np.linspace(0,10,10)
    
    def test_linear(self):
        vq = interpolator_wrap(self.x,self.v,self.xq,method='linear')
        correct_vq = np.array([
        [ 0.00000000e+00,  0.00000000e+00],
        [ 2.65987620e-01,  2.77777778e+00],
        [ 5.31975239e-01,  5.55555556e+00],
        [ 7.93400045e-02,  1.25000000e+01],
        [-6.12836182e-01,  2.08333333e+01],
        [-5.37385552e-01,  3.19444444e+01],
        [ 3.05691893e-01,  4.58333333e+01],
        [ 7.73330967e-01,  6.11111111e+01],
        [ 1.14654928e-01,  8.05555556e+01],
        [-5.44021111e-01,  1.00000000e+02]
        ])
        array_equal = np.isclose(vq,correct_vq).all()
        self.assertTrue(array_equal)


    def test_pchip(self):
        vq = interpolator_wrap(self.x,self.v,self.xq,method='pchip')
        correct_vq = np.array([
        [  0.        ,   0.        ],
        [  0.47952836,   1.57750343],
        [  0.59634519,   5.21262003],
        [  0.1947027 ,  10.76388889],
        [ -0.76238042,  19.843107  ],
        [ -0.71953191,  30.69415866],
        [  0.4462048 ,  44.48302469],
        [  0.92197998,  60.40237769],
        [  0.50904307,  78.89803384],
        [ -0.54402111, 100.        ]
        ])
        array_equal = np.isclose(vq,correct_vq).all()
        self.assertTrue(array_equal)
        pass

    def test_spline(self):
        vq = interpolator_wrap(self.x,self.v,self.xq,method='spline')
        correct_vq = np.array([
        [  0.        ,   0.        ],
        [  1.09041923,   1.2345679 ],
        [  0.79794503,   4.9382716 ],
        [ -0.11178833,  11.11111111],
        [ -0.87314655,  19.75308642],
        [ -0.75463378,  30.86419753],
        [  0.19006212,  44.44444444],
        [  1.12206955,  60.49382716],
        [  1.16837848,  79.01234568],
        [ -0.54402111, 100.        ]])
        array_equal = np.isclose(vq,correct_vq).all()
        self.assertTrue(array_equal)
        pass


if __name__ == "__main__":
    unittest.main()
