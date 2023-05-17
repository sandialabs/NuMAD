import unittest
from os.path import join

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.utils.mesh_utils import getShellMesh
from pynumad.objects.Blade import Blade


class TestMesh(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        yaml_file = join("data", "BAR1_SNL_1_18_2021.yaml")
        self.blade = Blade(yml_file)
        self.adhes = 1

    def test_shell_mesh(self):
        meshData = self.blade.getShellMesh(self.adhes); 