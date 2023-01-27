import unittest

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.utils.mesh_utils import getShellMesh
from pynumad.objects.Blade import Blade


class TestMesh(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        yml_file = 'data/BAR1_SNL_1_18_2021.yaml'
        self.yml_blade = Blade(yml_file)
        self.adhes = 1

    def test_shell_mesh(self):
        nodes,elements,OSSets,SWSets,adNds,adEls = self.yml_blade.getShellMesh(self.adhes); 