import unittest
from os.path import abspath, dirname, join

from pynumad.io.ansys import writeANSYSshellModel
from pynumad.objects.Blade import Blade

testdir = dirname(abspath(str(__file__)))
test_datadir = join(testdir, "test_data")

class TestANSYSWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.yamlfile = join(test_datadir,"blade_yamls","myBlade.yaml")
        
    def step1_read_yaml(self):
        self.blade = Blade(self.yamlfile)

    def step2_build_mesh(self):
        self.blade.mesh = 0.2
        adhes = 1
        self.meshData = self.blade.getShellMesh(includeAdhesive=adhes)

    def step3_write_ANSYS_src(self):
        config = {}
        config["BoundaryCondition"] = 'cantilevered'
        config["elementType"] = '181'
        config["MultipleLayerBehavior"] = 'multiply'
        config["dbgen"] = 1
        config["dbname"] = 'master'

        filename = "myblade_ansys.src"
        includeAdhesive = 0

        writeANSYSshellModel(
            self.blade,
            filename,
            self.meshData,
            config,
            includeAdhesive
        )
    
    def get_steps(self):
        for name in dir(self): #look at all functions
            if name.startswith("step"):
                yield name, getattr(self, name)

    def test_workflow(self):
        for name, step in self.get_steps():
            try:
                step()
            except Exception as e:
                self.fail("{} failed ({}: {})".format(step, type(e), e))


        
