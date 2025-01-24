import unittest
import PhotonObject
import RayTrace
import Materials
import numpy as np


class MyTestRT(unittest.TestCase):
    def setUp(self):
        self.RT = RayTrace.RayTrace()
        self.RT.ParsePath("./Job/")

    def test_ParsePath(self):
        

        self.assertEqual(7,self.RT.NumberOfPhotons)
        self.assertIsInstance(self.RT.air,Materials.MaterialObj)
        self.assertIsInstance(self.RT.blue,Materials.MaterialObj)

    def test_GeneratePhoton(self):
        self.RT.Photon.GeneratePhoton()

        self.assertEqual(self.RT.Photon.Location[0],0.0)
        self.assertEqual(self.RT.Photon.Location[1],0.0)
        self.assertEqual(self.RT.Photon.Location[2],0.1)

    def test_FindIntersections(self):
        self.RT.Photon.GeneratePhoton()

        self.assertEqual(set(self.RT.OpticalElements[2].Geometry.FindIntersections()),set([0.0,0.0,0.2]))
        self.assertEqual(set(self.RT.OpticalElements[7].Geometry.FindIntersections()),set([0.0,0.0,0.0]))
        self.assertEqual(set(self.RT.OpticalElements[0].Geometry.FindIntersections()),set([np.nan, np.nan, np.nan]))


    def test_ReflectOnPlane(self):
        self.RT.Photon.GeneratePhoton()        
        self.RT.Photon.BeenAbsorbed = [False, 0]

        self.RT.OldPoint = self.RT.Photon.Location
        self.RT.Photon.ResetPhoton(self.RT.Elements[0])

        self.RT.CheckElementsForCrossingPoint()
        self.RT.Photon.CheckForAbsorption(self.RT.OldPoint, self.RT.SmallestDistance)
        self.RT.OldPoint = self.RT.Photon.OldPoint

        self.RT.Photon.OperationPlane = self.RT.OldElement.Geometry.OperationPlane()

        self.RT.Photon.ReflectOnPlane()

        self.assertEqual(set(self.RT.Photon.Location),set([0.0,0.0,0.0]))
        self.assertEqual(set(self.RT.Photon.Direction),set([0.0,0.0,-1.0]))



        self.RT.Photon.ResetPhoton(self.RT.Elements[0])
        self.RT.Photon.Direction = np.array([1.0,1.0,-1.0])
        self.RT.CheckElementsForCrossingPoint()
        self.RT.Photon.CheckForAbsorption(self.RT.OldPoint, self.RT.SmallestDistance)
        self.RT.OldPoint = self.RT.Photon.OldPoint
        self.RT.Photon.OperationPlane = self.RT.OldElement.Geometry.OperationPlane()

        self.RT.Photon.ReflectOnPlane()

        self.assertAlmostEqual(self.RT.Photon.Location[0],0.09090909,0)
        self.assertAlmostEqual(self.RT.Photon.Location[1],0.09090909,0)
        self.assertAlmostEqual(self.RT.Photon.Location[2],0.09090909,0)

        self.assertTrue((set(self.RT.Photon.Direction) == set([-1.0,1.0,1.0])) or (set(self.RT.Photon.Direction) == set([1.0,-1.0,-1.0])))
        

        # self.RT.OldElement.Type.ExecuteOperation()
        # self.RT.OldPoint = self.RT.Photon.OldPoint





if __name__=="__main__":
    unittest.main(verbosity=2)