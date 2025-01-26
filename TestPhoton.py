import unittest
import PhotonObject
import RayTrace
import Materials
import numpy as np


class MyTestRT(unittest.TestCase):
    def setUp(self):
        """
        Initializes a RayTrace instance and parses data from the './Job/' directory.
        """
        self.RT = RayTrace.RayTrace()
        self.RT.ParsePath("./Job/")

    def test_ParsePath(self):
        """
        Test if the parsing of the job path is correct, checking the number
        of photons and certain materials.
        """
        self.assertEqual(7, self.RT.NumberOfPhotons)
        self.assertIsInstance(self.RT.air, Materials.MaterialObj)
        self.assertIsInstance(self.RT.blue, Materials.MaterialObj)

    def test_GeneratePhoton(self):
        """
        Test if photon generation uses expected initial coordinates.
        """
        self.RT.Photon.GeneratePhoton()

        self.assertEqual(self.RT.Photon.Location[0], 0.0)
        self.assertEqual(self.RT.Photon.Location[1], 0.0)
        self.assertEqual(self.RT.Photon.Location[2], 0.1)

    def test_FindIntersections(self):
        """
        Test intersection calculations with specific optical elements.
        """
        self.RT.Photon.GeneratePhoton()

        # Check intersection sets match expected results
        self.assertEqual(
            set(self.RT.OpticalElements[2].Geometry.FindIntersections()),
            set([0.0, 0.0, 0.2])
        )
        self.assertEqual(
            set(self.RT.OpticalElements[7].Geometry.FindIntersections()),
            set([0.0, 0.0, 0.0])
        )
        self.assertEqual(
            set(self.RT.OpticalElements[0].Geometry.FindIntersections()),
            set([np.nan, np.nan, np.nan])
        )

    def test_ReflectOnPlane(self):
        """
        Test reflection of photon on a plane with different directions.
        """
        self.RT.Photon.GeneratePhoton()
        self.RT.Photon.BeenAbsorbed = [False, 0]

        # Prepare for first reflection test
        self.RT.OldPoint = self.RT.Photon.Location
        self.RT.Photon.ResetPhoton(self.RT.Elements[0])

        self.RT.CheckElementsForCrossingPoint()
        self.RT.Photon.CheckForAbsorption(self.RT.OldPoint, self.RT.SmallestDistance)
        self.RT.OldPoint = self.RT.Photon.OldPoint

        self.RT.Photon.OperationPlane = self.RT.OldElement.Geometry.OperationPlane()
        self.RT.Photon.ReflectOnPlane()

        # Check final location and direction match expectations
        self.assertEqual(set(self.RT.Photon.Location), set([0.0, 0.0, 0.0]))
        self.assertEqual(set(self.RT.Photon.Direction), set([0.0, 0.0, -1.0]))

        # Second reflection test with a different direction
        self.RT.Photon.ResetPhoton(self.RT.Elements[0])
        self.RT.Photon.Direction = np.array([1.0, 1.0, -1.0])
        self.RT.CheckElementsForCrossingPoint()
        self.RT.Photon.CheckForAbsorption(self.RT.OldPoint, self.RT.SmallestDistance)
        self.RT.OldPoint = self.RT.Photon.OldPoint
        self.RT.Photon.OperationPlane = self.RT.OldElement.Geometry.OperationPlane()
        self.RT.Photon.ReflectOnPlane()

        # Check final location is near the expected coordinates
        self.assertAlmostEqual(self.RT.Photon.Location[0], 0.09090909, 0)
        self.assertAlmostEqual(self.RT.Photon.Location[1], 0.09090909, 0)
        self.assertAlmostEqual(self.RT.Photon.Location[2], 0.09090909, 0)

        # Check direction is flipped properly
        # The set of components is either [-1,1,1] or [1,-1,-1].
        direction_set = set(self.RT.Photon.Direction)
        self.assertTrue(
            direction_set == set([-1.0, 1.0, 1.0]) or
            direction_set == set([1.0, -1.0, -1.0])
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
