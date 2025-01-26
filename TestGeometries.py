import unittest
import PhotonObject
import RayTrace
import Materials
import Geometries
import numpy as np
import DrawGraph
import OpticalObject
import Elements


class MyTestRT(unittest.TestCase):
    def setUp(self):
        """
        Set up a RayTrace instance and parse an existing job folder.
        Then generate a photon for tests.
        """
        self.RT = RayTrace.RayTrace()
        self.RT.ParsePath("./Job/")
        self.RT.Photon.GeneratePhoton()

    def test_IntersectRectangle(self):
        """
        FOUND A BUG: Intersection with 'crooked' rectangles not working properly.
        That is known behavior though...
        """
        # First rectangle with normal vector [0.0, 1.0, 0.0],
        # no custom NormalVector argument to Rectangle
        self.Rectangle = Geometries.Rectangle(
            1.5,
            2.5,
            [1.0, 2.0, 3.0],
            [0.0, 1.0, 0.0],
            self.RT.Photon,
            [0.0, 0.0, 0.0]
        )

        Mirror = Elements.Mirror(self.RT.Photon)
        OpticalObj = OpticalObject.OpticalObject(Mirror, self.Rectangle)
        OpticalElements = [OpticalObj]

        # Case 1: The ray should NOT intersect the rectangle
        self.RT.Photon.Location = np.array([0.0, -2.0, 2.0])
        self.RT.Photon.Direction = np.array([2.0, 2.0, 0.0])
        point = OpticalElements[0].Geometry.FindIntersections()

        # The intersection is expected to be [nan, nan, nan]
        self.assertEqual(set(point), set([np.nan, np.nan, np.nan]))

        # Case 2: The ray should intersect at [0.4, 2.0, 2.2]
        self.RT.Photon.Location = np.array([0.0, -2.0, 2.0])
        self.RT.Photon.Direction = np.array([0.2, 2.0, 0.1])
        point = OpticalElements[0].Geometry.FindIntersections()

        # Visualization (draw) code
        self.DrawClass = DrawGraph.DrawElements(OpticalElements)
        self.DrawClass.DrawElementsPoly()
        self.DrawClass.DrawPoint(self.RT.Photon.Location)
        print(point)
        self.DrawClass.DrawPoint(self.RT.Photon.Location + 10 * self.RT.Photon.Direction)
        self.DrawClass.DrawLine(self.RT.Photon.Location, np.array(point), "b")

        self.assertEqual(set([0.4, 2.0, 2.2]), set(point))

        # Second rectangle with a custom NormalVector cross([0.0, 1.0, 0.0], [1, 2, 3])
        self.Rectangle = Geometries.Rectangle(
            1.5,
            2.5,
            [1.0, 2.0, 3.0],
            [0.0, 1.0, 0.0],
            self.RT.Photon,
            list(np.cross([0.0, 1.0, 0.0], [1, 2, 3]))
        )

        Mirror = Elements.Mirror(self.RT.Photon)
        OpticalObj = OpticalObject.OpticalObject(Mirror, self.Rectangle)
        OpticalElements = [OpticalObj]

        self.RT.Photon.Location = np.array([0.0, -2.0, 2.0])
        self.RT.Photon.Direction = np.array([0.2, 2.0, 0.1])
        point = OpticalElements[0].Geometry.FindIntersections()

        self.DrawClass = DrawGraph.DrawElements(OpticalElements)
        self.DrawClass.DrawElementsPoly()
        self.DrawClass.DrawPoint(self.RT.Photon.Location)
        print(point)
        self.DrawClass.DrawPoint(self.RT.Photon.Location + 10 * self.RT.Photon.Direction)
        self.DrawClass.DrawLine(self.RT.Photon.Location, np.array(point), "b")

        self.assertEqual(set([0.4, 2.0, 2.2]), set(point))

        # Uncomment if you want to pause execution to view the plots in an interactive session
        # input("Press Enter to continue...")


if __name__ == "__main__":
    unittest.main(verbosity=2)
