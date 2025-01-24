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
        self.RT = RayTrace.RayTrace()
        self.RT.ParsePath("./Job/")
        self.RT.Photon.GeneratePhoton()


    def test_IntersectRectangle(self):
        #FOUNDABUG! Intersection with crooked rectangles not working porperly. That is known behavior though...
        #         
        self.Rectangle = Geometries.Rectangle(1.5,2.5,[1.0,2.0,3.0],[0.0,1.0,0.0], self.RT.Photon,[0.0,0.0,0.0])
        
        Mirror = Elements.Mirror(self.RT.Photon)
        OpticalObj = OpticalObject.OpticalObject(Mirror,self.Rectangle)
        OpticalElements = [OpticalObj]

        self.RT.Photon.Location = np.array([0.0,-2.0,2.0])
        self.RT.Photon.Direction = np.array([2.0,2.0,0.0])
        point = OpticalElements[0].Geometry.FindIntersections()

        self.assertEqual(set(point),set([np.nan, np.nan, np.nan]))


        self.RT.Photon.Location = np.array([0.0,-2.0,2.0])
        self.RT.Photon.Direction = np.array([0.2,2.0,0.1])
        point = OpticalElements[0].Geometry.FindIntersections()

        self.DrawClass = DrawGraph.DrawElements(OpticalElements)
        self.DrawClass.DrawElementsPoly()        
        self.DrawClass.DrawPoint(self.RT.Photon.Location)
        print point
        self.DrawClass.DrawPoint(self.RT.Photon.Location + 10* self.RT.Photon.Direction)
        self.DrawClass.DrawLine(self.RT.Photon.Location,np.array(point),"b")

        self.assertEqual(set([0.4,2.0,2.2]),set(point))

        #raw_input()



        self.Rectangle = Geometries.Rectangle(1.5,2.5,[1.0,2.0,3.0],[0.0,1.0,0.0], self.RT.Photon,list(np.cross([0.0,1.0,0.0],[1,2,3])))
        
        Mirror = Elements.Mirror(self.RT.Photon)
        OpticalObj = OpticalObject.OpticalObject(Mirror,self.Rectangle)
        OpticalElements = [OpticalObj]

        self.RT.Photon.Location = np.array([0.0,-2.0,2.0])
        self.RT.Photon.Direction = np.array([0.2,2.0,0.1])
        point = OpticalElements[0].Geometry.FindIntersections()


        self.DrawClass = DrawGraph.DrawElements(OpticalElements)
        self.DrawClass.DrawElementsPoly()        
        self.DrawClass.DrawPoint(self.RT.Photon.Location)
        print point
        self.DrawClass.DrawPoint(self.RT.Photon.Location + 10* self.RT.Photon.Direction)
        self.DrawClass.DrawLine(self.RT.Photon.Location,np.array(point),"b")

        self.assertEqual(set([0.4,2.0,2.2]),set(point))

        #raw_input()


if __name__=="__main__":
    unittest.main(verbosity=2)