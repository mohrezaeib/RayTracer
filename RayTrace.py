import numpy as np
import scipy as sp
import pickle  # Python 3 uses pickle instead of cPickle
import copy as cp
import RTUtils as RTU
import os
import datetime
import gc
import Elements
import Geometries
import OpticalObject
import PhotonObject as PhObj
import DrawGraph
import Materials
import IOOperation as IOO

class RayTrace():
    def __init__(self):
        """
        Initializes the Variables, Loads some Data and Documents the Structure
        of the Arrays used for saving. Also loads a simple testing geometry
        that is later overwritten by loading an actual geometry.
        """
        # Define media
        air = Materials.MaterialObj("air", 1.0, color='red')
        pva = Materials.MaterialObj(
            "pva",
            1.51,
            IsAbsorber=False,
            conc_a=1e-3,
            conc_b=5e-3,
            color='green',
            IsAligned=True,
            AlignVect=[0, 0, 1.0],
            AlignS=0.785
        )

        """
        Define Elements, only planes are possible up to this point
        Line: Support Vector and Direction Vector
        Plane: Support Vector and NormalVector
        Cylinder: Two points
        """
        BigBox = RTU.MakeBox(M=[0, 0, 0], Size=[1, 1, 1], n_x=[1, 0, 0])
        SmallBox = RTU.MakeBox(
            M=[0.5, 0.5, 0.0],
            Size=[0.1, 0.2, 0.4],
            n_x=[1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]
        )
        Triang = RTU.MakeTriangle(
            [-0.7, -0.1, 0.5],
            [-0.8, 0.6, -0.4],
            [-0.4, 0.7, -0.2]
        )

        self.Elements = np.append(BigBox, SmallBox, axis=0)
        self.Elements = np.append(
            self.Elements,
            np.array([
                [[-0.1, 0, 0], [1, 0, 0]],
                [[-0.2, 0, 0], [1, 0, 0]],
                [[0.3, -0.8, -0.5], [0.8, -0.2, 0.5]],
                Triang,
                [[-0.5, -0.3, 0.0], [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]]
            ]),
            axis=0
        )

        """
        Properties of elements [Type, m1, m2]
        Type: 0: NoElement,JustPosition 1:Mirror 2:Interface 3:ThinLens 4:Dump
        m1/m2: Media on each side of interface if interface
        If Lens, second parameter is focal length
        If Dump, second parameter is # of dumped Photons (0 at the start),
        third parameter = weighted number of photons (accounting for
        lost energy by fluorescence), fourth parameter = spectrum for QE
        """
        self.ElementProperties = [
            [1], [1], [1], [1], [1], [1],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [2, air, pva],
            [4, 0, 0, 0, 1],
            [1],
            [3, 0.3]
        ]

        """
        Element geometry:
        1: Rectangular plane
        2: Circular plane
        3: Cylinder
        4: Triangle
        """
        self.ElementGeometry = [
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
            1, 1, 3, 4, 2
        ]

        """
        Element size parameters:
         - For rectangular plane: [height, width] (plus optional normal vector)
         - For circle or cylinder: [radius]
         - For triangle: [Point1, Point2, Point3]
        """
        self.ElementSize = [
            [1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1],
            [0.4, 0.2, [1.0 / np.sqrt(2), -1.0 / np.sqrt(2), 0]],
            [0.4, 0.2, [1.0 / np.sqrt(2), -1.0 / np.sqrt(2), 0]],
            [0.1, 0.4, [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]],
            [0.1, 0.4, [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]],
            [0.1, 0.2, [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]],
            [0.1, 0.2, [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]],
            [1, 1],
            [1, 1],
            [0.3],
            [[-0.7, -0.1, 0.5], [-0.8, 0.6, -0.4], [-0.4, 0.7, -0.2]],
            [0.2]
        ]

        print(
            np.shape(self.Elements),
            np.shape(self.ElementProperties),
            np.shape(self.ElementGeometry),
            np.shape(self.ElementSize)
        )

        """
        Photon array data structure (legacy):
         - [0]: [Location, Direction]
         - [1]: Photon status
         - [2]: Current wavelength
         - [3]: Current medium
         - [4]: Distance left to travel in absorber
         - [5]: Electric field vector
         - [6]: Starting wavelength of photon
        """
        self.PhArray = [
            np.array([
                [-0.5, -0.5, -0.1],
                RTU.GenerateDirectionSine()
            ]),
            'Start',
            450.0,
            air,
            0,
            [1, 0, 0],
            450.0
        ]
        self.PhArray[0][1] /= np.linalg.norm(self.PhArray[0][1])

        self.Photon = PhObj.Photon(self.PhArray, StartingMaterial=air)

        self.Photon.sgn = 1
        # Number of photons to simulate
        self.NumberOfPhotons = 1

        # Assign colors and opacities for plotting
        self.lineclr = self.Photon.StartingMaterial.color

        # Initialize IOOperation class
        self.IOOperator = IOO.LoadSaveRT("", self)

        # Dictionary of materials
        self.MaterialDict = {}

        self.CreateOpticalElementList()

        # Initialize random seed
        np.random.seed()

    def CreateOpticalElementList(self):
        """
        Create the optical object list from the parsed array data.
        """
        self.OpticalElements = []

        shapes = ["None", "Rectangle", "Circle", "Cylinder", "Triangle"]
        types = ["None", "Mirror", "Interface", "Lens", "Dump"]

        for i in range(len(self.ElementGeometry)):
            SupportVector = self.Elements[i][0]
            DirectionVector = self.Elements[i][1]

            Shape = shapes[self.ElementGeometry[i]]
            if Shape == "Circle":
                Geometry = Geometries.Circle(
                    self.ElementSize[i][0],
                    SupportVector,
                    DirectionVector,
                    self.Photon
                )
            elif Shape == "Cylinder":
                Geometry = Geometries.Cylinder(
                    self.ElementSize[i][0],
                    SupportVector,
                    DirectionVector,
                    self.Photon
                )
            elif Shape == "Rectangle":
                if len(self.ElementSize[i]) == 2:
                    NV = [0, 0, 0]
                else:
                    NV = self.ElementSize[i][2]
                Geometry = Geometries.Rectangle(
                    self.ElementSize[i][0],
                    self.ElementSize[i][1],
                    SupportVector,
                    DirectionVector,
                    self.Photon,
                    NV
                )
            elif Shape == "Triangle":
                Geometry = Geometries.Triangle(
                    self.ElementSize[i][0],
                    self.ElementSize[i][1],
                    self.ElementSize[i][2],
                    self.Photon
                )
            else:
                Geometry = None

            Type = types[self.ElementProperties[i][0]]
            ElementType = Elements.Mirror(self.Photon)
            if Type == "Interface":
                ElementType = Elements.Interface(
                    self.ElementProperties[i][1],
                    self.ElementProperties[i][2],
                    self.Photon
                )
            elif Type == "Lens":
                ElementType = Elements.Lens(
                    self.ElementProperties[i][1],
                    self.Photon
                )
            elif Type == "Dump":
                ElementType = Elements.Dump(
                    self.ElementProperties[i][4],
                    self.Photon
                )

            self.OpticalElements.append(
                OpticalObject.OpticalObject(ElementType, Geometry)
            )

    def CheckElementsForCrossingPoint(self):
        """
        Finds the closest valid element that intersects with the current trajectory.
        Sets self.FirstEl to the intersected element and updates the photon's
        intersection point (PointCy).
        """
        self.SmallestDistance = np.inf
        n = 0
        self.el_n = 0
        self.Photon.PointCy = [0, 0, 0]
        for el in self.OpticalElements:
            point = el.Geometry.FindIntersections()
            # If intersection is valid
            if np.all(point != [np.nan, np.nan, np.nan]):
                dist = np.linalg.norm(self.Photon.Location - point)
                # Accept intersection if it's smaller than current best
                # and not the exact current location.
                if dist < self.SmallestDistance and dist > 1e-8:
                    # If last position was on an element,
                    # check if the new intersection is valid.
                    if self.Photon.WasElement:
                        # Check if the photon is supposed to cross or not.
                        if self.OldElement.Type.IsValidNextPoint(
                            self.OldPoint, point
                        ):
                            self.SmallestDistance = dist
                            self.FirstEl = el
                            self.Photon.PointCy = point
                    else:
                        # For newly emitted photons
                        nzcomp = np.argwhere(
                            np.array(self.Photon.Direction) != 0
                        )[0][0]
                        # Ensure intersection is in forward direction
                        if (
                            (point[nzcomp] - self.Photon.Location[nzcomp])
                            / self.Photon.Direction[nzcomp]
                        ) > 0:
                            self.SmallestDistance = dist
                            self.FirstEl = el
                            self.Photon.PointCy = point
            n += 1

        self.OldElement = self.FirstEl

    def runSimulation(self):
        """
        Simulation Main Loop:
          1) Checks elements for crossing points
          2) Checks if the photon is absorbed
          3) Executes element operations (reflection, refraction, etc.)
        """
        Draw = DrawGraph.DrawElements(self.OpticalElements)
        Draw.DrawElementsPoly()

        # Generate the first photon
        self.Photon.GeneratePhoton()
        self.Photon.BeenAbsorbed = [False, 0]

        for c_Ph in range(self.NumberOfPhotons):
            # Clear plot for each new photon after 7 photons, to avoid clutter
            print("#################  " + str(c_Ph + 1) +
                  " / " + str(self.NumberOfPhotons) +
                  "  ###################")

            self.OldPoint = self.Photon.Location
            self.Photon.ResetPhoton(self.Elements[0])

            if c_Ph % 7 == 0:
                Draw.DrawElementsPoly()

            active = True
            i = 0
            while active:
                self.CheckElementsForCrossingPoint()
                self.Photon.CheckForAbsorption(self.OldPoint, self.SmallestDistance)
                self.OldPoint = self.Photon.OldPoint

                # If the photon isn't absorbed, proceed
                if self.Photon.Status == "Start":
                    self.Photon.OperationPlane = self.OldElement.Geometry.OperationPlane()
                    self.OldElement.Type.ExecuteOperation()
                    self.OldPoint = self.Photon.OldPoint

                # Draw the photon's path
                Draw.DrawLine(self.OldPoint, self.Photon.Location, self.Photon.LineColor)

                # Stop if photon is lost/converted/dumped
                if self.Photon.Status in ["InternalConverted", "Dumped"] or i >= 1e4:
                    active = False
                i += 1

            # Check if photon was absorbed
            if self.Photon.BeenAbsorbed[0]:
                self.Photon.BeenAbsorbed[1] += 1
            self.Photon.BeenAbsorbed[0] = False

            gc.collect()

        self.IOOperator.SaveLog()
