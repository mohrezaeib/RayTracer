import numpy as np
import scipy as sp
import cPickle
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
        '''Initializes the Variables, Loads some Data and Documents the Structure of the Arrays used for saving. Also Loads a Simple Testing Geometry that is later overwritten by Loading an actual Geometry'''
        #define Media
        air = Materials.MaterialObj("air", 1.0, color='red')
        pva = Materials.MaterialObj("pva", 1.51, IsAbsorber=False, conc_a=1e-3, conc_b=5e-3, color='green', IsAligned=True, AlignVect=[0, 0, 1.0], AlignS=0.785)

        '''Define Elements, only planes are possible up to this point
        Line: Support Vector and Direction Vector
        Plane: Support Vector and NormalVector
        Cylinder: Two points
        '''
        BigBox = RTU.MakeBox(M=[0, 0, 0], Size=[1, 1, 1], n_x=[1, 0, 0])
        SmallBox = RTU.MakeBox(M=[0.5, 0.5, 0.0], Size=[0.1, 0.2, 0.4], n_x=[1.0/np.sqrt(2), 1.0/np.sqrt(2), 0])
        Triang = RTU.MakeTriangle([-0.7, -0.1, 0.5], [-0.8, 0.6, -0.4], [-0.4, 0.7, -0.2])

        self.Elements = np.append(BigBox, SmallBox, axis=0)
        self.Elements = np.append(self.Elements, np.array([[[-0.1, 0, 0], [1, 0, 0]], [[-0.2, 0, 0], [1, 0, 0]], [[0.3, -0.8, -0.5], [0.8, -0.2, 0.5]], Triang, [[-0.5, -0.3, 0.0], [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0]]]), axis=0)
        '''Properties of elements [Type,m1,m2]
        Type: 0: NoElement,JustPosition 1:Mirror 2:interface 3:ThinLens 4:Dump
        m1/m2: Media on each side of interface if interface
        If Lens, second parameter is focal length
        If Dump, second parameter is number of dumped Photons (0 at the start) third parameter is weighted number of photons taking into account the lost energy by fluorescence etc, fourth parameter is spectrum for quantum efficiency
        '''
        self.ElementProperties = [[1], [1], [1], [1], [1], [1], [2, air, pva], [2, air, pva], [2, air, pva], [2, air, pva], [2, air, pva], [2, air, pva], [2, air, pva], [2, air, pva], [4, 0, 0, 0, 1], [1], [3, 0.3]]
        '''Define Element geometry, at this point: 1:rectancular plane 2:circular plane 3:Cylinder 4:Triangle'''
        self.ElementGeometry = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 2]
        '''Define the Size of the planes in two dimensions from the center point (the axis for the first direction can be defined additionaly for crooked stuff); if circle or cylinder, just one value-The radius; If triangle, pass all three points'''
        self.ElementSize = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [0.4, 0.2, [1.0/np.sqrt(2), -1.0/np.sqrt(2), 0]], [0.4, 0.2, [1.0/np.sqrt(2), -1.0/np.sqrt(2), 0]], [0.1, 0.4, [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0]], [0.1, 0.4, [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0]], [0.1, 0.2, [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0]], [0.1, 0.2, [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0]], [1, 1], [1, 1], [0.3], [[-0.7, -0.1, 0.5], [-0.8, 0.6, -0.4], [-0.4, 0.7, -0.2]], [0.2]]

        print np.shape(self.Elements), np.shape(self.ElementProperties), np.shape(self.ElementGeometry), np.shape(self.ElementSize)

        '''First Element is line equation of photon trajectory, second is the Current Photon status, third is the wavelength, fourth is the current medium and fifth distance left to travel in the absorber,sixth: Electric Field Vector,seventh:Starting Wavelength of photon '''
        self.PhArray = [np.array([[-0.5, -0.5, -0.1], RTU.GenerateDirectionSine()]), 'Start', 450.0, air, 0, [1, 0, 0], 450.0]
        self.PhArray[0][1] /= np.linalg.norm(self.PhArray[0][1])

        self.Photon = PhObj.Photon(self.PhArray, StartingMaterial=air)

        self.Photon.sgn = 1
        '''Number of photons to simulate'''
        self.NumberOfPhotons = 1

        '''Assign Colors and opacities for plotting'''
        self.lineclr = self.Photon.StartingMaterial.color

        self.IOOperator = IOO.LoadSaveRT("",self)

        self.MaterialDict = {}

        self.CreateOpticalElementList()

        np.random.seed()

    def CreateOpticalElementList(self):
        '''Create Object List'''
        self.OpticalElements = []


        '''Translate Parsed Data to Readable Variables, This should be changed once the parser works right or simply be made part of the parsing code'''


        shapes = ["None", "Rectangle", "Circle", "Cylinder", "Triangle"]
        types = ["None", "Mirror", "Interface", "Lens", "Dump"]

        for i in range(len(self.ElementGeometry)):

            SupportVector = self.Elements[i][0]
            DirectionVector = self.Elements[i][1]


            Shape = shapes[self.ElementGeometry[i]]

            if Shape == "Circle":
                 Geometry = Geometries.Circle(self.ElementSize[i][0], SupportVector, DirectionVector,self.Photon)
            elif Shape == "Cylinder":
                Geometry = Geometries.Cylinder(self.ElementSize[i][0], SupportVector, DirectionVector,self.Photon)
            elif Shape == "Rectangle":
                if len(self.ElementSize[i]) == 2:
                    NV = [0,0,0]
                else:
                    NV = self.ElementSize[i][2]
                Geometry = Geometries.Rectangle(self.ElementSize[i][0], self.ElementSize[i][1], SupportVector, DirectionVector,self.Photon,NV)
            elif Shape == "Triangle":
            
                Geometry = Geometries.Triangle(self.ElementSize[i][0],self.ElementSize[i][1],self.ElementSize[i][2],self.Photon)


            Type = types[self.ElementProperties[i][0]]
            ElementType = Elements.Mirror(self.Photon)
            if Type == "Interface":
                ElementType = Elements.Interface(self.ElementProperties[i][1],self.ElementProperties[i][2], self.Photon)
            elif Type == "Lens":
                ElementType = Elements.Lens(self.ElementProperties[i][1], self.Photon)
            elif Type == "Dump":
                ElementType = Elements.Dump(self.ElementProperties[i][4], self.Photon)


            self.OpticalElements.append(OpticalObject.OpticalObject(ElementType,Geometry))

            # self.OpticalElements.append(OptEl.OpticalElement(self.Elements[i][0], self.Elements[i][1], self.ElementProperties[i], self.ElementGeometry[i], self.ElementSize[i], self.Photon))
            # self.OpticalElements[-1].Photon = self.Photon


    def CheckElementsForCrossingPoint(self):
        '''Finds the closest valid element that intersects with the current trajectory'''
        #REFACTOR ME!
        
        self.SmallestDistance = np.inf
        #Calculate Intersection of trajectory with each element. Afterwards only the first intersected element is used for the next step
        n = 0
        self.el_n = 0
        self.Photon.PointCy = [0, 0, 0]
        for el in self.OpticalElements:
            point = el.Geometry.FindIntersections()
            #Only Check stuff if there is an intersection with the specified element
            if np.all(point != [np.nan, np.nan, np.nan]):
                #Check if points are at the smallest distance and not the original point, leaving a small tolerance for round-off errors and the like. could be even smaller than 1e-8 in small systems.
                if np.linalg.norm(self.Photon.Location-point) < self.SmallestDistance and np.linalg.norm(self.Photon.Location-point) > 1e-8:
                    #Check if last position is element
                    if self.Photon.WasElement:
                        #Check if points are at the same side of operation plane depending on the performed operation (Check if photon is supposed to cross the plane or not)
                        if self.OldElement.Type.IsValidNextPoint(self.OldPoint, point):
                            self.SmallestDistance = np.linalg.norm(self.Photon.Location - point)
                            self.FirstEl = el
                            self.Photon.PointCy = point
                    else:
                        #For newly emitted photons not originating from a previous element:
                        #check if ray propagates in the direction of the direction vector, else, do not allow intersection
                        #check which component of direction vector is nonzero
                        nzcomp = np.argwhere(np.array(self.Photon.Direction) != 0)[0][0]
                        if (point[nzcomp] - self.Photon.Location[nzcomp])/self.Photon.Direction[nzcomp] > 0:
                            self.SmallestDistance = np.linalg.norm(self.Photon.Location - point)
                            self.FirstEl = el
                            self.Photon.PointCy = point

            n += 1
        #print self.Photon.PointCy
        self.OldElement = self.FirstEl


    def runSimulation(self):
        '''Simulation Main Loop: Checks Elements for crossing Points, Checks if Photon is absorbed, also executes Elements, reflects, refracts etc.'''
        # self.DrawElements()
        Draw = DrawGraph.DrawElements(self.OpticalElements)
        Draw.DrawElementsPoly()

        self.Photon.GeneratePhoton()
        self.Photon.BeenAbsorbed = [False, 0]

        for c_Ph in xrange(self.NumberOfPhotons):
            #Clear Plot for each new photon, avoid cluttering
            print "#################  " + str(c_Ph+1) + " / " + str(self.NumberOfPhotons) + "  ###################" # CHANGED
            self.OldPoint = self.Photon.Location
            self.Photon.ResetPhoton(self.Elements[0])

            #Only draw 7 photons per plot
            if c_Ph%7 == 0:
                Draw.DrawElementsPoly()

            active = True
            i = 0
            while active:
                #print i
                self.CheckElementsForCrossingPoint()
                self.Photon.CheckForAbsorption(self.OldPoint, self.SmallestDistance)
                self.OldPoint = self.Photon.OldPoint
                #If no absorption just go on as usual
                #if photon is thermally destroyed, absorption
                #If reemitted, just go on in the loop
                if self.Photon.Status == "Start":
                    self.Photon.OperationPlane = self.OldElement.Geometry.OperationPlane()
                    self.OldElement.Type.ExecuteOperation()
                    self.OldPoint = self.Photon.OldPoint

                Draw.DrawLine(self.OldPoint, self.Photon.Location, self.Photon.LineColor)

                if self.Photon.Status == "InternalConverted" or self.Photon.Status == "Dumped" or i>=1e4:
                    active = False
                i += 1

            if self.Photon.BeenAbsorbed[0]:
                self.Photon.BeenAbsorbed[1] += 1
            self.Photon.BeenAbsorbed[0] = False

            gc.collect()


        self.IOOperator.SaveLog()






