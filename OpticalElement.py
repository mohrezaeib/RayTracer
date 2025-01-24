import copy as cp
import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF
import Geometries
import Elements

class OpticalElement():
    def __init__(self, SupportVector, DirectionVector, ElementProperties, ElementGeometry, ElementSize, Photon):
        '''Initialize elements,  would be nicer to directly parse elements the correct way'''
        self.SupportVector = SupportVector
        self.DirectionVector = DirectionVector/np.linalg.norm(DirectionVector)
        self.ElementGeometry = ElementGeometry
        self.ElementSize = ElementSize
        self.ElementProperties = ElementProperties
        self.Photon = Photon
        self.Photon.Reflected = False

        '''Translate Parsed Data to Readable Variables, This should be changed once the parser works right or simply be made part of the parsing code'''

        shapes = ["None", "Rectangle", "Circle", "Cylinder", "Triangle"]
        self.Shape = shapes[self.ElementGeometry]

        if self.Shape == "Circle":
             self.Geometry = Geometries.Circle(self.ElementSize[0], self.SupportVector, self.DirectionVector,self.Photon)
        elif self.Shape == "Cylinder":
            self.Geometry = Geometries.Cylinder(self.ElementSize[0], self.SupportVector, self.DirectionVector,self.Photon)
        elif self.Shape == "Rectangle":
            if len(self.ElementSize) == 2:
                NV = [0,0,0]
            else:
                NV = self.ElementSize[2]
            self.Geometry = Geometries.Rectangle(self.ElementSize[0], self.ElementSize[1], self.SupportVector, self.DirectionVector,self.Photon,NV)
        elif self.Shape == "Triangle":
            self.Geometry = Geometries.Triangle(self.ElementSize[0],self.ElementSize[1],self.ElementSize[2],self.Photon)


        types = ["None", "Mirror", "Interface", "Lens", "Dump"]
        self.Type = types[self.ElementProperties[0]]
        self.ElementType = Elements.Mirror(self.Photon)
        if self.Type == "Interface":
            self.ElementType = Elements.Interface(self.ElementProperties[1],self.ElementProperties[2], self.Photon)
        elif self.Type == "Lens":
            self.ElementType = Elements.Lens(self.ElementProperties[1], self.Photon)
        elif self.Type == "Dump":
            self.ElementType = Elements.Dump(self.ElementProperties[4], self.Photon)


    def FindIntersections(self):
        '''Find intersection points of Ray with Plane,  possibly a bug if the size of planes is different in both directions,  then needs a check of t1/t2 before return'''
        #Look up if Ray and Plane are parallel...
        if np.dot(self.Photon.Direction, self.DirectionVector) != 0:
            return self.Geometry.FindIntersections()
        else:
            return [np.nan, np.nan, np.nan]


    def FindSign(self, Plane):
        '''Find the sign of the normal vector for reflection,  refraction etc'''
        self.Photon.sgn = np.sign(-np.dot(self.Photon.Direction, Plane[1]))


    def ExecuteCurrentElement(self):
        '''Look at the element where the photon is currently located and act accordingly'''
        self.Photon.OperationPlane = self.Geometry.OperationPlane()

        self.FindSign(self.Photon.OperationPlane)

        self.OldPoint = cp.deepcopy(self.Photon.Location)
        self.OldDirection = cp.deepcopy(self.Photon.Direction)
        self.LineColor = self.Photon.CurrentMaterial.color

        #Check which kind of element was hit and act accordingly
        self.ElementType.ExecuteOperation()

        self.Photon.UpdatePolarization(self.OldDirection)
        self.Photon.WasElement = True



    def CheckForAbsorption(self, OldPoint, SmallestDistance):
        '''Check if current material should absorb parts of the photons,  if so,  it should end the Photon trajectory,  randomly'''

        self.OldPoint = OldPoint
        self.LineColor = "blue"

        if self.Photon.CurrentMaterial.IsAbsorber:
            if self.Photon.DistanceLeftToAbsorption > SmallestDistance:
                self.Photon.DistanceLeftToAbsorption -= SmallestDistance
                self.Photon.Status = 'Start'
            else:
                print "ABSORBED"
                self.Photon.BeenAbsorbed[0] = True
                self.Photon.WasElement = False
                #If Photon is absorbed,  move along trajectory for rest of the way
                self.OldPoint = cp.deepcopy(self.Photon.Location)
                self.Photon.Status = 'Absorbed'

                restdirect_v = self.Photon.PointCy - self.Photon.Location
                restdirect_v /= np.linalg.norm(restdirect_v)
                newp = (restdirect_v)*self.Photon.DistanceLeftToAbsorption + self.Photon.Location

                # Parent.ax.plot([self.Photon.Location[0], newp[0]], [self.Photon.Location[1], newp[1]], [self.Photon.Location[2], newp[2]],  color ='blue')
                self.Photon.Location = newp

                #Calculate probabilities to find out if the photon was absorbed by donor or acceptor

                b_d = self.Photon.CurrentMaterial.conc_b * 1e-3 *np.log(10)* self.Photon.CurrentMaterial.emax_D * np.maximum(1e-99, self.Photon.CurrentMaterial.sp_abs_D[RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_abs_D[:, 0], self.Photon.Wavelength), 1])
                b_a = self.Photon.CurrentMaterial.conc_a * 1e-3 *np.log(10)* self.Photon.CurrentMaterial.emax_A * np.maximum(1e-99, self.Photon.CurrentMaterial.sp_abs_A[RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_abs_A[:, 0], self.Photon.Wavelength), 1])

                if self.Photon.CurrentMaterial.IsAligned:
                    print "Overlap Integral", self.Photon.Psi
                    #The 3 comes from the integration over space. For isotropic alignment the PSI integral is 1/3 but it has to be 1,  thus the 3.
                    b_a *= self.Photon.Psi * 3

                #if Photon is absorbed by donor (here b) (Not sure if this is correct...)
                if np.random.rand() < (b_d)/(b_d + b_a):
                    #if photon is in donor
                    print "absorbed by donor"
                    fret = MCF.MCFRET(self.Photon, True, self.Photon.RespectPolarization)
                    self.Photon = fret.RunSimulation()
                    #Reset free distance in absorber
                    self.LamBeerDistance()
                else:
                    #If absorbed by acceptor
                    print "absorbed by acceptor"

                    fret = MCF.MCFRET(self.Photon, False, self.Photon.RespectPolarization)
                    self.Photon = fret.RunSimulation()
                    #Reset free distance in absorber
                    self.LamBeerDistance()

                    print self.Photon.Wavelength

        else:
            self.Photon.Status = 'Start'
            self.Photon.WasElement = True



    def LambertBeerFullIntegral(self, n_E):
        '''Calculate the absorption integral for lambert beer and aligned dyes, integral has been solved and is here computed with hard coded numerical values.'''
        x = n_E[0]
        y = n_E[1]
        z = n_E[2]
        ol = np.absolute(1/self.Photon.CurrentMaterial.AlignS * (np.exp(-0.5*self.Photon.CurrentMaterial.AlignS) * ((-0.2349964007466563j) * (x**2 + y**2 + (2/3.0)* z**2) * special.erf(2.221441469079183/self.Photon.CurrentMaterial.AlignS - (0.7071067811865475j) * self.Photon.CurrentMaterial.AlignS) + (0.2349964007466563j) * (x**2 + y**2 + (2/3.0) * z**2) * special.erf(2.221441469079183/self.Photon.CurrentMaterial.AlignS + (0.7071067811865475j)*self.Photon.CurrentMaterial.AlignS) + (0.4699928014933126 * x**2 + 0.4699928014933126 * y**2 + 0.313328534328875 * z**2) * special.erfi(0.7071067811865475 * self.Photon.CurrentMaterial.AlignS)) +    np.exp(-4.5 * self.Photon.CurrentMaterial.AlignS**2) * ((0.07833213358221876j) * (x**2 + y**2 - 2 * z**2) * special.erf(2.221441469079183/self.Photon.CurrentMaterial.AlignS - (2.1213203435596424j) * self.Photon.CurrentMaterial.AlignS) - (0.07833213358221876j) * (x**2 +  y**2 - 2* z**2)*special.erf(2.221441469079183/self.Photon.CurrentMaterial.AlignS + (2.1213203435596424j)*self.Photon.CurrentMaterial.AlignS) + (-0.15666426716443752*x**2 - 0.15666426716443752 * y**2 + 0.31332853432887503 * z**2) * special.erfi(2.1213203435596424 * self.Photon.CurrentMaterial.AlignS))))

        '''Leave some tolerance here as well'''
        if ol < 1e-10:
            return 0.0
        else:
            return ol

    def LamBeerDistance(self):
        '''Randomly (exponential distribution) assign a length to the photon that it can go within the medium,  only executed when an interface is crossed'''
        if self.Photon.CurrentMaterial.IsAbsorber:

            b_a = self.Photon.CurrentMaterial.conc_a * 1e-3 * self.Photon.CurrentMaterial.emax_A * np.log(10) * np.maximum(1e-99, self.Photon.CurrentMaterial.sp_abs_A[RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_abs_A[:, 0], self.Photon.Wavelength), 1]) #convert concentration so free length is in cm
            b_d = self.Photon.CurrentMaterial.conc_b * 1e-3 * self.Photon.CurrentMaterial.emax_D * np.log(10) * np.maximum(1e-99, self.Photon.CurrentMaterial.sp_abs_D[RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_abs_D[:, 0], self.Photon.Wavelength), 1])

            if self.Photon.CurrentMaterial.IsAligned:
                #Numerical Integrator has problems if distribution is too small so the vectors are rotated on the alignment vectors
                #The theta angle can be used to tell the integrator where to look
                #So the align vector is counted as the z-vector [0, 0, 1] and the n_E vector is transformed accordingly
                tm = np.roll(RTU.MakeVectAxis(self.Photon.CurrentMaterial.AlignVect), 2, axis=1)
                n_E = np.dot(self.Photon.ElectricFieldVector, tm)
                n_E /= np.linalg.norm(n_E)
                #Using analytic expression
                self.Photon.Psi = self.LambertBeerFullIntegral(n_E)
                print "Overlap Integral", self.Photon.Psi
                #The 3 comes from the integration over space. For isotropic alignment the PSI integral is 1/3 but it has to be 1,  thus the 3.
                b_a *= self.Photon.Psi * 3.0
            self.Photon.DistanceLeftToAbsorption = (1/self.Photon.ScaleFactor) * np.random.exponential(scale=1.0/(b_a + b_d))
        else:
            self.Photon.DistanceLeftToAbsorption = 0.0
        print "NEW FREE LENGTH", self.Photon.DistanceLeftToAbsorption






#
