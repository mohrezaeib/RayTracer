import numpy as np
import RTUtils as RTU
import copy as cp
import scipy.special as special
import MonteCarloFRET as MCF

class Photon():
    #REFACTOR ME
    def __init__(self, Location=[0, 0, 0], EmissionMode="Directional Isotropic", EmissionDirection=[0, 0, 0], EmissionWavelengthType="Sun spectrum AG15g", StartingWavelength=450.0, StartingMaterial=None, RespectPolarization=True, StartingLocation = [0,0,0], ScaleFactor = 1.0):

        self.Location = Location
        self.EmissionMode = EmissionMode
        self.EmissionDirection = EmissionDirection
        self.EmissionWavelengthType = EmissionWavelengthType
        self.StartingWavelength = StartingWavelength
        self.StartingMaterial = StartingMaterial
        self.RespectPolarization = RespectPolarization
        self.StartingLocation = StartingLocation
        self.ScaleFactor = ScaleFactor

        self.PointCy = [0, 0, 0]
        self.OperationPlane = [[0, 0, 0], [0, 0, 0]]
        self.WasElement = False
        self.BeenAbsorbed = [False, 0]
        self.sgn = 1
        self.Reflected = False
        self.Psi = 0
        self.Direction = RTU.GenerateDirectionSine()
        self.CurrentMaterial = StartingMaterial
        self.Status = "Start"
        self.Wavelength = 450.0
        self.DistanceLeftToAbsorption = 0.0
        self.ElectricFieldVector = [1, 0, 0]

        self.OldPoint = [0,0,0]

        self.LineColor = "black"

    def EmitSun(self):
        '''Emits Photon from the Sun Spectrum'''
        am15 = np.loadtxt("AM15.txt")
        #This converts the W/m^2 spectrum to Photons/s spectrum (only proportional though,  normalized anyway)
        #https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
        am15[:, 1] *= am15[:, 0]
        am15 = RTU.RectifySpectrum(am15)
        return np.random.choice(am15[:, 0], p=am15[:, 1])


    def ResetPhoton(self,OperationPlane):
        self.GeneratePhoton()
        self.OperationPlane = cp.deepcopy(OperationPlane) ##[self.OpticalElements[0].SupportVector,self.OpticalElements[0].DirectionVector]
        self.StartingWavelength = cp.deepcopy(self.Wavelength)
        self.Reflected = False


    def FindSign(self):
        '''Find the sign of the normal vector for reflection,  refraction etc'''
        self.sgn = np.sign(-np.dot(self.Direction, self.OperationPlane[1]))


    def GeneratePhoton(self):
        '''Generate a new Photon
        [self.NumPhEdit.value(),  map(float, self.StartPhEdit.text().split(", ")),  self.EmDPhEdit.currentText(), self.WavelPhEdit.value(), self.StMatPhEdit.currentText()]'''
        self.WasElement = False
        # self.Photon = [np.array([[-0.5, -0.5, -0.1], RTU.GenerateDirectionSine()]), 'Start', 450.0, "MaterialHere", 0, [1, 0, 0], 450.0]
        self.Location = cp.deepcopy(self.StartingLocation)

        if self.EmissionMode == "Directional Isotropic":
            self.Direction = RTU.GenerateDirectionSine()
            self.Direction[2] = np.absolute(self.Direction[2])
            tm = np.roll(RTU.MakeVectAxis(self.EmissionDirection), 2, axis=1)
            self.Direction = np.dot(self.Direction, tm)

            #Generate photons within rectangle
            #2x2cm
            #self.Location = [np.random.rand()*2 - 1.0,  np.random.rand()*2 - 1.0, self.Location[2]]
            #Hexagon
            #self.Location = [np.random.rand()*2 - 1.0, np.random.rand()*(2.0/11.0) - (1.0/11.0), self.Location[2]]
            # self.Photon[0][0] = [np.random.rand()*5 - 2.5,  np.random.rand()*5 - 2.5, self.Photon[0][0][2]]
        elif self.EmissionMode == "Directional Strict":
            self.Direction = cp.deepcopy(self.EmissionDirection)
            #Generate photons within rectangle
            #self.Location = [np.random.rand()*2 - 1.0, np.random.rand()*(2.0/11.0) - (1.0/11.0), self.Location[2]]
            #self.Location = [np.random.rand()*2 - 1.0, np.random.rand()*(2.0) - (1.0), self.Location[2]]

        self.Direction /= np.linalg.norm(self.Direction)

        if self.EmissionWavelengthType == "Sun spectrum AG15g":
            self.Wavelength = self.EmitSun()

        else:
            self.Wavelength = cp.deepcopy(self.StartingWavelength)

        print "Wavelength", self.Wavelength

        self.CurrentMaterial = cp.deepcopy(self.StartingMaterial)
        self.UpdatePolarizationRandom(self.Direction)


    def UpdatePolarizationMolOrientation(self, MolOrient, EmissionDirection):
        '''Update the molecules electric field vector depending on the molecular orientation and the emission direction'''
        v1 = np.cross(MolOrient, EmissionDirection)
        v2 = np.cross(EmissionDirection, v1)
        v2 /= np.linalg.norm(v2)
        self.ElectricFieldVector = v2

    def UpdatePolarizationRandom(self, EmissionDirection):
        '''Find a random vector perpendicular to the emission direction and set it as the new electric field vector'''
        #Find random vector on sphere
        v = np.random.randn(3)
        v /= np.linalg.norm(v)
        #Cross product with  emission direction
        self.ElectricFieldVector = np.cross(EmissionDirection, v)
        self.ElectricFieldVector /= np.linalg.norm(self.ElectricFieldVector)


    def UpdatePolarization(self, OldDirection):
        '''Update the Photon's polarization vector depending on the element'''
        #Plot Stuff
        # newpol = self.Photon[0][0] + 0.2*self.Photon[5]
        # self.ax.plot([self.Photon[0][0][0], newpol[0]], [self.Photon[0][0][1], newpol[1]], [self.Photon[0][0][2], newpol[2]],  color = 'black')

        newdir = self.Direction
        CrossProduct = np.cross(newdir, OldDirection)
        if np.linalg.norm(CrossProduct) != 0:
            acsval = np.dot(newdir, OldDirection)/(np.linalg.norm(newdir)*np.linalg.norm(OldDirection))
            if acsval <= 1.0 and acsval >= -1.0:
                self.ElectricFieldVector = np.dot(RTU.RotationMatrix(CrossProduct, -np.arccos(acsval)), self.ElectricFieldVector)
                self.ElectricFieldVector /= np.linalg.norm(self.ElectricFieldVector)

        #plot more stuff
        # newpol = self.Photon[0][0] + 0.2*self.Photon[5]
        # self.ax.plot([self.Photon[0][0][0], newpol[0]], [self.Photon[0][0][1], newpol[1]], [self.Photon[0][0][2], newpol[2]],  color = 'yellow')

    def PrepareExecution(self):
        self.FindSign()
        self.OldPoint = cp.deepcopy(self.Location)
        self.OldDirection = cp.deepcopy(self.Direction)
        self.LineColor = self.CurrentMaterial.color

    def PostProcessExecution(self):
        self.UpdatePolarization(self.OldDirection)
        self.WasElement = True

    def ReflectOnPlane(self):
        '''Reflect the line IncidentRay on the Plane'''
        try:

            self.Location = (np.dot((self.OperationPlane[0]-self.Location), self.OperationPlane[1])/(np.dot(self.Direction, self.OperationPlane[1])))*self.Direction + self.Location
            self.Direction = 2*np.dot(self.OperationPlane[1], self.Direction)*self.OperationPlane[1] - self.Direction

        except:
            print "Ray and plane possibly do not intersect"
            raise



    def FocusThinLens(self, f):
        '''Assume Plane as an ideal thin lens and refract the IncidentRay accordingly'''
        try:
            #Intersection of parallel ray through the center of the lense and the focal plane

            BaseVect = (np.dot((self.OperationPlane[0]-self.Location), self.OperationPlane[1])/(np.dot(self.Direction, self.OperationPlane[1])))*self.Direction + self.Location
            #Check which focal plane should be used since Ray should always cross the lens
            FPlane1 = [self.OperationPlane[0] + f*self.OperationPlane[1], self.OperationPlane[1]]
            FPlane2 = [self.OperationPlane[0] - f*self.OperationPlane[1], self.OperationPlane[1]]
            #If FPlane1 and Photon are on the same side of the lens,  use FPlane2 to find intersection,  or other way around
                   
            if RTU.CheckSides(self.Location, FPlane1[0], self.OperationPlane) == True:
                FPlane = FPlane2
            else:
                FPlane = FPlane1
                
           
            #This Does an actual lens
            #ParallelRay = [self.OperationPlane[0], self.Direction]
            #this does a thin "cylindrical" lens
            ParallelRay = [[BaseVect[0],BaseVect[1],self.OperationPlane[0][2]], self.Direction]
            
            ImgP = (np.dot((FPlane[0]-ParallelRay[0]), FPlane[1])/(np.dot(ParallelRay[1], FPlane[1])))*ParallelRay[1] + ParallelRay[0]
            # self.ax.plot([ParallelRay[0][0], ImgP[0]], [ParallelRay[0][1], ImgP[1]], [ParallelRay[0][2], ImgP[2]])
            DirectionVector = -(ImgP-BaseVect)

            self.Location = BaseVect
            #print "OLDDIRECTION",self.Direction
            self.Direction = DirectionVector/np.linalg.norm(DirectionVector)
            #print "NEWDIRECTION",self.Direction

        except:
            print "Ray and plane possibly do not intersect"
            raise



    def RefractOnPlane(self, n1, n2):
        '''Uses Snells law to refract the IncidentRay through the interface at Plane. The Refraction indices n1, n2 correspond to the media at each side of the interface '''
        try:
            self.Location = (np.dot((self.OperationPlane[0]-self.Location), self.OperationPlane[1])/(np.dot(self.Direction, self.OperationPlane[1])))*self.Direction + self.Location
            self.Direction = (n1/n2)*self.Direction + ((n1/n2)*(-np.dot(self.Direction, self.OperationPlane[1])) - self.sgn*np.sqrt(1-((n1/n2)**2)*(1-(-np.dot(self.Direction, self.OperationPlane[1]))**2))) * self.OperationPlane[1]

        except:
            print "Ray and plane possibly do not intersect"
            raise


    def CheckForAbsorption(self, OldPoint, SmallestDistance):
        '''Check if current material should absorb parts of the photons,  if so,  it should end the Photon trajectory,  randomly'''

        self.OldPoint = OldPoint
        # self.LineColor = "blue"

        #Somwhere here possibly set oldpoint correlctly

        if self.CurrentMaterial.IsAbsorber:
            if self.DistanceLeftToAbsorption > SmallestDistance:
                self.DistanceLeftToAbsorption -= SmallestDistance
                self.Status = 'Start'
            else:
                print "ABSORBED"
                self.BeenAbsorbed[0] = True
                self.WasElement = False
                #If Photon is absorbed,  move along trajectory for rest of the way
                self.OldPoint = cp.deepcopy(self.Location)
                self.Status = 'Absorbed'

                restdirect_v = self.PointCy - self.Location
                restdirect_v /= np.linalg.norm(restdirect_v)
                newp = (restdirect_v)*self.DistanceLeftToAbsorption + self.Location

                # Parent.ax.plot([self.Location[0], newp[0]], [self.Photon.Location[1], newp[1]], [self.Photon.Location[2], newp[2]],  color ='blue')
                self.Location = newp

                #Calculate probabilities to find out if the photon was absorbed by donor or acceptor

                b_d = self.CurrentMaterial.conc_b * 1e-3 *np.log(10)* self.CurrentMaterial.emax_D * np.maximum(1e-99, self.CurrentMaterial.sp_abs_D[RTU.FindClosestListElement(self.CurrentMaterial.sp_abs_D[:, 0], self.Wavelength), 1])
                b_a = self.CurrentMaterial.conc_a * 1e-3 *np.log(10)* self.CurrentMaterial.emax_A * np.maximum(1e-99, self.CurrentMaterial.sp_abs_A[RTU.FindClosestListElement(self.CurrentMaterial.sp_abs_A[:, 0], self.Wavelength), 1])

                if self.CurrentMaterial.IsAligned:
                    print "Overlap Integral", self.Psi
                    #The 3 comes from the integration over space. For isotropic alignment the PSI integral is 1/3 but it has to be 1,  thus the 3.
                    b_a *= self.Psi * 3

                #if Photon is absorbed by donor (here b) (Not sure if this is correct...)
                if np.random.rand() < (b_d)/(b_d + b_a):
                    #if photon is in donor
                    print "absorbed by donor"
                    fret = MCF.MCFRET(self, True, self.RespectPolarization)
                    self = fret.RunSimulation()
                    #Reset free distance in absorber
                    self.LamBeerDistance()
                else:
                    #If absorbed by acceptor
                    print "absorbed by acceptor"

                    fret = MCF.MCFRET(self, False, self.RespectPolarization)
                    self = fret.RunSimulation()
                    #Reset free distance in absorber
                    self.LamBeerDistance()

                    print self.Wavelength

        else:
            self.Status = 'Start'
            self.WasElement = True


    def LambertBeerFullIntegral(self, n_E):
        '''Calculate the absorption integral for lambert beer and aligned dyes'''
        x = n_E[0]
        y = n_E[1]
        z = n_E[2]
        if self.CurrentMaterial.AlignS > 1e-10:
            ol = np.absolute(1/self.CurrentMaterial.AlignS * (np.exp(-0.5*self.CurrentMaterial.AlignS) * ((-0.2349964007466563j) * (x**2 + y**2 + (2/3.0)* z**2) * special.erf(2.221441469079183/self.CurrentMaterial.AlignS - (0.7071067811865475j) * self.CurrentMaterial.AlignS) + (0.2349964007466563j) * (x**2 + y**2 + (2/3.0) * z**2) * special.erf(2.221441469079183/self.CurrentMaterial.AlignS + (0.7071067811865475j)*self.CurrentMaterial.AlignS) + (0.4699928014933126 * x**2 + 0.4699928014933126 * y**2 + 0.313328534328875 * z**2) * special.erfi(0.7071067811865475 * self.CurrentMaterial.AlignS)) +    np.exp(-4.5 * self.CurrentMaterial.AlignS**2) * ((0.07833213358221876j) * (x**2 + y**2 - 2 * z**2) * special.erf(2.221441469079183/self.CurrentMaterial.AlignS - (2.1213203435596424j) * self.CurrentMaterial.AlignS) - (0.07833213358221876j) * (x**2 +  y**2 - 2* z**2)*special.erf(2.221441469079183/self.CurrentMaterial.AlignS + (2.1213203435596424j)*self.CurrentMaterial.AlignS) + (-0.15666426716443752*x**2 - 0.15666426716443752 * y**2 + 0.31332853432887503 * z**2) * special.erfi(2.1213203435596424 * self.CurrentMaterial.AlignS))))
        else:
            print "Fastlane", n_E
            ol = np.dot(n_E,[0,0,1])**2

        if ol < 1e-10:
            return 0.0
        else:
            return ol

    def LamBeerDistance(self):
        '''Randomly (exponential distribution) assign a length to the photon that it can go within the medium,  only executed when an interface is crossed'''
        #print self.CurrentMaterial.conc_a ,  self.CurrentMaterial.emax_A,  self.CurrentMaterial.conc_b ,  self.CurrentMaterial.emax_D #CHANGED
        if self.CurrentMaterial.IsAbsorber:

            b_a = self.CurrentMaterial.conc_a * 1e-3 * self.CurrentMaterial.emax_A * np.log(10) * np.maximum(1e-99, self.CurrentMaterial.sp_abs_A[RTU.FindClosestListElement(self.CurrentMaterial.sp_abs_A[:, 0], self.Wavelength), 1]) #convert concentration so free length is in cm
            b_d = self.CurrentMaterial.conc_b * 1e-3 * self.CurrentMaterial.emax_D * np.log(10) * np.maximum(1e-99, self.CurrentMaterial.sp_abs_D[RTU.FindClosestListElement(self.CurrentMaterial.sp_abs_D[:, 0], self.Wavelength), 1])

            if self.CurrentMaterial.IsAligned:
                #Numerical Integrator has problems if distribution is too small so the vectors are rotated on the alignment vectors
                #The theta angle can be used to tell the integrator where to look
                #So the align vector is counted as the z-vector [0, 0, 1] and the n_E vector is transformed accordingly
                tm = np.roll(RTU.MakeVectAxis(self.CurrentMaterial.AlignVect), 2, axis=1)
                n_E = np.dot(self.ElectricFieldVector, tm)
                n_E /= np.linalg.norm(n_E)
                #Using analytic expression
                self.Psi = self.LambertBeerFullIntegral(n_E)
                print "Overlap Integral", self.Psi
                #The 3 comes from the integration over space. For isotropic alignment the PSI integral is 1/3 but it has to be 1,  thus the 3.
                b_a *= self.Psi * 3.0
            self.DistanceLeftToAbsorption = (1/self.ScaleFactor) * np.random.exponential(scale=1.0/(b_a + b_d))
        else:
            self.DistanceLeftToAbsorption = 0.0
        #print "NEW FREE LENGTH", self.DistanceLeftToAbsorption #CHANGED










##
