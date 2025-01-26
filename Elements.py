import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF


class Mirror():
    '''Class defining behavior of mirrors'''
    def __init__(self, Photon):
        self.Photon = Photon
        self.color = "blue"
        self.alpha = 0.05
        self.BounceCounter = 0

    def ExecuteOperation(self):
        '''
        Do what a mirror does, in this case, reflect and confirm operation
        '''
        self.Photon.PrepareExecution()
        self.Photon.ReflectOnPlane()
        self.BounceCounter += 1
        self.Photon.PostProcessExecution()

    def IsValidNextPoint(self, OldPoint, Point):
        '''
        Check if new trajectory goes towards the correct side.
        In case of a mirror the photon has to stay on the same side.
        '''
        return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) is True


class Interface():
    '''Class defining behavior of interfaces'''
    def __init__(self, Material1, Material2, Photon):
        self.MaterialSide1 = Material1
        self.MaterialSide2 = Material2
        self.Photon = Photon
        self.color = "yellow"
        self.alpha = 0.1
        self.BounceCounter = 0

    def ExecuteOperation(self):
        '''
        Take care of interfaces: randomly reflect or refract or totally reflect
        '''
        self.Photon.PrepareExecution()

        # Find out in which material the photon is traveling, then assign refractive properties
        if self.Photon.CurrentMaterial.name == self.MaterialSide1.name:
            mat = self.MaterialSide2
        else:
            mat = self.MaterialSide1

        # Check for total reflection or refraction
        if self.IsTotalReflect(self.Photon.CurrentMaterial.n, mat.n):
            # Total reflection
            self.Photon.ReflectOnPlane()
            self.Photon.Reflected = True
            self.BounceCounter += 1
        else:
            # Fresnel-based reflection vs. refraction
            if np.random.rand() > self.Reflectance(self.Photon.CurrentMaterial.n, mat.n):
                # Crossing of interface
                self.Photon.RefractOnPlane(self.Photon.CurrentMaterial.n, mat.n)
                # Also change wavelength when interface is crossed; code not necessary:
                # self.Photon.Wavelength = self.Photon.Wavelength * (self.Photon.CurrentMaterial.n / mat.n)
                # #print(self.Photon.CurrentMaterial.name, mat.name)  #CHANGED
                self.Photon.CurrentMaterial = mat
                self.Photon.Reflected = False
                # If the interface is crossed: assign a new free length for the photon to travel in the new medium
                self.Photon.LamBeerDistance()
            else:
                # If reflected
                self.Photon.ReflectOnPlane()
                self.Photon.Reflected = True
                self.BounceCounter += 1

        self.Photon.PostProcessExecution()

    def IsTotalReflect(self, n1, n2):
        '''
        Checks if total reflection happens
        '''
        thetai = np.arccos(
            np.dot(
                self.Photon.Direction,
                -self.Photon.sgn * self.Photon.OperationPlane[1]
            ) / (
                np.linalg.norm(self.Photon.Direction) *
                np.linalg.norm(self.Photon.OperationPlane[1])
            )
        )
        r = False
        if n1 > n2 and np.abs(thetai) > np.arcsin(n2 / n1):
            r = True
        return r

    def Reflectance(self, n1, n2):
        '''
        Calculate the Reflectance at an interface between n1, n2 using Fresnel Equations
        '''
        thetai = np.arccos(
            -np.dot(
                self.Photon.Direction,
                self.Photon.sgn * self.Photon.OperationPlane[1]
            )
        )
        thetat = np.arcsin((n1 / n2) * np.sin(thetai))
        R_perp = (
            (n1 * np.cos(thetai) - n2 * np.cos(thetat)) /
            (n1 * np.cos(thetai) + n2 * np.cos(thetat))
        ) ** 2
        R_par = (
            (n2 * np.cos(thetai) - n1 * np.cos(thetat)) /
            (n2 * np.cos(thetai) + n1 * np.cos(thetat))
        ) ** 2

        # #print("Reflectance", R)
        '''
        Calculate the incident plane and the projection of the electric field vector
        onto the incident plane, then weigh components
        '''
        ElectricVect = self.Photon.ElectricFieldVector
        ElectricVect /= np.linalg.norm(ElectricVect)
        n_EE = np.cross(
            self.Photon.Direction,
            self.Photon.sgn * self.Photon.OperationPlane[1]
        )
        nrm = np.linalg.norm(n_EE)

        if nrm != 0:
            n_EE /= nrm
            dp = np.dot(ElectricVect, n_EE)
            sproj_v = dp * n_EE
            parproj_v = ElectricVect - sproj_v
            parproj = np.linalg.norm(parproj_v)
            sproj = np.abs(dp)
            R = (parproj ** 2) * R_par + (sproj ** 2) * R_perp
        else:
            R = ((n1 - n2) / (n1 + n2)) ** 2

        return R

    def IsValidNextPoint(self, OldPoint, Point):
        '''
        Check if the next point is valid based on whether the photon is reflected or not
        '''
        return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) == self.Photon.Reflected


class Lens():
    '''Class defining behavior of lenses'''
    def __init__(self, FocalLength, Photon):
        self.FocalLength = FocalLength
        self.Photon = Photon
        self.color = "red"
        self.alpha = 0.1

    def ExecuteOperation(self):
        '''
        Do what a thin lens does, in this case: focus
        '''
        self.Photon.PrepareExecution()
        self.Photon.FocusThinLens(self.FocalLength)
        self.Photon.PostProcessExecution()

    def IsValidNextPoint(self, OldPoint, Point):
        '''
        For a lens, check if we are going to the correct side
        '''
        return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) is False


class Dump():
    '''Class defining behavior of detectors'''
    def __init__(self, DetectionSpectrum, Photon):
        self.Photon = Photon
        self.Counter = 0
        self.EnergyCounter = 0
        self.SpectralCounter = 0

        self.IncidentSpectrum = np.zeros((3800, 2))
        self.IncidentSpectrum[:, 0] = np.arange(200, 4000)

        self.DetectionSpectrum = DetectionSpectrum

        self.color = "black"
        self.alpha = 0.2

    def ExecuteOperation(self):
        '''
        Stop photon and log its wavelength to the respective detector
        '''
        self.Photon.PrepareExecution()

        # Mark photon as dumped
        self.Photon.Status = 'Dumped'
        self.Counter += 1
        self.EnergyCounter += 1 * self.Photon.StartingWavelength / self.Photon.Wavelength

        # Increment the photon count for the detected wavelength
        index = RTU.FindClosestListElement(
            self.IncidentSpectrum[:, 0],
            self.Photon.Wavelength
        )
        self.IncidentSpectrum[index, 1] += 1

        # Decide if the photon is counted in the detection spectrum
        if isinstance(self.DetectionSpectrum, int):
            self.SpectralCounter += 1
        else:
            detection_index = RTU.FindClosestListElement(
                self.DetectionSpectrum[:, 0],
                self.Photon.Wavelength
            )
            if self.DetectionSpectrum[detection_index, 1] > np.random.rand():
                self.SpectralCounter += 1

        # Move photon location to end
        self.Photon.Location = self.Photon.PointCy

        self.Photon.PostProcessExecution()
