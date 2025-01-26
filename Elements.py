

import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF


class Mirror():
    '''Class defining behavior of mirrors'''
    def __init__(self,Photon):        
        self.Photon = Photon
        self.color = "blue"
        self.alpha = 0.05
        self.BounceCounter = 0


    def ExecuteOperation(self):
        '''Do what a mirror does, in this case, reflect and confirm operation'''
        self.Photon.PrepareExecution()
        self.Photon.ReflectOnPlane()
        self.BounceCounter += 1
        self.Photon.PostProcessExecution()


    def IsValidNextPoint(self,OldPoint, Point):
        '''Check if new trajectory goes towards the correct side. In case of a mirror the photon has to stay on the same side'''
        return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) == True


class Interface():
    '''Class defining behavior of interfaces'''
    def __init__(self,Material1,Material2, Photon):
        self.MaterialSide1 = Material1
        self.MaterialSide2 = Material2
        self.Photon = Photon
        self.color = "yellow"
        self.alpha = 0.1
        self.BounceCounter = 0

    def ExecuteOperation(self):
        '''take care of interfaces, randomly reflect or refract or total reflect'''
        #find out in which Material the photon is traveling and assign refractive properties accordingly
        self.Photon.PrepareExecution()        
        if self.Photon.CurrentMaterial.name == self.MaterialSide1.name:
            mat = self.MaterialSide2
        else:
            mat = self.MaterialSide1

        #Check For total reflection or refraction
        if self.IsTotalReflect(self.Photon.CurrentMaterial.n, mat.n):
            #TotalReflection
            self.Photon.ReflectOnPlane()
            self.Photon.Reflected = True
            self.BounceCounter += 1
        else:
            if np.random.rand() > self.Reflectance(self.Photon.CurrentMaterial.n, mat.n):
                #Crossing of interface
                self.Photon.RefractOnPlane(self.Photon.CurrentMaterial.n, mat.n)
                #Also change wavelength when interface is crossed,  this code is not necessary
                # self.Photon.Wavelength = self.Photon.Wavelength * (self.Photon.CurrentMaterial.n/mat.n)
                #print self.Photon.CurrentMaterial.name,mat.name #CHANGED
                self.Photon.CurrentMaterial = mat
                self.Photon.Reflected = False
                #If the interface is crossed: Assign a new free length for the photon to travel within the medium,  can be done regarless of medium since absorption is only checked if medium is defined as absorber
                self.Photon.LamBeerDistance()
            else:
                #If reflected
                self.Photon.ReflectOnPlane()
                self.Photon.Reflected = True
                self.BounceCounter += 1

        self.Photon.PostProcessExecution()

    def IsTotalReflect(self, n1, n2):
        '''Checks if TotalReflection happens'''
        thetai = np.arccos(np.dot(self.Photon.Direction, -self.Photon.sgn*self.Photon.OperationPlane[1])/(np.linalg.norm(self.Photon.Direction)*np.linalg.norm(self.Photon.OperationPlane[1])))
        r = False
        if n1 > n2 and np.absolute(thetai) > np.arcsin(n2/n1):
            r = True
        return r


    def Reflectance(self, n1, n2):
        '''Calculate the Reflectance at a interface between n1 n2 using Fresnel Equations'''

        thetai = np.arccos(-np.dot(self.Photon.Direction, self.Photon.sgn*self.Photon.OperationPlane[1]))
        thetat = np.arcsin((n1/n2) * np.sin(thetai))
        R_perp = ((n1*np.cos(thetai)- n2*np.cos(thetat))/(n1*np.cos(thetai) + n2*np.cos(thetat)))**2
        R_par = ((n2*np.cos(thetai)- n1*np.cos(thetat))/(n2*np.cos(thetai) + n1*np.cos(thetat)))**2
        #R = (R_perp+R_par)/2.0
        '''Calculate the incident plane and the projection of the electric field vector onto the incident plane,  then weigh components'''
        ElectricVect = self.Photon.ElectricFieldVector
        ElectricVect /= np.linalg.norm(ElectricVect)
        n_EE = np.cross(self.Photon.Direction, self.Photon.sgn*self.Photon.OperationPlane[1])
        nrm = np.linalg.norm(n_EE)

        if nrm != 0:
            n_EE /= np.linalg.norm(n_EE)
            dp = (np.dot(ElectricVect, n_EE))
            sproj_v = dp * n_EE
            parproj_v = ElectricVect - sproj_v
            parproj = np.linalg.norm(parproj_v)
            sproj = np.absolute(dp)
            R = (parproj**2)*R_par + (sproj**2)*R_perp
        else:
            R = ((n1-n2)/(n1+n2))**2
        #print "Reflectance", R
        return R



    def IsValidNextPoint(self,OldPoint, Point):
         return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) == self.Photon.Reflected


class Lens():
    '''Class defining behavior of lenses'''
    def __init__(self,FocalLength, Photon):
        self.FocalLength = FocalLength
        self.Photon = Photon
        self.color = "red"
        self.alpha = 0.1

    def ExecuteOperation(self):
        '''Do what a thin lens does, in this case: focus'''
        self.Photon.PrepareExecution()
        self.Photon.FocusThinLens(self.FocalLength)
        self.Photon.PostProcessExecution()

    def IsValidNextPoint(self,OldPoint, Point):         
         return RTU.CheckSides(OldPoint, Point, self.Photon.OperationPlane) == False


class Dump():
    '''Class defining behavior of detectors'''
    def __init__(self, DetectionSpectrum, Photon):
        self.Photon = Photon
        self.Counter = 0
        self.EnergyCounter = 0
        self.SpectralCounter = 0

        self.IncidentSpectrum = np.zeros((3800,2))
        self.IncidentSpectrum[:,0] = np.arange(200,4000)

        self.DetectionSpectrum = DetectionSpectrum


        self.color = "black"
        self.alpha = 0.2

    def ExecuteOperation(self):
        '''Stop photon and log its wavelength to the respective detector'''
        self.Photon.PrepareExecution()
        #Check if photon comes from active side of dump
        self.Photon.Status = 'Dumped'
        self.Counter += 1
        self.EnergyCounter += 1 * self.Photon.StartingWavelength / self.Photon.Wavelength
        self.IncidentSpectrum[RTU.FindClosestListElement(self.IncidentSpectrum[:, 0], self.Photon.Wavelength), 1] += 1

        if type(self.DetectionSpectrum) == int:
            self.SpectralCounter += 1
        else:
            if  self.DetectionSpectrum[RTU.FindClosestListElement(self.DetectionSpectrum[:, 0], self.Photon.Wavelength), 1] > np.random.rand():
                self.SpectralCounter += 1

        self.Photon.Location = self.Photon.PointCy

        self.Photon.PostProcessExecution()






#
