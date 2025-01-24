import numpy as np
import cPickle
import RayTrace
import Elements
import datetime
import os
import Materials
import sys
import PhotonObject as PhObj
import __builtin__


class LoadSaveRT():
    def __init__(self,Path,RayTracerInstance):
        self.RTI = RayTracerInstance
        self.ActivePath = Path

    
    def LoadArrays(self):
        f = open(self.ActivePath + "/Elements.prt", 'rb')
        self.Elements = np.array(cPickle.load(f))
        f = open(self.ActivePath + "/ElementProperties.prt", 'rb')
        self.ElementProperties = cPickle.load(f)
        f = open(self.ActivePath + "/ElementGeometry.prt", 'rb')
        self.ElementGeometry = cPickle.load(f)
        f = open(self.ActivePath + "/ElementSize.prt", 'rb')
        self.ElementSize = cPickle.load(f)
        f = open(self.ActivePath + "/Materials.prt", 'rb')
        self.Materials = cPickle.load(f)
        f = open(self.ActivePath + "/Photon.prt", 'rb')
        self.PhArray = cPickle.load(f)


    def SaveLog(self):
        '''Save the jobfiles as a copy, also save results'''
        foldername = self.ActivePath + "/" + datetime.datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss') + "/"

        if not os.path.exists(foldername):
            os.makedirs(foldername)


        f = open(foldername + "/Elements.prt", 'wb')
        cPickle.dump(self.RTI.Elements, f)
        f = open(foldername + "/ElementProperties.prt", 'wb')
        cPickle.dump(self.RTI.ElementProperties, f)
        f = open(foldername + "/ElementGeometry.prt", 'wb')
        cPickle.dump(self.RTI.ElementGeometry, f)
        f = open(foldername + "/ElementSize.prt", 'wb')
        cPickle.dump(self.RTI.ElementSize, f)
        f = open(foldername + "/Materials.prt", 'wb')
        cPickle.dump(self.Materials, f)
        f = open(foldername + "/Photon.prt", 'wb')
        cPickle.dump(self.RTI.PhArray, f)

        #save logfile only
        loglist = np.zeros((len(self.RTI.OpticalElements), 4))
        for i in range(len(self.RTI.OpticalElements)):
            if isinstance(self.RTI.OpticalElements[i].Type,Elements.Dump) :
                loglist[i, 0] = self.RTI.OpticalElements[i].Type.Counter
                loglist[i, 1] = self.RTI.OpticalElements[i].Type.EnergyCounter
                loglist[i, 2] = self.RTI.OpticalElements[i].Type.SpectralCounter
                np.savetxt(foldername + "/DetectorSpectrum" + str(i) + ".txt", self.RTI.OpticalElements[i].Type.IncidentSpectrum)

            elif isinstance(self.RTI.OpticalElements[i].Type,Elements.Mirror) or isinstance(self.RTI.OpticalElements[i].Type,Elements.Interface):
                loglist[i, 3] = self.RTI.OpticalElements[i].Type.BounceCounter

        np.savetxt(foldername + "/Logfile.txt", loglist)

        np.savetxt(foldername + "/AbsorbedPhotons.txt", np.array([self.RTI.Photon.BeenAbsorbed[1]]))


    
    def ParsePath(self, Path):
        '''Load Data from saved path'''
        
        '''Unpickle the data from the Gui'''
        self.ActivePath = Path
        self.LoadArrays()

        for m in self.Materials:
            self.RTI.MaterialDict.setdefault(m[0],Materials.MaterialObj( name = m[0], n = m[1]  ,IsAbsorber = m[2], conc_a = m[4], conc_b = m[3]  , color =   m[5]  , IsAligned =  m[6] , AlignVect =   m[7] , AlignS =  m[8] , R_0=   m[9] , Q_D =   m[10]  , Q_A =   m[11]  , tau_D =  m[12]  , tau_A =   m[13]  , emax_D =   m[14]  , emax_A =   m[15]  , sp_abs_D =  m[16]  , sp_em_D =  m[17]  , sp_abs_A =   m[18]  , sp_em_A =    m[19]  , AbsEmAngle =   m[20]  , FixDist =   m[21]  , FixDistR =   m[22]  , FixDistNeighbors =   m[23]   ))


        for i in range(len(self.ElementGeometry)):
            if self.ElementProperties[i][0] == 2:
                self.ElementProperties[i] = [ 2, self.RTI.MaterialDict[self.ElementProperties[i][1]] , self.RTI.MaterialDict[self.ElementProperties[i][2]]]
            elif self.ElementProperties[i][0] == 4:
                if self.ElementProperties[i][3] == '':
                    self.ElementProperties[i] = [4, 0, 0, 0, 1]
                else:
                    self.ElementProperties[i] = [4, 0, 0, 0, np.loadtxt(self.ElementProperties[i][3])]

        self.RTI.Elements = self.Elements
        self.RTI.ElementGeometry = self.ElementGeometry
        self.RTI.ElementSize = self.ElementSize
        self.RTI.ElementProperties =  self.ElementProperties

        '''Create array of objects from parsed data arrays'''

        self.RTI.NumberOfPhotons = self.PhArray[0]
        print self.PhArray
        StartingMaterial = self.RTI.MaterialDict[self.PhArray[5]]

        if len(self.PhArray) == 9:
            ScaleFactor = self.PhArray[8]
        else:
            ScaleFactor = 1.0

        self.RTI.Photon = PhObj.Photon(Location=self.PhArray[1], EmissionMode=self.PhArray[2], EmissionDirection=self.PhArray[6], EmissionWavelengthType=self.PhArray[3], StartingWavelength=self.PhArray[4], StartingMaterial = StartingMaterial, RespectPolarization=self.PhArray[7], StartingLocation=self.PhArray[1], ScaleFactor = ScaleFactor)

        self.RTI.lineclr = self.RTI.Photon.StartingMaterial.color


        self.RTI.CreateOpticalElementList()
        self.LoadArrays()

