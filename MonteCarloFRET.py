import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import copy as cp
# from RayTrace import Material
import RTUtils as RTU

'''Monte Carlo Simulation for FRET dynamics. A single excited molecule is analyzed and the energy transfer probabilities are calculated'''

class MCFRET(object):
    def __init__(self,Photon,InDonor = True, EmPol = True):

        #self.R_0 = Photon.CurrentMaterial.R_0 #FoersterRadius in nm

        self.c_D = Photon.CurrentMaterial.conc_b #Donor concentration in mmol/L
        self.c_A = Photon.CurrentMaterial.conc_a #Acceptor concentration in mmol/L
        self.Q_D = Photon.CurrentMaterial.Q_D #Donor fluorescence quantum yield
        self.Q_A = Photon.CurrentMaterial.Q_A #Acceptor fluorescence quantum yield
        self.tau_A = Photon.CurrentMaterial.tau_A #Donor fluorescence lifetime in ns
        self.tau_D = Photon.CurrentMaterial.tau_D #Donor Internal conversion lifetime in ns

        self.AlignS = Photon.CurrentMaterial.AlignS
        self.AlignVect = Photon.CurrentMaterial.AlignVect
        self.IsAligned = Photon.CurrentMaterial.IsAligned

        self.EmPol = EmPol

        self.k_D = self.Q_D/self.tau_D
        self.k_A = self.Q_A/self.tau_A
        self.k_D_IC = (1-self.Q_D)/self.tau_D
        self.k_A_IC = (1-self.Q_A)/self.tau_A

        self.Photon = Photon
        self.InDonor = InDonor


        self.overlDA = self.Photon.CurrentMaterial.ol_DA
        self.overlAD = self.Photon.CurrentMaterial.ol_AD

        #self.overlDA = self.CalcOverlap(self.Photon.CurrentMaterial.sp_em_D,self.Photon.CurrentMaterial.sp_abs_A,self.Photon.CurrentMaterial.emax_A)
        #self.overlAD = self.CalcOverlap(self.Photon.CurrentMaterial.sp_em_A,self.Photon.CurrentMaterial.sp_abs_D,self.Photon.CurrentMaterial.emax_D)
        #self.TESTEMISSION()

        # '''Plot Stuff'''
        # self.fig = plt.figure()
        # self.ax = self.fig.add_subplot(111, projection='3d')
        # self.ax.set_xlabel('x')
        # self.ax.set_ylabel('y')
        # self.ax.set_zlabel('z')
        # self.ax.set_xlim([-1,1])
        # self.ax.set_ylim([-1,1])
        # self.ax.set_zlim([-1,1])

    def CalcOverlap(self,EmD,AbsA,eA):
        overl = 0
        for l in range(int(np.amin(EmD[:,0])),int(np.amax(AbsA[:,0]))):
            '''Integration here is correct as it is integrated over the integers l with equal interval width'''
            overl += EmD[RTU.FindClosestListElement(EmD[:,0],l),1] * (eA * AbsA[RTU.FindClosestListElement(AbsA[:,0],l),1]) * ((l)**4)
        return overl

    def RunSimulation(self):
        if self.Photon.CurrentMaterial.conc_a > 0:
            self.GenerateMolecules()
            #Calculate the probabilities for each event, first being fluorescence, second internal converison, afterwards energy transfers to acceptors
            p = np.zeros(3)
            p[0] = (1-self.P_ET_DA)*self.Q_D
            p[1] = (1-self.P_ET_DA)*(1-self.Q_D)
            p[2] = self.P_ET_DA

            # print p

            pA = np.zeros(3)
            pA[0] = (1-self.P_ET_AD)*self.Q_A
            pA[1] = (1-self.P_ET_AD)*(1-self.Q_A)
            pA[2] = self.P_ET_AD

            print "probabilities Donor: Fl:" + str(p[0]) + "IC:"+str(p[1]) + "FRET" + str(np.sum(p[2:])) +"Overlap Integral:" + str(self.overlDA)
            print "probabilities Acceptor: Fl:" + str(pA[0]) + "IC:"+str(pA[1]) + "FRET" + str(np.sum(pA[2:])) +"Overlap Integral:" + str(self.overlAD)
        else:
            p = np.array([self.Photon.CurrentMaterial.Q_D,1-self.Photon.CurrentMaterial.Q_D])



        #Now select one Pathway depending on the probabilities
        #Pathway ==0: Fluorescence
        #Pathway ==1: IC
        #Pathway ==2: FRET
        Pathway = np.random.choice(p.shape[0],1,p = p)[0]

        if Pathway == 0 and self.InDonor:
            #fluorescence
            #generate new direction of photon, isotropic since donors are not aligned
            self.Photon.Direction = RTU.GenerateDirectionSine()
            #print self.Photon.Direction # CHANGED

            #Set New Wavelength on Photon
            self.Photon.Wavelength = self.EmitWavelength(True)

            self.Photon.Status = 'ReemittedDonor'
            print 'ReemittedDonor'
            self.Photon.UpdatePolarizationRandom(self.Photon.Direction)

            print self.Photon.Status
        elif Pathway == 1 and self.InDonor:
            #Internal conversion, all is lost
            self.Photon.Status = 'InternalConverted'
            print 'InternalConverted Donor'
        elif Pathway > 1 or not self.InDonor:
            if Pathway > 1 and self.InDonor:
                print "FRET from Donor to Acceptor"
            else:
                print "In Acceptor"


            #Energy transfer to acceptor
            #Choose Pathway in acceptor
            PathwayAcceptor = np.random.choice(pA.shape[0],1,p = pA)[0]
            # PathwayAcceptor = np.random.choice(2,1,p=[self.Q_A,1-self.Q_A])[0]
            if PathwayAcceptor == 0:
                #fluorescence
                print "Reemitted Acceptor"
                if self.IsAligned:
                    #If the acceptor is algined: emit with a sinesquared function
                    self.Photon.Direction = RTU.GenerateDirectionSineSq()
                    #Make Molecule transition dipole moment axis the new z axis so the emission occurs around the equator
                    #tm = np.roll(self.MakeVectAxis(self.AlignVect),2,axis=1)
                    # if self.InDonor:
                    #     remdir = self.GenerateDirectionSine() ##self.direction_A[Pathway - 2]
                    # else:
                    '''reemission direction for aligned molecule'''
                    remdir = RTU.GenerateDirectionGauss(self.Photon)

                    '''If there is an angle between absorption and emission dipole moment, add that angle here, else make the alignment axis the 'new z axis' '''
                    #Create Matrix to make gaussian orientation new z axis
                    tm = np.roll(RTU.MakeVectAxis(remdir),2,axis=1)

                    if self.Photon.CurrentMaterial.AbsEmAngle > 0:
                        th = (self.Photon.CurrentMaterial.AbsEmAngle/180.0) * np.pi
                        ph = np.random.rand() * 2*np.pi

                        v = [0,0,0]
                        v[0] = np.sin(th)*np.cos(ph)
                        v[1] = np.sin(th)*np.sin(ph)
                        v[2] = np.cos(th)

                        etm = np.roll(RTU.MakeVectAxis(v),2,axis=1)
                        self.Photon.Direction = np.dot(etm,self.Photon.Direction)

                    self.Photon.Direction = np.dot(tm,self.Photon.Direction)
                    self.Photon.Direction /= np.linalg.norm(self.Photon.Direction)

                    if self.EmPol:
                        self.Photon.UpdatePolarizationMolOrientation(remdir,self.Photon.Direction)
                    else:
                        self.Photon.UpdatePolarizationRandom(self.Photon.Direction)

                    #Set New Wavelength on Photon
                    self.Photon.Wavelength = self.EmitWavelength(False)

                    self.Photon.Status = 'ReemittedAcceptor'

                    '''Plot molecule direction and emisison direction'''
                    # self.ax.plot([0,self.Photon.Direction[0]],[0,self.Photon.Direction[1]],[0,self.Photon.Direction[2]])
                    # self.ax.plot([0,self.direction_A[Pathway - 2][0]],[0,self.direction_A[Pathway - 2][1]],[0,self.direction_A[Pathway - 2][2]])
                    # self.ax.plot([0,self.Photon[5][0]],[0,self.Photon[5][1]],[0,self.Photon[5][2]])


                else:
                    #if acceptor is not aligned: emit isotropically
                    self.Photon.Direction = RTU.GenerateDirectionSine()
                    self.Photon.UpdatePolarizationRandom(self.Photon.Direction)
                    self.Photon.Wavelength = self.EmitWavelength(False)
                    self.Photon.Status = 'ReemittedAcceptor'

                    '''Plot molecule direction and emisison direction'''
                    # self.ax.plot([0,self.Photon.Direction[0]],[0,self.Photon.Direction[1]],[0,self.Photon.Direction[2]])
                    # self.ax.plot([0,self.Photon[5][0]],[0,self.Photon[5][1]],[0,self.Photon[5][2]])

                print self.Photon.Wavelength
            elif PathwayAcceptor == 1:
                #IC
                print "Internal Converted Acceptor"
                self.Photon.Status = 'InternalConverted'

            else:
                print "FRET from Acceptor to Donor"
                self.InDonor = True
                return self.RunSimulation()

        print "Wavelength:",self.Photon.Wavelength
        return self.Photon

    def GenerateMolecules(self):
        '''Fret Taken from paper Th. Foerster 1959'''
        '''Here FRET from Donor to Acceptor, Calculate FRET probabilities if there is a fixed distance between donor and acceptor or not'''
        if self.Photon.CurrentMaterial.FixDistR == 0:
            raise ValueError("FixDistR is zero, which will cause a division by zero error.")
            invR6 = 1/((self.Photon.CurrentMaterial.FixDistR*1e-9)**6)
            k_ET = ((9*(2/3.0)*self.Q_D*np.log(10.0))/(128*(self.Photon.CurrentMaterial.n**4) * 6.022e23 * (np.pi**5) * self.tau_D) ) * invR6 * (1e-54) * self.overlDA * 1e17
            self.P_ET_DA = k_ET/(1.0/self.tau_D + k_ET)
        else:
            #Fret Radius
            #1e17 comes from: 1 (L/(mol*cm))*nm^4 in (nm^6/mol)
            R_0 = ((9*(2/3.0)*self.Q_D*np.log(10.0))/(128*(self.Photon.CurrentMaterial.n**4) * 6.022e23 * (np.pi**5))) * self.overlDA * 1e17
            R_0 = R_0**(1/6.0)
            print "R_0_DA:",R_0
            #Convert to cm
            R_0 *= 1e-7
            #concentration in molecules/cm^3
            c_A_mpcm = self.c_A * 6.022e17
            #critical concentration also molecules/cm^3
            c_0 = 1.0/((4.0/3)*np.pi*(R_0**3))
            self.P_ET_DA = 1.0 - np.exp(-1.42 * c_A_mpcm/c_0)


        '''Here FRET from Acceptor to donor, Calculate FRET probabilities if there is a fixed distance between donor and acceptor or not'''
        if self.Photon.CurrentMaterial.FixDist:
            invR6 = 1/((self.Photon.CurrentMaterial.FixDistR*1e-9)**6)
            k_ET = ((9*(2/3.0)*self.Q_A*np.log(10.0))/(128*(self.Photon.CurrentMaterial.n**4) * 6.022e23 * (np.pi**5) * self.tau_A) ) * invR6 * (1e-54) * self.overlAD * 1e17
            self.P_ET_AD = k_ET/(1.0/self.tau_A + k_ET)
        else:
            #Fret Radius
            R_0 = ((9*(2/3.0)*self.Q_A*np.log(10.0))/(128*(self.Photon.CurrentMaterial.n**4) * 6.022e23 * (np.pi**5))) * self.overlAD * 1e17
            R_0 = R_0**(1/6.0)
            #Convert to cm
            print "R_0_AD:",R_0
            R_0 *= 1e-7
            #concentration in molecules/cm^3
            c_D_mpcm = self.c_D * 6.022e17
            #critical concentration also molecules/cm^3
            if R_0 >0:
                c_0 = 1.0/((4.0/3)*np.pi*(R_0**3))
                self.P_ET_AD = 1.0 - np.exp(-1.42 * c_D_mpcm/c_0)
            else:
                self.P_ET_AD = 0.0


        # ax.quiver(coord_A[:,0],coord_A[:,1],coord_A[:,2],self.direction_A[:,0],self.direction_A[:,1],self.direction_A[:,2], length = 0.1)
        # plt.show()

        # plt.ion()
        '''Plot Stuff'''
        # self.fig = plt.figure()
        # self.ax = self.fig.add_subplot(111, projection='3d')
        # self.ax.set_xlabel('x')
        # self.ax.set_ylabel('y')
        # self.ax.set_zlabel('z')
        #
        # for i in self.direction_A:
        #     self.ax.scatter(i[0],i[1],i[2])

        #Generate orientations



    def EmitWavelength(self,Donor = True):
        '''Find a random emission wavelength depending on the emission spectra taking into account that energy can not be created, just cut emission spectrum at the wavelength'''
        if Donor:
            # ind = RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_em_D[:,0],self.Photon.Wavelength)
            # emspect = cp.deepcopy(self.Photon.CurrentMaterial.sp_em_D[ind:,:])
            # emspect[:,1] /= np.sum(emspect[:,1])
            # print ind,self.Photon.Wavelength
            # print self.Photon.CurrentMaterial.sp_em_D

            emspect = RTU.RectifySpectrum(self.Photon.CurrentMaterial.sp_em_D)

            # print emspect, np.sum(emspect[:,1])
            # print emspect
            return np.random.choice(emspect[:,0], p = emspect[:,1])
        else:
            # ind = RTU.FindClosestListElement(self.Photon.CurrentMaterial.sp_em_A[:,0],self.Photon.Wavelength)
            # emspect = cp.deepcopy(self.Photon.CurrentMaterial.sp_em_A[ind:,:])
            # emspect[:,1] /= np.sum(emspect[:,1])

            emspect = RTU.RectifySpectrum(self.Photon.CurrentMaterial.sp_em_A)

            return np.random.choice(emspect[:,0], p = emspect[:,1])



def TESTEMISSION():
    #plt.ion()
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')
   # ax.set_xlim([-1,1])
  #  ax.set_ylim([-1,1])
 #   ax.set_zlim([-1,1])
    count = 0
    nmb = 100000
    for i in range(nmb):
        #v = RTU.GenerateDirectionSineSq()
        v = RTU.GenerateDirectionSine()
        # tm = np.roll(RTU.MakeVectAxis([0,0,1]),2,axis=1)
        # v = np.dot(tm,v)
        # v /= np.linalg.norm(v)
        #ax.scatter(v[0],v[1],v[2])
        #print (np.arccos(v[2])/np.pi)*180.0
        if (np.arccos(v[2]) > (41.81/180.0)*np.pi) and (np.arccos(v[2]) < ((180.0 - 41.81)/180.0)*np.pi):
            count += 1
    print count, nmb- count

#    plt.show()
    raw_input()

# air = Material(1.0, color = 'red')
# MCF = MCFRET([np.array([[-0.5,-0.5,-0.1],[1,2,3]]),'Start',488,air,0,[1,0,0]])
#TESTEMISSION()

# MCF.GenerateMolecules()
