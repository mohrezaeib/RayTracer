import numpy as np
import RTUtils as RTU
import __builtin__

class MaterialObj():
    def __init__(self, name, n, IsAbsorber=False, conc_a=0, conc_b=0, color='blue', IsAligned=False, AlignVect=[1,0,0], AlignS=0.1, R_0=5, Q_D=0.9, Q_A=0.9, tau_D=8, tau_A=8, emax_D=23500, emax_A=54000, sp_abs_D="Coumarin1_Abs.txt", sp_em_D="Coumarin1_Em.txt", sp_abs_A="Coumarin6_Abs.txt", sp_em_A="Coumarin6_Em.txt", AbsEmAngle=0, FixDist=False, FixDistR=1, FixDistNeighbors=1):
        '''alignment thus far only for dye a "acceptor"'''
        self.name = name
        self.n = n
        self.IsAbsorber = IsAbsorber
        self.conc_a = conc_a
        self.conc_b = conc_b

        self.R_0 = R_0

        self.Q_D = Q_D # Donor fluorescence quantum yield
        self.Q_A = Q_A # Acceptor fluorescence quantum yield
        self.tau_A = tau_A # Donor fluorescence lifetime in ns
        self.tau_D = tau_D # Donor Internal conversion lifetime in ns

        self.color = color
        self.IsAligned = IsAligned
        self.AlignVect = AlignVect
        self.AlignS = AlignS
        self.AbsEmAngle = AbsEmAngle

        # https://en.wikipedia.org/wiki/Molar_attenuation_coefficient und tummeltshammer eq2
        self.emax_D = emax_D
        self.emax_A = emax_A

        self.FixDist = FixDist
        self.FixDistR  = FixDistR if FixDistR != 0 else 1  # Ensure FixDistR is not zero by default
        self.FixDistNeighbors = FixDistNeighbors

        setattr(__builtin__, name, self)

        # Read Spectra
        if IsAbsorber:
            '''This part is really important. depending on the spectra the reabsorption process is hugely influenced by the treatment of negative values'''

            self.sp_abs_D = np.maximum(0, np.loadtxt(sp_abs_D))
            self.sp_em_D = np.maximum(0, np.loadtxt(sp_em_D))
            self.sp_abs_A = np.maximum(0, np.loadtxt(sp_abs_A))
            self.sp_em_A = np.maximum(0, np.loadtxt(sp_em_A))

            # Normalize spectra
            self.sp_abs_D[:, 1] /= np.amax(self.sp_abs_D[:, 1])
            self.sp_abs_A[:, 1] /= np.amax(self.sp_abs_A[:, 1])

            # Put in 2e-4 as minimum value
            self.sp_abs_D[:, 1] = np.where(np.logical_and(self.sp_abs_D[:, 1] < 2e-4, self.sp_abs_D[:, 1] > 0), 2e-4, self.sp_abs_D[:, 1])
            self.sp_abs_A[:, 1] = np.where(np.logical_and(self.sp_abs_A[:, 1] < 2e-4, self.sp_abs_A[:, 1] > 0), 2e-4, self.sp_abs_A[:, 1])

            # Kill Values that are too far off
            self.sp_abs_D[:, 1] = np.where(self.sp_abs_D[:, 0] > self.sp_abs_D[np.argmax(self.sp_abs_D[:, 1]), 0] * 1.5, 0.0, self.sp_abs_D[:, 1])
            self.sp_abs_D = self.sp_abs_D[self.sp_abs_D[:, 1] > 0]
            self.sp_abs_D = np.append(self.sp_abs_D, [[self.sp_abs_D[-1, 0] + 1, 0.0]], axis=0)

            self.sp_abs_A[:, 1] = np.where(self.sp_abs_A[:, 0] > self.sp_abs_A[np.argmax(self.sp_abs_A[:, 1]), 0] * 1.5, 0.0, self.sp_abs_A[:, 1])
            self.sp_abs_A = self.sp_abs_A[self.sp_abs_A[:, 1] > 0]
            self.sp_abs_A = np.append(self.sp_abs_A, [[self.sp_abs_A[-1, 0] + 1, 0.0]], axis=0)

            self.sp_em_D[:, 1] = np.where(self.sp_em_D[:, 0] > 2000, 0.0, self.sp_em_D[:, 1])
            self.sp_em_D = self.sp_em_D[self.sp_em_D[:, 1] > 0]
            self.sp_em_A[:, 1] = np.where(self.sp_em_A[:, 0] > 2000, 0.0, self.sp_em_A[:, 1])
            self.sp_em_A = self.sp_em_A[self.sp_em_A[:, 1] > 0]

            # Normalize emission spectra
            self.sp_em_D[:, 1] /= np.sum(self.sp_em_D[:, 1])
            self.sp_em_A[:, 1] /= np.sum(self.sp_em_A[:, 1])

            # Calculate Overlap Integrals of spectra for FRET
            self.ol_DA = self.CalcOverlap(self.sp_em_D, self.sp_abs_A, self.emax_A)
            self.ol_AD = self.CalcOverlap(self.sp_em_A, self.sp_abs_D, self.emax_D)

            print("DID LOAD SPECTRA")
            print(name, "Donor", self.sp_abs_D)
            print(name, "Acceptor", self.sp_abs_A)

    def CalcOverlap(self, EmD, AbsA, eA):
        overl = 0
        for l in range(int(np.amin(EmD[:, 0])), int(np.amax(AbsA[:, 0]))):
            overl += EmD[RTU.FindClosestListElement(EmD[:, 0], l), 1] * (eA * AbsA[RTU.FindClosestListElement(AbsA[:, 0], l), 1]) * ((l)**4)
        return overl
