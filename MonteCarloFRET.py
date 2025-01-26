import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# The following import is fine in Python 3, though Axes3D is often optional:
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import copy as cp
# from RayTrace import Material
import RTUtils as RTU

"""
Monte Carlo Simulation for FRET (Förster Resonance Energy Transfer) dynamics.
A single excited molecule is analyzed, and the energy transfer probabilities
are calculated.
"""

class MCFRET(object):
    def __init__(self, Photon, InDonor=True, EmPol=True):
        # self.R_0 = Photon.CurrentMaterial.R_0  # FoersterRadius in nm
        self.c_D = Photon.CurrentMaterial.conc_b  # Donor concentration in mmol/L
        self.c_A = Photon.CurrentMaterial.conc_a  # Acceptor concentration in mmol/L
        self.Q_D = Photon.CurrentMaterial.Q_D     # Donor fluorescence quantum yield
        self.Q_A = Photon.CurrentMaterial.Q_A     # Acceptor fluorescence quantum yield
        self.tau_A = Photon.CurrentMaterial.tau_A # Acceptor fluorescence lifetime in ns
        self.tau_D = Photon.CurrentMaterial.tau_D # Donor fluorescence lifetime in ns

        self.AlignS = Photon.CurrentMaterial.AlignS
        self.AlignVect = Photon.CurrentMaterial.AlignVect
        self.IsAligned = Photon.CurrentMaterial.IsAligned

        self.EmPol = EmPol

        # Rate constants
        self.k_D = self.Q_D / self.tau_D
        self.k_A = self.Q_A / self.tau_A
        self.k_D_IC = (1 - self.Q_D) / self.tau_D  # Internal conversion (Donor)
        self.k_A_IC = (1 - self.Q_A) / self.tau_A  # Internal conversion (Acceptor)

        self.Photon = Photon
        self.InDonor = InDonor

        self.overlDA = self.Photon.CurrentMaterial.ol_DA  # Overlap integral Donor->Acceptor
        self.overlAD = self.Photon.CurrentMaterial.ol_AD  # Overlap integral Acceptor->Donor

    def CalcOverlap(self, EmD, AbsA, eA):
        """
        Numerical approximation of the overlap integral by summing discrete points.
        integrated over integer wavelengths from min(EmD) to max(AbsA).
        """
        overl = 0
        for l in range(int(np.amin(EmD[:, 0])), int(np.amax(AbsA[:, 0]))):
            idx_em = RTU.FindClosestListElement(EmD[:, 0], l)
            idx_ab = RTU.FindClosestListElement(AbsA[:, 0], l)
            # Weighted by λ^4
            overl += EmD[idx_em, 1] * (eA * AbsA[idx_ab, 1]) * (l ** 4)
        return overl

    def RunSimulation(self):
        """
        Perform one "event" of the photon in a donor/acceptor medium.
        Possible events:
          1) Fluorescence
          2) Internal Conversion
          3) FRET (Donor -> Acceptor or Acceptor -> Donor)
        """
        if self.Photon.CurrentMaterial.conc_a > 0:
            self.GenerateMolecules()

            # Probabilities for Donor
            #  0: Fluorescence
            #  1: Internal Conversion
            #  2: FRET
            p = np.zeros(3)
            p[0] = (1 - self.P_ET_DA) * self.Q_D
            p[1] = (1 - self.P_ET_DA) * (1 - self.Q_D)
            p[2] = self.P_ET_DA

            # Probabilities for Acceptor
            pA = np.zeros(3)
            pA[0] = (1 - self.P_ET_AD) * self.Q_A
            pA[1] = (1 - self.P_ET_AD) * (1 - self.Q_A)
            pA[2] = self.P_ET_AD

            print("probabilities Donor: Fl: " + str(p[0]) +
                  " IC: " + str(p[1]) +
                  " FRET: " + str(np.sum(p[2:])) +
                  " Overlap Integral: " + str(self.overlDA))
            print("probabilities Acceptor: Fl: " + str(pA[0]) +
                  " IC: " + str(pA[1]) +
                  " FRET: " + str(np.sum(pA[2:])) +
                  " Overlap Integral: " + str(self.overlAD))
        else:
            # If there's no acceptor, only two possible events for Donor
            # 0 -> Fluorescence
            # 1 -> Internal Conversion
            p = np.array([self.Photon.CurrentMaterial.Q_D,
                          1 - self.Photon.CurrentMaterial.Q_D])

        # Now choose the outcome for the donor (if InDonor)
        # or reuse p if there's no acceptor concentration
        Pathway = np.random.choice(p.shape[0], 1, p=p)[0]

        # Pathway 0: Fluorescence
        # Pathway 1: Internal Conversion
        # Pathway 2: FRET
        if Pathway == 0 and self.InDonor:
            # Fluorescence from Donor
            self.Photon.Direction = RTU.GenerateDirectionSine()
            self.Photon.Wavelength = self.EmitWavelength(True)
            self.Photon.Status = 'ReemittedDonor'
            print("ReemittedDonor")
            self.Photon.UpdatePolarizationRandom(self.Photon.Direction)

        elif Pathway == 1 and self.InDonor:
            # Internal conversion (energy lost)
            self.Photon.Status = 'InternalConverted'
            print("InternalConverted Donor")

        elif Pathway > 1 or not self.InDonor:
            # FRET event (Donor -> Acceptor) or Already in Acceptor
            if Pathway > 1 and self.InDonor:
                print("FRET from Donor to Acceptor")
            else:
                print("In Acceptor")

            # Next step: choose acceptor pathway
            PathwayAcceptor = np.random.choice(3, 1, p=pA)[0]
            if PathwayAcceptor == 0:
                # Fluorescence from Acceptor
                print("Reemitted Acceptor")
                if self.IsAligned:
                    # Emit with a sine-squared distribution around the acceptor alignment
                    remdir = RTU.GenerateDirectionGauss(self.Photon)
                    tm = np.roll(RTU.MakeVectAxis(remdir), 2, axis=1)

                    if self.Photon.CurrentMaterial.AbsEmAngle > 0:
                        th = (self.Photon.CurrentMaterial.AbsEmAngle / 180.0) * np.pi
                        ph = np.random.rand() * 2 * np.pi
                        v = [0, 0, 0]
                        v[0] = np.sin(th) * np.cos(ph)
                        v[1] = np.sin(th) * np.sin(ph)
                        v[2] = np.cos(th)

                        etm = np.roll(RTU.MakeVectAxis(v), 2, axis=1)
                        self.Photon.Direction = np.dot(etm, self.Photon.Direction)

                    self.Photon.Direction = np.dot(tm, self.Photon.Direction)
                    self.Photon.Direction /= np.linalg.norm(self.Photon.Direction)

                    if self.EmPol:
                        self.Photon.UpdatePolarizationMolOrientation(remdir, self.Photon.Direction)
                    else:
                        self.Photon.UpdatePolarizationRandom(self.Photon.Direction)

                    self.Photon.Wavelength = self.EmitWavelength(False)
                    self.Photon.Status = 'ReemittedAcceptor'

                else:
                    # Isotropic emission from the acceptor
                    self.Photon.Direction = RTU.GenerateDirectionSine()
                    self.Photon.UpdatePolarizationRandom(self.Photon.Direction)
                    self.Photon.Wavelength = self.EmitWavelength(False)
                    self.Photon.Status = 'ReemittedAcceptor'

                print(self.Photon.Wavelength)

            elif PathwayAcceptor == 1:
                # Internal conversion in Acceptor
                print("Internal Converted Acceptor")
                self.Photon.Status = 'InternalConverted'

            else:
                # FRET from Acceptor back to Donor
                print("FRET from Acceptor to Donor")
                self.InDonor = True
                return self.RunSimulation()

        print("Wavelength:", self.Photon.Wavelength)
        return self.Photon

    def GenerateMolecules(self):
        """
        FRET probabilities from:
        T. Förster (1959).
        For fixed or variable distances between donor and acceptor.
        """
        if self.Photon.CurrentMaterial.FixDistR == 0:
            # Edge case protection
            raise ValueError("FixDistR is zero, which will cause a division by zero error.")

        # Donor -> Acceptor
        if not self.Photon.CurrentMaterial.FixDist:
            # R0 calculation in nm, then converted to cm
            R_0 = (
                (9 * (2 / 3.0) * self.Q_D * np.log(10.0)) /
                (128 * (self.Photon.CurrentMaterial.n ** 4) * 6.022e23 * (np.pi ** 5))
            ) * self.overlDA * 1e17
            R_0 = R_0 ** (1 / 6.0)
            print("R_0_DA:", R_0)
            # Convert R_0 from nm to cm
            R_0 *= 1e-7

            # Concentration in molecules/cm^3
            c_A_mpcm = self.c_A * 6.022e17
            c_0 = 1.0 / ((4.0 / 3) * np.pi * (R_0 ** 3))
            self.P_ET_DA = 1.0 - np.exp(-1.42 * c_A_mpcm / c_0)
        else:
            # Fixed distance approach for Donor -> Acceptor
            invR6 = 1 / ((self.Photon.CurrentMaterial.FixDistR * 1e-9) ** 6)
            k_ET = (
                (9 * (2 / 3.0) * self.Q_D * np.log(10.0)) /
                (
                    128
                    * (self.Photon.CurrentMaterial.n ** 4)
                    * 6.022e23
                    * (np.pi ** 5)
                    * self.tau_D
                )
            ) * invR6 * (1e-54) * self.overlDA * 1e17
            self.P_ET_DA = k_ET / (1.0 / self.tau_D + k_ET)

        # Acceptor -> Donor
        if not self.Photon.CurrentMaterial.FixDist:
            R_0 = (
                (9 * (2 / 3.0) * self.Q_A * np.log(10.0)) /
                (128 * (self.Photon.CurrentMaterial.n ** 4) * 6.022e23 * (np.pi ** 5))
            ) * self.overlAD * 1e17
            R_0 = R_0 ** (1 / 6.0)
            print("R_0_AD:", R_0)
            R_0 *= 1e-7
            c_D_mpcm = self.c_D * 6.022e17
            if R_0 > 0:
                c_0 = 1.0 / ((4.0 / 3) * np.pi * (R_0 ** 3))
                self.P_ET_AD = 1.0 - np.exp(-1.42 * c_D_mpcm / c_0)
            else:
                self.P_ET_AD = 0.0
        else:
            invR6 = 1 / ((self.Photon.CurrentMaterial.FixDistR * 1e-9) ** 6)
            k_ET = (
                (9 * (2 / 3.0) * self.Q_A * np.log(10.0)) /
                (
                    128
                    * (self.Photon.CurrentMaterial.n ** 4)
                    * 6.022e23
                    * (np.pi ** 5)
                    * self.tau_A
                )
            ) * invR6 * (1e-54) * self.overlAD * 1e17
            self.P_ET_AD = k_ET / (1.0 / self.tau_A + k_ET)

    def EmitWavelength(self, Donor=True):
        """
        Randomly select a new emission wavelength from the relevant emission
        spectrum. The spectrum is assumed to be normalized. If you wish to
        account for energy conservation (shorter than absorbed wavelength),
        you could modify the spectrum accordingly.
        """
        if Donor:
            emspect = RTU.RectifySpectrum(self.Photon.CurrentMaterial.sp_em_D)
            return np.random.choice(emspect[:, 0], p=emspect[:, 1])
        else:
            emspect = RTU.RectifySpectrum(self.Photon.CurrentMaterial.sp_em_A)
            return np.random.choice(emspect[:, 0], p=emspect[:, 1])


def TESTEMISSION():
    """
    Simple test function that samples random directions and checks if
    they fall within a certain range of angles.
    """
    count = 0
    nmb = 100000
    for _ in range(nmb):
        v = RTU.GenerateDirectionSine()
        # Check polar angle
        angle = np.arccos(v[2])  # in radians
        deg = (angle / np.pi) * 180.0
        if 41.81 < deg < (180.0 - 41.81):
            count += 1
    print(count, nmb - count)
    input("Press Enter to continue...")
