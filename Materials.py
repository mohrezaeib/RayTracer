import numpy as np
import builtins  # In Python 3, use 'builtins' instead of '__builtin__'
import RTUtils as RTU


class MaterialObj():
    def __init__(
        self,
        name,
        n,
        IsAbsorber=False,
        conc_a=0,
        conc_b=0,
        color='blue',
        IsAligned=False,
        AlignVect=[1, 0, 0],
        AlignS=0.1,
        R_0=5,
        Q_D=0.9,
        Q_A=0.9,
        tau_D=8,
        tau_A=8,
        emax_D=23500,
        emax_A=54000,
        sp_abs_D="Coumarin1_Abs.txt",
        sp_em_D="Coumarin1_Em.txt",
        sp_abs_A="Coumarin6_Abs.txt",
        sp_em_A="Coumarin6_Em.txt",
        AbsEmAngle=0,
        FixDist=False,
        FixDistR=1,
        FixDistNeighbors=1
    ):
        """
        MaterialObj class defines parameters and behaviors of a material,
        especially dyes in FRET simulations.

        name: Material name (string)
        n: Refractive index
        IsAbsorber: If True, treat this material as a dye absorber
        conc_a: Concentration of acceptor dyes
        conc_b: Concentration of donor dyes
        color: Color for visualization
        IsAligned: If True, there's an orientation for the absorber
        AlignVect: Orientation vector for the dye
        AlignS: Order parameter for alignment
        R_0: Förster radius
        Q_D: Donor fluorescence quantum yield
        Q_A: Acceptor fluorescence quantum yield
        tau_D: Donor fluorescence lifetime in ns
        tau_A: Acceptor fluorescence lifetime in ns
        emax_D: Molar extinction coefficient for donor (maximum)
        emax_A: Molar extinction coefficient for acceptor (maximum)
        sp_abs_D: Filename for donor absorption spectrum
        sp_em_D: Filename for donor emission spectrum
        sp_abs_A: Filename for acceptor absorption spectrum
        sp_em_A: Filename for acceptor emission spectrum
        AbsEmAngle: Angle between absorption and emission dipoles
        FixDist: If True, a fixed distance is enforced in FRET calculations
        FixDistR: Distance to be used if FixDist is True
        FixDistNeighbors: Number of nearest neighbors to consider in FRET
        """

        self.name = name
        self.n = n
        self.IsAbsorber = IsAbsorber
        self.conc_a = conc_a
        self.conc_b = conc_b

        self.R_0 = R_0
        self.Q_D = Q_D  # Donor fluorescence quantum yield
        self.Q_A = Q_A  # Acceptor fluorescence quantum yield
        self.tau_A = tau_A  # Acceptor fluorescence lifetime in ns
        self.tau_D = tau_D  # Donor fluorescence lifetime in ns

        self.color = color
        self.IsAligned = IsAligned
        self.AlignVect = AlignVect
        self.AlignS = AlignS
        self.AbsEmAngle = AbsEmAngle

        self.emax_D = emax_D
        self.emax_A = emax_A

        self.FixDist = FixDist
        # Ensure FixDistR is never zero
        self.FixDistR = FixDistR if FixDistR != 0 else 1
        self.FixDistNeighbors = FixDistNeighbors

        # In Python 3, use builtins instead of __builtin__
        setattr(builtins, name, self)

        # If material is an absorber, load and normalize relevant spectra
        if IsAbsorber:
            # Load data from text files, ensuring no negative values
            self.sp_abs_D = np.maximum(0, np.loadtxt(sp_abs_D))
            self.sp_em_D = np.maximum(0, np.loadtxt(sp_em_D))
            self.sp_abs_A = np.maximum(0, np.loadtxt(sp_abs_A))
            self.sp_em_A = np.maximum(0, np.loadtxt(sp_em_A))

            # Normalize absorption (Donor, Acceptor)
            self.sp_abs_D[:, 1] /= np.amax(self.sp_abs_D[:, 1])
            self.sp_abs_A[:, 1] /= np.amax(self.sp_abs_A[:, 1])

            # Enforce a small minimum for non-zero absorption
            self.sp_abs_D[:, 1] = np.where(
                np.logical_and(self.sp_abs_D[:, 1] < 2e-4, self.sp_abs_D[:, 1] > 0),
                2e-4,
                self.sp_abs_D[:, 1]
            )
            self.sp_abs_A[:, 1] = np.where(
                np.logical_and(self.sp_abs_A[:, 1] < 2e-4, self.sp_abs_A[:, 1] > 0),
                2e-4,
                self.sp_abs_A[:, 1]
            )

            # Suppress absorption far beyond peak (>1.5x peak wavelength)
            peak_D = self.sp_abs_D[np.argmax(self.sp_abs_D[:, 1]), 0]
            self.sp_abs_D[:, 1] = np.where(
                self.sp_abs_D[:, 0] > peak_D * 1.5,
                0.0,
                self.sp_abs_D[:, 1]
            )
            # Filter out zeroed rows and add a trailing zero row
            self.sp_abs_D = self.sp_abs_D[self.sp_abs_D[:, 1] > 0]
            self.sp_abs_D = np.append(
                self.sp_abs_D,
                [[self.sp_abs_D[-1, 0] + 1, 0.0]],
                axis=0
            )

            peak_A = self.sp_abs_A[np.argmax(self.sp_abs_A[:, 1]), 0]
            self.sp_abs_A[:, 1] = np.where(
                self.sp_abs_A[:, 0] > peak_A * 1.5,
                0.0,
                self.sp_abs_A[:, 1]
            )
            self.sp_abs_A = self.sp_abs_A[self.sp_abs_A[:, 1] > 0]
            self.sp_abs_A = np.append(
                self.sp_abs_A,
                [[self.sp_abs_A[-1, 0] + 1, 0.0]],
                axis=0
            )

            # Limit emission data to below 2000 nm
            self.sp_em_D[:, 1] = np.where(
                self.sp_em_D[:, 0] > 2000,
                0.0,
                self.sp_em_D[:, 1]
            )
            self.sp_em_D = self.sp_em_D[self.sp_em_D[:, 1] > 0]
            self.sp_em_A[:, 1] = np.where(
                self.sp_em_A[:, 0] > 2000,
                0.0,
                self.sp_em_A[:, 1]
            )
            self.sp_em_A = self.sp_em_A[self.sp_em_A[:, 1] > 0]

            # Normalize emission spectra
            self.sp_em_D[:, 1] /= np.sum(self.sp_em_D[:, 1])
            self.sp_em_A[:, 1] /= np.sum(self.sp_em_A[:, 1])

            # Calculate overlap integrals for FRET
            self.ol_DA = self.CalcOverlap(self.sp_em_D, self.sp_abs_A, self.emax_A)
            self.ol_AD = self.CalcOverlap(self.sp_em_A, self.sp_abs_D, self.emax_D)

            print("DID LOAD SPECTRA")
            print(name, "Donor", self.sp_abs_D)
            print(name, "Acceptor", self.sp_abs_A)

    def CalcOverlap(self, EmD, AbsA, eA):
        """
        Calculate the spectral overlap for the emission of one dye and
        absorption of another, weighted by wavelength^4.
        """
        overl = 0
        # Iterate over integer wavelengths from the min of EmD to the max of AbsA
        for lam in range(int(np.amin(EmD[:, 0])), int(np.amax(AbsA[:, 0]))):
            em_idx = RTU.FindClosestListElement(EmD[:, 0], lam)
            ab_idx = RTU.FindClosestListElement(AbsA[:, 0], lam)
            overl += EmD[em_idx, 1] * (eA * AbsA[ab_idx, 1]) * (lam ** 4)
        return overl
