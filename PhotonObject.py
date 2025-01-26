import copy as cp
import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF


class Photon():
    """
    Photon class responsible for photon initialization, emission, absorption,
    reflection, refraction, focusing, polarization, etc.
    """
    def __init__(
        self,
        Location=[0, 0, 0],
        EmissionMode="Directional Isotropic",
        EmissionDirection=[0, 0, 0],
        EmissionWavelengthType="Sun spectrum AG15g",
        StartingWavelength=450.0,
        StartingMaterial=None,
        RespectPolarization=True,
        StartingLocation=[0, 0, 0],
        ScaleFactor=1.0
    ):
        self.Location = Location
        self.EmissionMode = EmissionMode
        self.EmissionDirection = EmissionDirection
        self.EmissionWavelengthType = EmissionWavelengthType
        self.StartingWavelength = StartingWavelength
        self.StartingMaterial = StartingMaterial
        self.RespectPolarization = RespectPolarization
        self.StartingLocation = StartingLocation
        self.ScaleFactor = ScaleFactor

        # Defaults and placeholders
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

        self.OldPoint = [0, 0, 0]
        self.LineColor = "black"

    def EmitSun(self):
        """
        Emits a photon wavelength sampled from the AM1.5 solar spectrum.
        The AM15.txt file should contain two columns: wavelength (nm), intensity.
        """
        am15 = np.loadtxt("AM15.txt")
        # Convert W/m^2 to a photon flux-like weighting (proportional only).
        # https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
        am15[:, 1] *= am15[:, 0]
        am15 = RTU.RectifySpectrum(am15)
        return np.random.choice(am15[:, 0], p=am15[:, 1])

    def ResetPhoton(self, OperationPlane):
        """
        Reset the photon before a new optical event or run.
        """
        self.GeneratePhoton()
        self.OperationPlane = cp.deepcopy(OperationPlane)
        self.StartingWavelength = cp.deepcopy(self.Wavelength)
        self.Reflected = False

    def FindSign(self):
        """
        Determine the direction sign relative to the plane normal
        (used in reflection/refraction).
        """
        self.sgn = np.sign(-np.dot(self.Direction, self.OperationPlane[1]))

    def GeneratePhoton(self):
        """
        Generate a fresh photon based on the initial parameters like emission mode,
        wavelength type, etc.
        """
        self.WasElement = False
        self.Location = cp.deepcopy(self.StartingLocation)

        # Direction generation
        if self.EmissionMode == "Directional Isotropic":
            self.Direction = RTU.GenerateDirectionSine()
            self.Direction[2] = abs(self.Direction[2])
            tm = np.roll(RTU.MakeVectAxis(self.EmissionDirection), 2, axis=1)
            self.Direction = np.dot(self.Direction, tm)
        elif self.EmissionMode == "Directional Strict":
            self.Direction = cp.deepcopy(self.EmissionDirection)

        self.Direction /= np.linalg.norm(self.Direction)

        # Wavelength generation
        if self.EmissionWavelengthType == "Sun spectrum AG15g":
            self.Wavelength = self.EmitSun()
        else:
            self.Wavelength = cp.deepcopy(self.StartingWavelength)

        print("Wavelength", self.Wavelength)

        # Set current material and initialize polarization
        self.CurrentMaterial = cp.deepcopy(self.StartingMaterial)
        self.UpdatePolarizationRandom(self.Direction)

    def UpdatePolarizationMolOrientation(self, MolOrient, EmissionDirection):
        """
        Update the photon's electric field vector based on a given molecular
        orientation and the emission direction.
        """
        v1 = np.cross(MolOrient, EmissionDirection)
        v2 = np.cross(EmissionDirection, v1)
        v2 /= np.linalg.norm(v2)
        self.ElectricFieldVector = v2

    def UpdatePolarizationRandom(self, EmissionDirection):
        """
        Assign a random polarization vector perpendicular to the emission direction.
        """
        # Find a random unit vector
        v = np.random.randn(3)
        v /= np.linalg.norm(v)
        # Cross product to ensure perpendicularity to EmissionDirection
        self.ElectricFieldVector = np.cross(EmissionDirection, v)
        self.ElectricFieldVector /= np.linalg.norm(self.ElectricFieldVector)

    def UpdatePolarization(self, OldDirection):
        """
        Update the photon's polarization vector to account for changes in direction
        after reflection, refraction, etc.
        """
        newdir = self.Direction
        cross_product = np.cross(newdir, OldDirection)
        norm_cp = np.linalg.norm(cross_product)
        if norm_cp != 0:
            dot_val = np.dot(newdir, OldDirection) / (
                np.linalg.norm(newdir) * np.linalg.norm(OldDirection)
            )
            # Ensure dot_val stays within [-1, 1]
            if -1.0 <= dot_val <= 1.0:
                rot_angle = -np.arccos(dot_val)
                rot_matrix = RTU.RotationMatrix(cross_product, rot_angle)
                self.ElectricFieldVector = np.dot(rot_matrix, self.ElectricFieldVector)
                self.ElectricFieldVector /= np.linalg.norm(self.ElectricFieldVector)

    def PrepareExecution(self):
        """
        Prepare the photon for an operation on the current plane or element.
        """
        self.FindSign()
        self.OldPoint = cp.deepcopy(self.Location)
        self.OldDirection = cp.deepcopy(self.Direction)
        self.LineColor = self.CurrentMaterial.color

    def PostProcessExecution(self):
        """
        Final updates after executing a reflection, refraction, or another operation.
        """
        self.UpdatePolarization(self.OldDirection)
        self.WasElement = True

    def ReflectOnPlane(self):
        """
        Reflect the photon on its current operation plane (like a mirror reflection).
        """
        try:
            denom = np.dot(self.Direction, self.OperationPlane[1])
            num = np.dot((self.OperationPlane[0] - self.Location), self.OperationPlane[1])
            self.Location = (num / denom) * self.Direction + self.Location
            self.Direction = 2 * np.dot(self.OperationPlane[1], self.Direction) * self.OperationPlane[1] - self.Direction
        except ZeroDivisionError:
            print("Ray and plane possibly do not intersect (ReflectOnPlane).")
            raise
        except:
            print("Unknown error in ReflectOnPlane.")
            raise

    def FocusThinLens(self, f):
        """
        Model an ideal thin lens at the plane (self.OperationPlane) with focal length f.
        The photon is refracted as if passing through a thin lens, focusing to a point
        on the focal plane.
        """
        try:
            plane_den = np.dot(self.Direction, self.OperationPlane[1])
            if plane_den == 0:
                print("Ray is parallel to the lens plane; no intersection.")
                return

            # Find intersection with lens plane
            BaseVect = (
                np.dot((self.OperationPlane[0] - self.Location), self.OperationPlane[1])
                / plane_den
            ) * self.Direction + self.Location

            # Two possible focal planes: f in front or behind the lens
            FPlane1 = [self.OperationPlane[0] + f * self.OperationPlane[1], self.OperationPlane[1]]
            FPlane2 = [self.OperationPlane[0] - f * self.OperationPlane[1], self.OperationPlane[1]]

            # Check sides to decide which focal plane to use
            if RTU.CheckSides(self.Location, FPlane1[0], self.OperationPlane):
                FPlane = FPlane2
            else:
                FPlane = FPlane1

            # Project a "parallel" ray from BaseVect
            ParallelRay = [[BaseVect[0], BaseVect[1], self.OperationPlane[0][2]], self.Direction]
            denom_pr = np.dot(ParallelRay[1], FPlane[1])
            if denom_pr == 0:
                print("ParallelRay is parallel to focal plane; no intersection.")
                return

            # Intersection with the chosen focal plane
            num_pr = np.dot((FPlane[0] - ParallelRay[0]), FPlane[1])
            ImgP = (num_pr / denom_pr) * ParallelRay[1] + ParallelRay[0]

            DirectionVector = -(ImgP - BaseVect)
            self.Location = BaseVect
            self.Direction = DirectionVector / np.linalg.norm(DirectionVector)

        except ZeroDivisionError:
            print("Ray and plane possibly do not intersect (FocusThinLens).")
            raise
        except:
            print("Unknown error in FocusThinLens.")
            raise

    def RefractOnPlane(self, n1, n2):
        """
        Use Snell's law to refract the photon through the interface plane
        from material with refractive index n1 to material with refractive index n2.
        """
        try:
            denom = np.dot(self.Direction, self.OperationPlane[1])
            if denom == 0:
                print("Ray and plane possibly do not intersect (RefractOnPlane).")
                return
            num = np.dot((self.OperationPlane[0] - self.Location), self.OperationPlane[1])
            self.Location = (num / denom) * self.Direction + self.Location

            ratio = n1 / n2
            dotval = -np.dot(self.Direction, self.OperationPlane[1])
            under_sqrt = 1 - (ratio ** 2) * (1 - (dotval ** 2))
            # Check for total internal reflection
            if under_sqrt < 0:
                print("Total internal reflection in RefractOnPlane.")
                # If total internal reflection, reflect
                self.ReflectOnPlane()
            else:
                self.Direction = ratio * self.Direction + (
                    ratio * dotval - self.sgn * np.sqrt(under_sqrt)
                ) * self.OperationPlane[1]

        except ZeroDivisionError:
            print("Ray and plane possibly do not intersect (RefractOnPlane) - ZeroDivisionError.")
            raise
        except:
            print("Unknown error in RefractOnPlane.")
            raise

    def CheckForAbsorption(self, OldPoint, SmallestDistance):
        """
        Check if the photon is absorbed by the current material over SmallestDistance.
        If so, handle the absorption event, including potential FRET emission.
        """
        self.OldPoint = OldPoint

        if self.CurrentMaterial.IsAbsorber:
            if self.DistanceLeftToAbsorption > SmallestDistance:
                self.DistanceLeftToAbsorption -= SmallestDistance
                self.Status = 'Start'
            else:
                print("ABSORBED")
                self.BeenAbsorbed[0] = True
                self.WasElement = False
                self.Status = 'Absorbed'

                # Move photon location the remaining distance
                restdirect_v = self.PointCy - self.Location
                restdirect_v /= np.linalg.norm(restdirect_v)
                newp = restdirect_v * self.DistanceLeftToAbsorption + self.Location

                self.Location = newp

                # Probability distribution for donor vs. acceptor absorption
                b_d = (
                    self.CurrentMaterial.conc_b
                    * 1e-3
                    * np.log(10)
                    * self.CurrentMaterial.emax_D
                    * np.maximum(
                        1e-99,
                        self.CurrentMaterial.sp_abs_D[
                            RTU.FindClosestListElement(
                                self.CurrentMaterial.sp_abs_D[:, 0],
                                self.Wavelength
                            ),
                            1
                        ]
                    )
                )
                b_a = (
                    self.CurrentMaterial.conc_a
                    * 1e-3
                    * np.log(10)
                    * self.CurrentMaterial.emax_A
                    * np.maximum(
                        1e-99,
                        self.CurrentMaterial.sp_abs_A[
                            RTU.FindClosestListElement(
                                self.CurrentMaterial.sp_abs_A[:, 0],
                                self.Wavelength
                            ),
                            1
                        ]
                    )
                )

                if self.CurrentMaterial.IsAligned:
                    print("Overlap Integral", self.Psi)
                    b_a *= self.Psi * 3  # Factor of 3 for aligned absorption

                # Absorbed by donor or acceptor?
                if np.random.rand() < (b_d / (b_d + b_a)):
                    print("absorbed by donor")
                    fret = MCF.MCFRET(self, True, self.RespectPolarization)
                    self = fret.RunSimulation()
                    self.LamBeerDistance()
                else:
                    print("absorbed by acceptor")
                    fret = MCF.MCFRET(self, False, self.RespectPolarization)
                    self = fret.RunSimulation()
                    print(self.Wavelength)
                    self.LamBeerDistance()

        else:
            self.Status = 'Start'
            self.WasElement = True

    def LambertBeerFullIntegral(self, n_E):
        """
        Calculate the absorption integral for Lambert-Beer with aligned dyes.
        This is a direct numeric/symbolic expression. For extremely small AlignS,
        fallback to cos^2 alignment if needed.
        """
        x, y, z = n_E
        if self.CurrentMaterial.AlignS > 1e-10:
            ol = np.abs(
                1
                / self.CurrentMaterial.AlignS
                * (
                    np.exp(-0.5 * self.CurrentMaterial.AlignS)
                    * (
                        (-0.2349964007466563j)
                        * (x ** 2 + y ** 2 + (2 / 3.0) * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.CurrentMaterial.AlignS
                            - (0.7071067811865475j) * self.CurrentMaterial.AlignS
                        )
                        + (0.2349964007466563j)
                        * (x ** 2 + y ** 2 + (2 / 3.0) * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.CurrentMaterial.AlignS
                            + (0.7071067811865475j) * self.CurrentMaterial.AlignS
                        )
                        + (
                            0.4699928014933126 * x ** 2
                            + 0.4699928014933126 * y ** 2
                            + 0.313328534328875 * z ** 2
                        )
                        * special.erfi(
                            0.7071067811865475 * self.CurrentMaterial.AlignS
                        )
                    )
                    + np.exp(-4.5 * self.CurrentMaterial.AlignS ** 2)
                    * (
                        (0.07833213358221876j)
                        * (x ** 2 + y ** 2 - 2 * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.CurrentMaterial.AlignS
                            - (2.1213203435596424j) * self.CurrentMaterial.AlignS
                        )
                        - (0.07833213358221876j)
                        * (x ** 2 + y ** 2 - 2 * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.CurrentMaterial.AlignS
                            + (2.1213203435596424j) * self.CurrentMaterial.AlignS
                        )
                        + (
                            -0.15666426716443752 * x ** 2
                            - 0.15666426716443752 * y ** 2
                            + 0.31332853432887503 * z ** 2
                        )
                        * special.erfi(
                            2.1213203435596424 * self.CurrentMaterial.AlignS
                        )
                    )
                )
            )
        else:
            # For extremely small alignment, approximate with cos^2
            print("Fastlane", n_E)
            ol = (np.dot(n_E, [0, 0, 1])) ** 2

        return 0.0 if ol < 1e-10 else ol

    def LamBeerDistance(self):
        """
        Randomly assign the photon's free path length in an absorbing medium,
        according to an exponential distribution. Called when crossing a new medium.
        """
        if self.CurrentMaterial.IsAbsorber:
            idx_absA = RTU.FindClosestListElement(
                self.CurrentMaterial.sp_abs_A[:, 0], self.Wavelength
            )
            b_a = (
                self.CurrentMaterial.conc_a
                * 1e-3
                * self.CurrentMaterial.emax_A
                * np.log(10)
                * np.maximum(1e-99, self.CurrentMaterial.sp_abs_A[idx_absA, 1])
            )

            idx_absD = RTU.FindClosestListElement(
                self.CurrentMaterial.sp_abs_D[:, 0], self.Wavelength
            )
            b_d = (
                self.CurrentMaterial.conc_b
                * 1e-3
                * self.CurrentMaterial.emax_D
                * np.log(10)
                * np.maximum(1e-99, self.CurrentMaterial.sp_abs_D[idx_absD, 1])
            )

            if self.CurrentMaterial.IsAligned:
                tm = np.roll(RTU.MakeVectAxis(self.CurrentMaterial.AlignVect), 2, axis=1)
                n_E = np.dot(self.ElectricFieldVector, tm)
                n_E /= np.linalg.norm(n_E)
                self.Psi = self.LambertBeerFullIntegral(n_E)
                print("Overlap Integral", self.Psi)
                b_a *= self.Psi * 3.0

            total_abs_coeff = b_a + b_d
            if total_abs_coeff == 0:
                scale_val = 1e12  # effectively infinite
            else:
                scale_val = 1.0 / total_abs_coeff

            self.DistanceLeftToAbsorption = (1 / self.ScaleFactor) * np.random.exponential(
                scale=scale_val
            )
        else:
            self.DistanceLeftToAbsorption = 0.0

        print("NEW FREE LENGTH", self.DistanceLeftToAbsorption)
