import copy as cp
import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF
import Geometries
import Elements


class OpticalElement():
    def __init__(
        self,
        SupportVector,
        DirectionVector,
        ElementProperties,
        ElementGeometry,
        ElementSize,
        Photon
    ):
        """
        Initialize elements. It would be nicer to directly parse elements in the
        correct way or make this part of the parser code.

        Parameters
        ----------
        SupportVector : array-like
            Position vector (x, y, z) of the element's support/reference point.
        DirectionVector : array-like
            Orientation of the element; will be normalized.
        ElementProperties : list
            Properties describing the optical element's behavior (type, materials, etc.).
        ElementGeometry : int
            Index describing the geometry shape.
        ElementSize : array-like
            Size parameters for the geometry (radius, height, width, etc.).
        Photon : PhotonObject
            Reference to the Photon that interacts with this element.
        """
        self.SupportVector = SupportVector
        self.DirectionVector = DirectionVector / np.linalg.norm(DirectionVector)
        self.ElementGeometry = ElementGeometry
        self.ElementSize = ElementSize
        self.ElementProperties = ElementProperties
        self.Photon = Photon
        self.Photon.Reflected = False

        # Translate parsed data into readable variables
        shapes = ["None", "Rectangle", "Circle", "Cylinder", "Triangle"]
        self.Shape = shapes[self.ElementGeometry]

        # Create the geometry object
        if self.Shape == "Circle":
            self.Geometry = Geometries.Circle(
                self.ElementSize[0],
                self.SupportVector,
                self.DirectionVector,
                self.Photon
            )
        elif self.Shape == "Cylinder":
            self.Geometry = Geometries.Cylinder(
                self.ElementSize[0],
                self.SupportVector,
                self.DirectionVector,
                self.Photon
            )
        elif self.Shape == "Rectangle":
            # If a normal vector is provided in ElementSize, use it
            if len(self.ElementSize) == 2:
                NV = [0, 0, 0]
            else:
                NV = self.ElementSize[2]
            self.Geometry = Geometries.Rectangle(
                self.ElementSize[0],
                self.ElementSize[1],
                self.SupportVector,
                self.DirectionVector,
                self.Photon,
                NV
            )
        elif self.Shape == "Triangle":
            self.Geometry = Geometries.Triangle(
                self.ElementSize[0],
                self.ElementSize[1],
                self.ElementSize[2],
                self.Photon
            )
        else:
            self.Geometry = None

        # Create the element type object
        types = ["None", "Mirror", "Interface", "Lens", "Dump"]
        self.Type = types[self.ElementProperties[0]]
        self.ElementType = Elements.Mirror(self.Photon)

        if self.Type == "Interface":
            self.ElementType = Elements.Interface(
                self.ElementProperties[1],
                self.ElementProperties[2],
                self.Photon
            )
        elif self.Type == "Lens":
            self.ElementType = Elements.Lens(
                self.ElementProperties[1],
                self.Photon
            )
        elif self.Type == "Dump":
            self.ElementType = Elements.Dump(
                self.ElementProperties[4],
                self.Photon
            )

    def FindIntersections(self):
        """
        Find intersection points of the Photon's ray with the element geometry.
        If the ray is parallel to the geometry's normal, return [np.nan, np.nan, np.nan].
        """
        if np.dot(self.Photon.Direction, self.DirectionVector) != 0:
            return self.Geometry.FindIntersections()
        else:
            return [np.nan, np.nan, np.nan]

    def FindSign(self, Plane):
        """
        Determine the sign of the normal vector for reflection/refraction.
        """
        self.Photon.sgn = np.sign(-np.dot(self.Photon.Direction, Plane[1]))

    def ExecuteCurrentElement(self):
        """
        Execute the optical operation (reflection, refraction, etc.) associated
        with this element, based on the Photon and the element type.
        """
        self.Photon.OperationPlane = self.Geometry.OperationPlane()

        self.FindSign(self.Photon.OperationPlane)

        self.OldPoint = cp.deepcopy(self.Photon.Location)
        self.OldDirection = cp.deepcopy(self.Photon.Direction)
        self.LineColor = self.Photon.CurrentMaterial.color

        # Execute the operation (Mirror, Interface, Lens, Dump, etc.)
        self.ElementType.ExecuteOperation()

        # Update photon polarization after the operation
        self.Photon.UpdatePolarization(self.OldDirection)
        self.Photon.WasElement = True

    def CheckForAbsorption(self, OldPoint, SmallestDistance):
        """
        Check if the Photon is absorbed in the current material. If so,
        stop the Photon and handle any resulting FRET emission processes.
        """
        self.OldPoint = OldPoint
        self.LineColor = "blue"

        if self.Photon.CurrentMaterial.IsAbsorber:
            if self.Photon.DistanceLeftToAbsorption > SmallestDistance:
                # Not absorbed yet
                self.Photon.DistanceLeftToAbsorption -= SmallestDistance
                self.Photon.Status = 'Start'
            else:
                print("ABSORBED")
                self.Photon.BeenAbsorbed[0] = True
                self.Photon.WasElement = False
                self.Photon.Status = 'Absorbed'

                # Move the photon location along the remaining distance
                restdirect_v = self.Photon.PointCy - self.Photon.Location
                restdirect_v /= np.linalg.norm(restdirect_v)
                newp = restdirect_v * self.Photon.DistanceLeftToAbsorption + self.Photon.Location

                self.Photon.Location = newp

                # Determine whether donor or acceptor absorbs
                b_d = (
                    self.Photon.CurrentMaterial.conc_b
                    * 1e-3
                    * np.log(10)
                    * self.Photon.CurrentMaterial.emax_D
                    * np.maximum(
                        1e-99,
                        self.Photon.CurrentMaterial.sp_abs_D[
                            RTU.FindClosestListElement(
                                self.Photon.CurrentMaterial.sp_abs_D[:, 0],
                                self.Photon.Wavelength
                            ),
                            1
                        ]
                    )
                )
                b_a = (
                    self.Photon.CurrentMaterial.conc_a
                    * 1e-3
                    * np.log(10)
                    * self.Photon.CurrentMaterial.emax_A
                    * np.maximum(
                        1e-99,
                        self.Photon.CurrentMaterial.sp_abs_A[
                            RTU.FindClosestListElement(
                                self.Photon.CurrentMaterial.sp_abs_A[:, 0],
                                self.Photon.Wavelength
                            ),
                            1
                        ]
                    )
                )

                if self.Photon.CurrentMaterial.IsAligned:
                    print("Overlap Integral", self.Photon.Psi)
                    # The factor of 3 accounts for isotropic vs. aligned absorption
                    b_a *= self.Photon.Psi * 3

                # Probability photon is absorbed by donor
                if np.random.rand() < (b_d / (b_d + b_a)):
                    print("absorbed by donor")
                    fret = MCF.MCFRET(self.Photon, True, self.Photon.RespectPolarization)
                    self.Photon = fret.RunSimulation()
                    self.LamBeerDistance()
                else:
                    print("absorbed by acceptor")
                    fret = MCF.MCFRET(self.Photon, False, self.Photon.RespectPolarization)
                    self.Photon = fret.RunSimulation()
                    print(self.Photon.Wavelength)
                    self.LamBeerDistance()
        else:
            self.Photon.Status = 'Start'
            self.Photon.WasElement = True

    def LambertBeerFullIntegral(self, n_E):
        """
        Calculate the absorption integral for Lambert-Beer under aligned dyes.
        This is a hard-coded expression from a symbolic/numerical solution.
        """
        x, y, z = n_E
        ol = abs(
            1 / self.Photon.CurrentMaterial.AlignS
            * (
                np.exp(-0.5 * self.Photon.CurrentMaterial.AlignS)
                * (
                    (
                        (-0.2349964007466563j)
                        * (x ** 2 + y ** 2 + (2 / 3.0) * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.Photon.CurrentMaterial.AlignS
                            - (0.7071067811865475j) * self.Photon.CurrentMaterial.AlignS
                        )
                    )
                    + (
                        (0.2349964007466563j)
                        * (x ** 2 + y ** 2 + (2 / 3.0) * z ** 2)
                        * special.erf(
                            2.221441469079183 / self.Photon.CurrentMaterial.AlignS
                            + (0.7071067811865475j) * self.Photon.CurrentMaterial.AlignS
                        )
                    )
                    + (
                        0.4699928014933126 * x ** 2
                        + 0.4699928014933126 * y ** 2
                        + 0.313328534328875 * z ** 2
                    )
                    * special.erfi(
                        0.7071067811865475 * self.Photon.CurrentMaterial.AlignS
                    )
                )
                + np.exp(-4.5 * self.Photon.CurrentMaterial.AlignS ** 2)
                * (
                    (0.07833213358221876j)
                    * (x ** 2 + y ** 2 - 2 * z ** 2)
                    * special.erf(
                        2.221441469079183 / self.Photon.CurrentMaterial.AlignS
                        - (2.1213203435596424j) * self.Photon.CurrentMaterial.AlignS
                    )
                    - (0.07833213358221876j)
                    * (x ** 2 + y ** 2 - 2 * z ** 2)
                    * special.erf(
                        2.221441469079183 / self.Photon.CurrentMaterial.AlignS
                        + (2.1213203435596424j) * self.Photon.CurrentMaterial.AlignS
                    )
                    + (
                        -0.15666426716443752 * x ** 2
                        - 0.15666426716443752 * y ** 2
                        + 0.31332853432887503 * z ** 2
                    )
                    * special.erfi(
                        2.1213203435596424 * self.Photon.CurrentMaterial.AlignS
                    )
                )
            )
        )

        # Provide a small cutoff
        if ol < 1e-10:
            return 0.0
        else:
            return ol

    def LamBeerDistance(self):
        """
        Randomly assign a distance the photon can travel before absorption
        in the current material, following an exponential (Lambert-Beer).
        Only executed when the photon crosses into an absorbing medium.
        """
        if self.Photon.CurrentMaterial.IsAbsorber:
            b_a = (
                self.Photon.CurrentMaterial.conc_a
                * 1e-3
                * self.Photon.CurrentMaterial.emax_A
                * np.log(10)
                * np.maximum(
                    1e-99,
                    self.Photon.CurrentMaterial.sp_abs_A[
                        RTU.FindClosestListElement(
                            self.Photon.CurrentMaterial.sp_abs_A[:, 0],
                            self.Photon.Wavelength
                        ),
                        1
                    ]
                )
            )
            b_d = (
                self.Photon.CurrentMaterial.conc_b
                * 1e-3
                * self.Photon.CurrentMaterial.emax_D
                * np.log(10)
                * np.maximum(
                    1e-99,
                    self.Photon.CurrentMaterial.sp_abs_D[
                        RTU.FindClosestListElement(
                            self.Photon.CurrentMaterial.sp_abs_D[:, 0],
                            self.Photon.Wavelength
                        ),
                        1
                    ]
                )
            )

            if self.Photon.CurrentMaterial.IsAligned:
                # Rotate the electric field vector into alignment space
                tm = np.roll(RTU.MakeVectAxis(self.Photon.CurrentMaterial.AlignVect), 2, axis=1)
                n_E = np.dot(self.Photon.ElectricFieldVector, tm)
                n_E /= np.linalg.norm(n_E)
                self.Photon.Psi = self.LambertBeerFullIntegral(n_E)
                print("Overlap Integral", self.Photon.Psi)
                b_a *= self.Photon.Psi * 3.0

            total_abs_coeff = b_a + b_d
            scale_param = 1.0 / total_abs_coeff if total_abs_coeff != 0 else 1e12

            # Distance left to absorption follows an exponential distribution
            self.Photon.DistanceLeftToAbsorption = (
                1 / self.Photon.ScaleFactor
            ) * np.random.exponential(scale=scale_param)
        else:
            self.Photon.DistanceLeftToAbsorption = 0.0

        print("NEW FREE LENGTH", self.Photon.DistanceLeftToAbsorption)
