import copy as cp
import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF


class Rectangle():
    def __init__(self, Height, Width, SupportVector, DirectionVector, Photon, NormalVector=[0, 0, 0]):
        self.Height = Height
        self.Width = Width
        self.DirectionVector = np.array(DirectionVector, dtype=float)
        # Normalize direction
        self.DirectionVector /= np.linalg.norm(self.DirectionVector)

        self.SupportVector = np.array(SupportVector, dtype=float)
        self.Photon = Photon

        # Decide how to build e1, e2
        if NormalVector == [0, 0, 0]:
            # Build orthogonal vectors e1 and e2
            # so that e1, e2, DirectionVector are mutually perpendicular
            self.e1 = [
                -self.DirectionVector[1],
                self.DirectionVector[0],
                0
            ]
            self.e2 = [
                0,
                -self.DirectionVector[2],
                self.DirectionVector[1]
            ]
            ealt = [
                -self.DirectionVector[2],
                0,
                self.DirectionVector[0]
            ]
            # Check for degenerate cases
            if self.e1 == [0, 0, 0]:
                self.e1 = ealt
            elif self.e2 == [0, 0, 0]:
                self.e2 = ealt
            elif (
                self.e1 == self.e2
                or np.all(np.array(self.e1) == -np.array(self.e2))
            ):
                self.e2 = ealt
        else:
            # If a normal vector is given, use it directly
            self.e1 = np.array(NormalVector, dtype=float)
            # e2 is the cross product with the main direction
            self.e2 = np.cross(self.e1, self.DirectionVector)

        # Convert e1, e2 to numpy arrays for future in-place operations
        self.e1 = np.array(self.e1, dtype=float)
        self.e2 = np.array(self.e2, dtype=float)

    def FindIntersections(self):
        """
        Find intersection of the Photon's ray with this (limited) rectangle.
        """
        denom = np.dot(self.Photon.Direction, self.DirectionVector)
        if denom == 0:
            # Ray is parallel to rectangle plane
            return [np.nan, np.nan, np.nan]

        # Intersection with the infinite plane
        BaseVect = (
            np.dot((self.SupportVector - self.Photon.Location), self.DirectionVector)
            / denom
        ) * self.Photon.Direction + self.Photon.Location

        # Check if intersection lies within rectangle boundaries
        d = BaseVect - self.SupportVector

        # Normalize e1 and e2 here to project properly
        if np.linalg.norm(self.e1) != 0:
            self.e1 /= np.linalg.norm(self.e1)
        if np.linalg.norm(self.e2) != 0:
            self.e2 /= np.linalg.norm(self.e2)

        t1 = np.dot(self.e1, d)
        t2 = np.dot(self.e2, d)

        if (abs(t1) > self.Height) or (abs(t2) > self.Width):
            # Intersection falls outside rectangle boundary
            BaseVect = [np.nan, np.nan, np.nan]

        return BaseVect

    def OperationPlane(self):
        """
        Return [support point on plane, normal vector of plane].
        """
        return [self.SupportVector, self.DirectionVector]


class Circle():
    def __init__(self, Radius, SupportVector, DirectionVector, Photon):
        self.Radius = Radius
        self.DirectionVector = np.array(DirectionVector, dtype=float)
        # Normalize direction
        self.DirectionVector /= np.linalg.norm(self.DirectionVector)

        self.SupportVector = np.array(SupportVector, dtype=float)
        self.Photon = Photon

    def FindIntersections(self):
        """
        Find intersection of Photon's ray with this circle (disk).
        """
        denom = np.dot(self.Photon.Direction, self.DirectionVector)
        if denom == 0:
            # Ray is parallel to circle's plane
            return [np.nan, np.nan, np.nan]

        BaseVect = (
            np.dot((self.SupportVector - self.Photon.Location), self.DirectionVector)
            / denom
        ) * self.Photon.Direction + self.Photon.Location

        # Check if intersection point is within the radius
        dist = np.linalg.norm(BaseVect - self.SupportVector)
        if dist > self.Radius:
            BaseVect = [np.nan, np.nan, np.nan]

        return BaseVect

    def OperationPlane(self):
        """
        Return [support point on plane, normal vector of plane].
        """
        return [self.SupportVector, self.DirectionVector]


class Cylinder():
    def __init__(self, Radius, SupportVector, DirectionVector, Photon):
        self.Radius = Radius
        self.SupportVector = np.array(SupportVector, dtype=float)
        self.DirectionVector = np.array(DirectionVector, dtype=float)
        self.Photon = Photon

    def FindIntersections(self):
        """
        Deal with cylinders: find intersection of the Photon's ray with a finite cylinder.
        Currently, the code treats the cylinder as a segment aligned between
        SupportVector and DirectionVector with a certain radius.
        """
        p1, p2 = self.IntersectLineCylinder()

        # If the photon has not yet hit an element, ensure direction is forward
        # (no negative t solutions)
        if not self.Photon.WasElement:
            # Find first non-zero component to check sign
            nzcomp = np.argwhere(np.array(self.Photon.Direction) != 0)[0][0]
            if (
                (p1[nzcomp] - self.Photon.Location[nzcomp])
                / self.Photon.Direction[nzcomp]
            ) < 0:
                p1 = [np.nan, np.nan, np.nan]

            if (
                (p2[nzcomp] - self.Photon.Location[nzcomp])
                / self.Photon.Direction[nzcomp]
            ) < 0:
                p2 = [np.nan, np.nan, np.nan]

        # Choose the correct intersection
        #   - If both p1 and p2 are valid and not the starting location,
        #     choose the one closest to the photon's location.
        #   - If only one is valid, return it.
        #   - Else return [np.nan, np.nan, np.nan].
        BaseVect = [np.nan, np.nan, np.nan]

        valid_p1 = np.all(np.array(p1) != np.array([np.nan, np.nan, np.nan]))
        valid_p2 = np.all(np.array(p2) != np.array([np.nan, np.nan, np.nan]))

        if valid_p1 and valid_p2:
            dist_p1 = np.linalg.norm(p1 - self.Photon.Location)
            dist_p2 = np.linalg.norm(p2 - self.Photon.Location)

            # If both intersection points differ from the starting point
            if dist_p1 > 1e-7 and dist_p2 > 1e-7:
                BaseVect = p1 if dist_p1 < dist_p2 else p2
            elif dist_p1 > 1e-7 and dist_p2 <= 1e-7:
                BaseVect = p1
            elif dist_p1 <= 1e-7 and dist_p2 > 1e-7:
                BaseVect = p2

        elif valid_p1 and not valid_p2:
            BaseVect = p1
        elif valid_p2 and not valid_p1:
            BaseVect = p2

        return BaseVect

    def FindNormalCylinder(self):
        """
        Find normal of cylinder at the photon's current intersection point.
        The cylinder is defined by a support vector and a direction vector.
        """
        a = self.SupportVector
        d = self.DirectionVector - self.SupportVector
        d /= np.linalg.norm(d)

        P = self.Photon.PointCy
        # Project P onto the cylinder's axis
        lfp = a + (np.dot(P - a, d)) * d
        # Normal is the direction from P to the axis
        n_p = lfp - P
        n_p /= np.linalg.norm(n_p)
        return n_p

    def IntersectLineCylinder(self):
        """
        Find intersection points of a ray with a finite cylinder of radius R.
        Cylinder axis: from self.SupportVector to self.DirectionVector.
        """
        r = self.Radius

        # Vector along the cylinder axis
        a = self.SupportVector - self.DirectionVector
        abs_a = np.linalg.norm(a)
        if abs_a == 0:
            # Degenerate cylinder
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]
        a /= abs_a

        # Build transformation (rotation) so that
        # the cylinder axis aligns with the new x-axis in transformed coordinates
        tm = RTU.MakeVectAxis(a)
        tm = np.transpose(tm)
        tmi = np.linalg.inv(tm)

        # Transform photon location and direction
        t_Supp = np.dot(tm, self.Photon.Location)
        t_Direc = np.dot(tm, self.Photon.Direction)

        # Transform cylinder endpoints
        Cy_p1 = np.dot(tm, self.SupportVector)
        Cy_p2 = np.dot(tm, self.DirectionVector)

        # Coefficients for intersection with circle in y-z plane (transformed space)
        a1 = t_Direc[1] ** 2 + t_Direc[2] ** 2
        a2 = 2 * (
            t_Direc[1] * (t_Supp[1] - Cy_p1[1])
            + t_Direc[2] * (t_Supp[2] - Cy_p1[2])
        )
        a3 = (
            (t_Supp[1] - Cy_p1[1]) ** 2
            + (t_Supp[2] - Cy_p1[2]) ** 2
            - r ** 2
        )

        # If a1 is zero, direction is purely x-axis in transformed space
        if a1 == 0.0:
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]

        dsc = a2 ** 2 - 4 * a1 * a3
        if dsc < 0:
            # No real intersection
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]

        # Two solutions
        sqrt_dsc = np.sqrt(dsc)
        t1 = (-a2 + sqrt_dsc) / (2 * a1)
        t2 = (-a2 - sqrt_dsc) / (2 * a1)

        v1 = t_Supp + t1 * t_Direc
        v2 = t_Supp + t2 * t_Direc

        p1 = np.dot(tmi, v1)
        p2 = np.dot(tmi, v2)

        # Clip by cylinder length in the x-axis (transformed space)
        min_x = min(Cy_p1[0], Cy_p2[0])
        max_x = max(Cy_p1[0], Cy_p2[0])

        if not (min_x <= v1[0] <= max_x):
            p1 = [np.nan, np.nan, np.nan]
        if not (min_x <= v2[0] <= max_x):
            p2 = [np.nan, np.nan, np.nan]

        return p1, p2

    def OperationPlane(self):
        """
        Return [current photon intersection point on cylinder,
                the local normal at that intersection].
        """
        n = self.FindNormalCylinder()
        return [self.Photon.PointCy, n]


class Triangle():
    def __init__(self, Point1, Point2, Point3, Photon):
        self.TrianglePoint1 = np.array(Point1, dtype=float)
        self.TrianglePoint2 = np.array(Point2, dtype=float)
        self.TrianglePoint3 = np.array(Point3, dtype=float)
        self.Photon = Photon

        # Build plane data from the 3 points
        self.SupportVector, self.DirectionVector = RTU.MakeTriangle(
            self.TrianglePoint1,
            self.TrianglePoint2,
            self.TrianglePoint3
        )

    def FindIntersections(self):
        """
        Find intersection of the Photon's ray with the triangle (treated as
        a finite area in a plane).
        """
        denom = np.dot(self.Photon.Direction, self.DirectionVector)
        if denom == 0:
            return [np.nan, np.nan, np.nan]

        BaseVect = (
            np.dot(
                (self.SupportVector - self.Photon.Location),
                self.DirectionVector
            )
            / denom
        ) * self.Photon.Direction + self.Photon.Location

        # Reject if point is not inside the triangle boundaries
        if not self.CheckPointInTriangle(BaseVect):
            BaseVect = [np.nan, np.nan, np.nan]

        return BaseVect

    def CheckPointInTriangle(self, point):
        """
        Checks if a given point lies inside this triangle (2D test in
        a local coordinate system).
        Borrowed logic from:
        https://stackoverflow.com/questions/2049582/
        """
        trng = [
            self.TrianglePoint1,
            self.TrianglePoint2,
            self.TrianglePoint3
        ]

        plane = self.OperationPlane()
        tm = RTU.MakeVectAxis(plane[1])
        tm = np.transpose(tm)

        # Transform the triangle corners and the point
        trng_trans = []
        for i in range(3):
            trng_trans.append(np.dot(tm, trng[i]))

        pt = np.dot(tm, point)
        pt_2D = pt[1:]  # ignore transformed x (axis along plane normal)

        # Evaluate sign of areas
        b1 = self.SclPrdSign(pt_2D, trng_trans[0][1:], trng_trans[1][1:]) < 0.0
        b2 = self.SclPrdSign(pt_2D, trng_trans[1][1:], trng_trans[2][1:]) < 0.0
        b3 = self.SclPrdSign(pt_2D, trng_trans[2][1:], trng_trans[0][1:]) < 0.0

        # If all signs are the same, point lies inside
        return (b1 == b2 == b3)

    def SclPrdSign(self, p1, p2, p3):
        """
        Find the sign of the cross product (2D) between the vectors
        (p1 - p3) and (p2 - p3). Assists in CheckPointInTriangle.
        """
        return (
            (p1[0] - p3[0]) * (p2[1] - p3[1])
            - (p2[0] - p3[0]) * (p1[1] - p3[1])
        )

    def OperationPlane(self):
        """
        Return [support vector, normal of plane].
        """
        return [self.SupportVector, self.DirectionVector]
