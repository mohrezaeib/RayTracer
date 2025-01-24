import copy as cp
import numpy as np
import RTUtils as RTU
import scipy.special as special
import MonteCarloFRET as MCF


class Rectangle():
    def __init__(self,Height,Width,SupportVector,DirectionVector,Photon,NormalVector=[0,0,0]):
        self.Height = Height
        self.Width = Width
        self.DirectionVector = DirectionVector
        self.DirectionVector /= np.linalg.norm(self.DirectionVector)
        self.SupportVector = SupportVector
        self.Photon = Photon

        if NormalVector == [0,0,0]:
            self.e1 = [-self.DirectionVector[1], self.DirectionVector[0], 0]
            self.e2 = [0, -self.DirectionVector[2], self.DirectionVector[1]]
            ealt = [-self.DirectionVector[2], 0, self.DirectionVector[0]]
            if self.e1 == [0, 0, 0]:
                self.e1 = ealt
            elif self.e2 == [0, 0, 0]:
                self.e2 = ealt
            elif self.e1 == self.e2 or np.all(np.array(self.e1) == -np.array(self.e2)):
                self.e2 = ealt
        else:
            self.e1 = NormalVector
            self.e2 = np.cross(self.e1, self.DirectionVector)


    def FindIntersections(self):
        if np.dot(self.Photon.Direction, self.DirectionVector) != 0:
            BaseVect = (np.dot((self.SupportVector-self.Photon.Location), self.DirectionVector)/(np.dot(self.Photon.Direction, self.DirectionVector)))*self.Photon.Direction + self.Photon.Location
            #Check if there is an intersection with the current,  limited plane
            #Calculate components of distance between middle of plane and the intersection points and reject points accordingly
            d = BaseVect - self.SupportVector
            #s = np.dot(self.DirectionVector, d)

            self.e1 /= np.linalg.norm(self.e1)
            self.e2 /= np.linalg.norm(self.e2)

            t1 = np.dot(self.e1, d)
            t2 = np.dot(self.e2, d)

            if self.Height < np.absolute(t1) or self.Width < np.absolute(t2):
                # print "Bypassed"
                BaseVect = [np.nan, np.nan, np.nan]

            return BaseVect
        else:
            return [np.nan, np.nan, np.nan]



    def OperationPlane(self):
        return [self.SupportVector, self.DirectionVector]


class Circle():
    def __init__(self, Radius, SupportVector, DirectionVector,Photon):
        self.Radius = Radius
        self.DirectionVector = DirectionVector
        self.DirectionVector /= np.linalg.norm(self.DirectionVector)
        self.SupportVector = SupportVector
        self.Photon = Photon

    def FindIntersections(self):

        if np.dot(self.Photon.Direction, self.DirectionVector) != 0:
            '''Circles'''
            BaseVect = (np.dot((self.SupportVector-self.Photon.Location), self.DirectionVector)/(np.dot(self.Photon.Direction, self.DirectionVector)))*self.Photon.Direction + self.Photon.Location
            # print "#######################", np.linalg.norm(BaseVect - self.SupportVector), BaseVect,  self.SupportVector
            '''Reject Points that are not within the circle'''
            if self.Radius < np.linalg.norm(BaseVect - self.SupportVector):
                BaseVect = [np.nan, np.nan, np.nan]

            return BaseVect
        else:
            return [np.nan, np.nan, np.nan]




    def OperationPlane(self):
        return [self.SupportVector, self.DirectionVector]



class Cylinder():
    def __init__(self, Radius, SupportVector, DirectionVector,Photon):
        self.Radius = Radius
        self.DirectionVector = DirectionVector
        self.SupportVector = SupportVector
        self.Photon = Photon

    def FindIntersections(self):
        '''Deal with Cylinders'''
        #Find Intersections
        p1, p2 = self.IntersectLineCylinder()

        if not self.Photon.WasElement:
            nzcomp = np.argwhere(np.array(self.Photon.Direction) != 0)[0][0]
            if (p1[nzcomp] - self.Photon.Location[nzcomp])/self.Photon.Direction[nzcomp] < 0:
                p1 = [np.nan, np.nan, np.nan]
            if (p2[nzcomp] - self.Photon.Location[nzcomp])/self.Photon.Direction[nzcomp] < 0:
                p2 = [np.nan, np.nan, np.nan]

        #print self.Photon.Location,p1,p2
        #print np.all(p1 != [np.nan, np.nan, np.nan]),np.all(p2 != [np.nan, np.nan, np.nan])


        #Reject intersections too far off
        #find first intersection with cylinder
        if np.all(p1 != [np.nan, np.nan, np.nan]) and np.all(p2 != [np.nan, np.nan, np.nan]):
            #If both intersection points are not the starting point take the closest point,  else take the one that is not the starting point
            if np.linalg.norm(p1-self.Photon.Location) > 1e-7 and np.linalg.norm(p2-self.Photon.Location) > 1e-7:
                if np.linalg.norm(p1-self.Photon.Location) < np.linalg.norm(p2-self.Photon.Location):
                    BaseVect = p1
                else:
                    BaseVect = p2

            if np.linalg.norm(p1-self.Photon.Location) > 1e-7 and np.linalg.norm(p2-self.Photon.Location) <= 1e-7:
                BaseVect = p1
            if np.linalg.norm(p1-self.Photon.Location) <= 1e-7 and np.linalg.norm(p2-self.Photon.Location) > 1e-7:
                BaseVect = p2

        elif np.all(p1 != [np.nan, np.nan, np.nan]) and np.all(p2 == [np.nan, np.nan, np.nan]):
            BaseVect = p1
        elif np.all(p2 != [np.nan, np.nan, np.nan]) and np.all(p1 == [np.nan, np.nan, np.nan]):
            BaseVect = p2
        else:
            BaseVect = [np.nan, np.nan, np.nan]

        #print BaseVect
        return BaseVect

    def FindNormalCylinder(self):
        '''Find nromal of cylinder at a point,  n is the element number of the cylinder while p is an arbitary point
        returns the normal
        '''
        a = self.SupportVector
        d = self.DirectionVector - self.SupportVector
        d /= np.linalg.norm(d)
        P = self.Photon.PointCy

        lfp = a + (np.dot(P-a, d)) *d
        n_p = lfp-P
        n_p /= np.linalg.norm(n_p)

        # plpt = point + n_p
        # self.ax.plot([point[0], plpt[0]], [point[1], plpt[1]], [point[2], plpt[2]], color='blue')
        return n_p


    def IntersectLineCylinder(self):
        '''Find intersection point of Ray and Cylinder'''
        r = self.Radius
        #find the vector from point p1 to point p2
        a = self.SupportVector - self.DirectionVector
        abs_a = np.linalg.norm(a)
        a /= abs_a
        #Transformation matrix for coordinates,  make cylinder axis the new x axis
        tm = RTU.MakeVectAxis(a)
        tm = np.transpose(tm)
        tmi = np.linalg.inv(tm)

        #t_cy_line = np.dot(tm,  a)
        t_Supp = np.dot(tm, self.Photon.Location)
        t_Direc = np.dot(tm, self.Photon.Direction)
        #t_End = np.dot(tm, self.Photon.Location+self.Photon.Direction)
        Cy_p1 = np.dot(tm, self.SupportVector)
        Cy_p2 = np.dot(tm, self.DirectionVector)

        #Find intersection of line with circle
        a1 = t_Direc[1]**2 + t_Direc[2]**2
        a2 = 2*(t_Direc[1]*(t_Supp[1]-Cy_p1[1]) + t_Direc[2]*(t_Supp[2]-Cy_p1[2]))
        a3 = (t_Supp[1]-Cy_p1[1])**2 + (t_Supp[2]-Cy_p1[2])**2 - r**2
        dsc = a2 ** 2 - 4*a1*a3
        # print dsc

        if a1 == 0.0:
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]

        #Check if there is an intersection
        if dsc >= 0:
            t1 = (-a2 + np.sqrt(dsc))/(2*a1)
            t2 = (-a2 - np.sqrt(dsc))/(2*a1)
            v1 = t_Supp + t1*t_Direc
            v2 = t_Supp + t2*t_Direc
            p1 = np.dot(tmi, t_Supp + t1*t_Direc)
            p2 = np.dot(tmi, t_Supp + t2*t_Direc)
            # print "TransformedPoints", Cy_p1[0], Cy_p2[0], "v1v2", v1, v2
            if v1[0] < np.minimum(Cy_p1[0], Cy_p2[0]) or v1[0] > np.maximum(Cy_p1[0], Cy_p2[0]):
                p1 = [np.nan, np.nan, np.nan]
            if v2[0] < np.minimum(Cy_p1[0], Cy_p2[0]) or v2[0] > np.maximum(Cy_p1[0], Cy_p2[0]):
                p2 = [np.nan, np.nan, np.nan]

            return p1, p2
            # self.ax.scatter(p1[0], p1[1], p1[2])
            # self.ax.scatter(p2[0], p2[1], p2[2], marker = 1)
        else:
            # print "no intersect with", n
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]


    def OperationPlane(self):
        n = self.FindNormalCylinder()
        return [self.Photon.PointCy, n]


class Triangle():
    def __init__(self, Point1,Point2,Point3, Photon):
        self.TrianglePoint1 = Point1
        self.TrianglePoint2 = Point2
        self.TrianglePoint3 = Point3
        self.Photon = Photon
        self.SupportVector, self.DirectionVector = RTU.MakeTriangle(Point1, Point2, Point3)

    def FindIntersections(self):
        #Triangles
        if np.dot(self.Photon.Direction, self.DirectionVector) != 0:
            BaseVect = (np.dot((self.SupportVector-self.Photon.Location), self.DirectionVector)/(np.dot(self.Photon.Direction, self.DirectionVector)))*self.Photon.Direction + self.Photon.Location
            #Reject if point is not in triangle
            if not self.CheckPointInTriangle(BaseVect):
                BaseVect = [np.nan, np.nan, np.nan]

            return BaseVect
        else:
            return [np.nan, np.nan, np.nan]

    def CheckPointInTriangle(self, point):
        '''checks if point is inside triangle,  partly stolen from http://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle'''
        trng = [self.TrianglePoint1, self.TrianglePoint2, self.TrianglePoint3]
        #plane = RTU.MakeTriangle(trng[0], trng[1], trng[2])
        plane = self.OperationPlane()
        tm = RTU.MakeVectAxis(plane[1])
        tm = np.transpose(tm)
        for i in range(3):
            trng[i] = np.dot(tm, trng[i])
        pt = cp.deepcopy(np.dot(tm, point)[1:])

        b1 = self.SclPrdSign(pt, trng[0][1:], trng[1][1:]) < 0.0
        b2 = self.SclPrdSign(pt, trng[1][1:], trng[2][1:]) < 0.0
        b3 = self.SclPrdSign(pt, trng[2][1:], trng[0][1:]) < 0.0
        # print "TRIANGLEBOOLS", b1, b2, b3
        return b1 == b2 == b3

    def SclPrdSign(self, p1, p2, p3):
        '''Find sign of dot product,  assist Checkpoint in triangle'''
        return (p1[0]-p3[0]) * (p2[1]-p3[1]) - (p2[0]-p3[0])*(p1[1]-p3[1])


    def OperationPlane(self):
        return [self.SupportVector, self.DirectionVector]


#
