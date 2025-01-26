import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import Geometries


class DrawElements():
    def __init__(self, OpticalElements):
        # Define Plot axis
        plt.ion()
        self.fig = plt.figure(figsize=plt.figaspect(1.0))
        self.ax = self.fig.add_subplot(111, projection='3d')

        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')

        self.OpticalElements = OpticalElements

    def DrawCylinder(self, Element, Type):
        """
        Draw a Cylinder from p1 to p2 with a radius r
        Adapted from here:
        http://stackoverflow.com/questions/32317247/how-to-draw-a-cylinder-using-matplotlib-along-length-of-point-x1-y1-and-x2-y2
        """
        r = Element.Radius
        # Find the vector from point p1 to point p2
        a = Element.DirectionVector - Element.SupportVector
        abs_a = np.linalg.norm(a)
        a /= abs_a

        # Make some vector not in the same direction as v
        not_a = np.array([1, 0, 0])
        if (a == not_a).all():
            not_a = np.array([0, 1, 0])
        # Make vector perpendicular to v
        n1 = np.cross(a, not_a)
        n1 /= np.linalg.norm(n1)
        # Make unit vector perpendicular to v and n1
        n2 = np.cross(a, n1)

        t = np.linspace(0, abs_a, 10)
        theta = np.linspace(0, 2 * np.pi, 10)
        t, theta = np.meshgrid(t, theta)

        # Generate coordinates for surface
        X, Y, Z = [
            Element.SupportVector[i]
            + a[i] * t
            + r * np.sin(theta) * n1[i]
            + r * np.cos(theta) * n2[i]
            for i in [0, 1, 2]
        ]
        self.ax.plot_surface(X, Y, Z, color=Type.color, alpha=Type.alpha)

    def DrawCircle(self, Element, Type):
        """
        Draw Circle with radius r
        """
        r = Element.Radius
        n = Element.DirectionVector
        # Make some vector not in the same direction as v
        not_n = np.array([1, 0, 0])
        if (n == not_n).all():
            not_n = np.array([0, 1, 0])
        # Make vector perpendicular to v
        n1 = np.cross(n, not_n)
        n1 /= np.linalg.norm(n1)
        # Make unit vector perpendicular to v and n1
        n2 = np.cross(Element.DirectionVector, n1)

        theta = np.linspace(0, 2 * np.pi, 10)
        X, Y, Z = [
            Element.SupportVector[i]
            + r * np.sin(theta) * n1[i]
            + r * np.cos(theta) * n2[i]
            for i in [0, 1, 2]
        ]

        # Convert zip object to list for Poly3DCollection in Python 3
        verts = [list(zip(X, Y, Z))]

        Collection = Poly3DCollection(verts, linewidths=1, alpha=Type.alpha)
        FaceColor = Type.color
        Collection.set_facecolor(FaceColor)
        Collection.set_edgecolor("black")
        self.ax.add_collection3d(Collection)

    def DrawElementsPoly(self):
        """
        Draw Stuff in 3D using polygons
        """
        self.ax.cla()
        i = 0
        for el in self.OpticalElements:
            # If element is a plane
            if isinstance(el.Geometry, Geometries.Rectangle):
                e1 = el.Geometry.e1 / np.linalg.norm(el.Geometry.e1)
                e2 = el.Geometry.e2 / np.linalg.norm(el.Geometry.e2)

                a = el.Geometry.SupportVector + e1 * el.Geometry.Height + e2 * el.Geometry.Width
                b = el.Geometry.SupportVector + e1 * el.Geometry.Height - e2 * el.Geometry.Width
                c = el.Geometry.SupportVector - e1 * el.Geometry.Height - e2 * el.Geometry.Width
                d = el.Geometry.SupportVector - e1 * el.Geometry.Height + e2 * el.Geometry.Width

                x = [a[0], b[0], c[0], d[0]]
                y = [a[1], b[1], c[1], d[1]]
                z = [a[2], b[2], c[2], d[2]]

                verts = [list(zip(x, y, z))]
                Collection = Poly3DCollection(verts, linewidths=1, alpha=el.Type.alpha)
                FaceColor = el.Type.color
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)

            if isinstance(el.Geometry, Geometries.Circle):
                self.DrawCircle(el.Geometry, el.Type)

            if isinstance(el.Geometry, Geometries.Cylinder):
                self.DrawCylinder(el.Geometry, el.Type)

            if isinstance(el.Geometry, Geometries.Triangle):
                print(el.Geometry.TrianglePoint1, el.Geometry.TrianglePoint2, el.Geometry.TrianglePoint3)
                x = [
                    el.Geometry.TrianglePoint1[0],
                    el.Geometry.TrianglePoint2[0],
                    el.Geometry.TrianglePoint3[0],
                ]
                y = [
                    el.Geometry.TrianglePoint1[1],
                    el.Geometry.TrianglePoint2[1],
                    el.Geometry.TrianglePoint3[1],
                ]
                z = [
                    el.Geometry.TrianglePoint1[2],
                    el.Geometry.TrianglePoint2[2],
                    el.Geometry.TrianglePoint3[2],
                ]
                verts = [list(zip(x, y, z))]
                Collection = Poly3DCollection(verts, linewidths=1, alpha=el.Type.alpha)
                FaceColor = el.Type.color
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)

            i += 1

        self.ax.set_axis_off()
        self.ax.grid(False)
        self.ax.set_autoscale_on(True)
        self.ax.set_xlim([-1, 1])
        self.ax.set_ylim([-1, 1])
        self.ax.set_zlim([-1, 1])

    def DrawLine(self, point1, point2, LineColor):
        """
        Draw a Line from one point to another using the specified color
        """
        self.ax.plot(
            [point1[0], point2[0]],
            [point1[1], point2[1]],
            [point1[2], point2[2]],
            color=LineColor
        )

    def DrawPoint(self, point):
        """
        Draw a single point
        """
        self.ax.scatter(point[0], point[1], point[2])
