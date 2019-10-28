import mpl_toolkits.mplot3d as a3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy as sp
import math
from random import randint
from point_base import point
from vector_base import vector
from basis_base import basis
from quaternion_base import quaternion
from ellipsoid_base import ellipsoid_base

class EFS(ellipsoid_base):
    """ Ellipsoidal fibonacci spiral.

    Create an ellipsoidal distribution of N points
    using the Fibonacci Sprial.
    """

    def __init__(self, rx=1, ry=1, rz=1, n_points=100, precision=3):
        super().__init__(rx,ry,rz)
        self.dec_prec = precision
        self.points = self._spiral(n_points)

    def _spiral(self, n_points):
        """ Create a 'regular' ellipsoidal ichosaedron.

        Use the radii defined in the constructor
        to shape the spiral to the desired ellipsoid.
        """

        # Create spiral on unit sphere
        dlong = math.pi*(3 - math.sqrt(5))
        dz = 2.0 / n_points
        long = 0
        z = 1 - dz/2
        points = []
        for i in range(n_points):
            r = math.sqrt(1-z*z)
            points.append(point(math.cos(long)*r, math.sin(long)*r, z))
            z -=  dz
            long += dlong
        # Adjust to ellipsoid radii
        for p in points:
            p.x *= self.rx
            p.y *= self.ry
            p.z *= self.rz
        # Return
        return points

    def _plot_points(self, x, y, z, crop_z=-10, plane=True):
        """ Plot a set of 3D points."""

        # Create figure with axes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Filter points lower than crop distance
        xc = [x[i] for i in range(len(x)) if z[i] >=crop_z]
        yc = [y[i] for i in range(len(y)) if z[i] >=crop_z]
        zc = [z[i] for i in range(len(z)) if z[i] >=crop_z]
        print("# points: " + str(len(xc)))
        # Add points to graph
        ax.scatter(xc, yc, zc, c='r', marker='o')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        # Add visual plane
        if plane:
            poly = [[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]]
            face = a3.art3d.Poly3DCollection([poly])
            face.set_alpha(0.1)
            face.set_color(colors.rgb2hex(sp.rand(3)))
            face.set_edgecolor('k')
            ax.add_collection3d(face)
        # Scale
        ax.auto_scale_xyz([-1, 1], [-1, 1], [-1, 1])
        # Call matplotlib plot()
        plt.show()

    def _plot_bases(self, bases, crop_z=-10, plane=True, vec_scale=5):
        """ Plot a set of 3D points with coordinate bases."""

        # Create figure with axes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Filter bases with z lower than crop distance
        filtered_bases = [i for i in bases if i.origin[2] >= crop_z]
        xc = [i.origin[0] for i in filtered_bases]
        yc = [i.origin[1] for i in filtered_bases]
        zc = [i.origin[2] for i in filtered_bases]
        vxc = [i.vx for i in filtered_bases]
        vyc = [i.vy for i in filtered_bases]
        vzc = [i.vz for i in filtered_bases]
        print("# points: " + str(len(filtered_bases)))
        # Add points to graph
        ax.scatter(xc, yc, zc, c='r', marker='o')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        # Add vectors to graph
        ax.quiver(xc, yc, zc, [vxc[i][0]/vec_scale for i in range(len(vxc))], [vxc[i][1]/vec_scale for i in range(len(vxc))], [vxc[i][2]/vec_scale for i in range(len(vxc))], color="red", arrow_length_ratio=0.5)
        ax.quiver(xc, yc, zc, [vyc[i][0]/vec_scale for i in range(len(vyc))], [vyc[i][1]/vec_scale for i in range(len(vyc))], [vyc[i][2]/vec_scale for i in range(len(vyc))], color="green", arrow_length_ratio=0.5)
        ax.quiver(xc, yc, zc, [vzc[i][0]/vec_scale for i in range(len(vzc))], [vzc[i][1]/vec_scale for i in range(len(vzc))], [vzc[i][2]/vec_scale for i in range(len(vzc))], color="blue", arrow_length_ratio=0.5)
        # Add visual plane
        if plane:
            poly = [[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]]
            face = a3.art3d.Poly3DCollection([poly])
            face.set_alpha(0.1)
            face.set_color(colors.rgb2hex(sp.rand(3)))
            face.set_edgecolor('k')
            ax.add_collection3d(face)
        # Scale
        ax.auto_scale_xyz([-1, 1], [-1, 1], [-1, 1])
        # Call matplotlib plot()
        plt.show()

    def _xyz_list(self, vs):
        """ Transform a nx3 array into 3 nx1 arrays."""
        
        x = [vi[0] for vi in vs]
        y = [vi[1] for vi in vs]
        z = [vi[2] for vi in vs]
        return x,y,z

    def tesselate(self, n_times=1):
        px, py, pz = self._xyz_list(self.points)
        self._plot_points(px, py, pz)
        vs = [basis(i, self._point_normal(i)) for i in self.points]
        self._plot_bases(vs, 0.15)
        for i in range(len(vs)):
            a = quaternion()
            b = a.from_matrix(vs[i].matrix())

ob = EFS(1,1,3)
ob.tesselate(2)

#Objectives:
#   -From points to robtargets
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
