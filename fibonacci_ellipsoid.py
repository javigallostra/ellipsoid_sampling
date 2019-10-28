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

    def __init__(self, rx=1, ry=1, rz=1, n_points=100):
        super().__init__(rx,ry,rz)
        self.points = self._spiral(n_points)
        self.basis = [basis(p, self._point_normal(p)) for p in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh

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

    def plot(self, points=True, basis=False, faces=False, crop_z=None):
        """ Plot the sampled ellipsoid."""

        # Create axes
        ax = a3.Axes3D(plt.figure())
        # Plot points if desired
        if points:
            # Crop
            if crop_z != None:
                points_crop = [p for p in self.points if p[2] >=crop_z]
            else:
                points_crop = self.points
            # Prepare
            px, py, pz = self._xyz_list(points_crop)
            # Plot
            ax.scatter(px, py, pz, c='r', marker='o')
        # Plot basis if desired
        if basis:
            # Crop
            if crop_z != None:
                basis_crop = [b for b in self.basis if b.origin[2] >= crop_z]
            else:
                basis_crop = self.basis
            # Prepare
            px = [b.origin[0] for b in basis_crop]
            py = [b.origin[1] for b in basis_crop]
            pz = [b.origin[2] for b in basis_crop]
            vx = [b.vx for b in basis_crop]
            vy = [b.vy for b in basis_crop]
            vz = [b.vz for b in basis_crop]
            # Plot
            scale = 5
            ax.quiver(px, py, pz, [v[0]/scale for v in vx], [v[1]/scale for v in vx], [v[2]/scale for v in vx], color="red", arrow_length_ratio=0.5)
            ax.quiver(px, py, pz, [v[0]/scale for v in vy], [v[1]/scale for v in vy], [v[2]/scale for v in vy], color="green", arrow_length_ratio=0.5)
            ax.quiver(px, py, pz, [v[0]/scale for v in vz], [v[1]/scale for v in vz], [v[2]/scale for v in vz], color="blue", arrow_length_ratio=0.5)
        # Plot faces if desired
        if faces:
            # Plot each face
            for f in self.faces:
                poly=[self.points[f[0]], self.points[f[1]], self.points[f[2]]]
                # Crop
                if crop_z != None:
                    if not(poly[0][2]>=crop_z and poly[1][2]>=crop_z and poly[2][2]>=crop_z):
                        continue
                # Plot
                face = a3.art3d.Poly3DCollection([poly])
                face.set_alpha(1)
                face.set_color(colors.rgb2hex(sp.rand(3)))
                face.set_edgecolor('k')
                ax.add_collection3d(face)
        # Add visual plane
        if crop_z != None:
            poly = [[self.rx,self.ry,crop_z],[-self.rx,self.ry,crop_z],[-self.rx,-self.ry,crop_z],[self.rx,-self.ry,crop_z]]
            face = a3.art3d.Poly3DCollection([poly])
            face.set_alpha(0.5)
            face.set_color(colors.rgb2hex(sp.rand(3)))
            face.set_edgecolor('k')
            ax.add_collection3d(face)
        # Add labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # Scale
        if crop_z != None:
            d = max([self.rx, self.ry, (self.rz - crop_z)/2])
            z_c = crop_z + (self.rz - crop_z)/2
            ax.auto_scale_xyz([-d, d], [-d, d], [z_c - d, z_c + d])
        else:
            d = max([self.rx, self.ry, self.rz])
            ax.auto_scale_xyz([-d, d], [-d, d], [-d, d])
        # Call matplotlib plot()
        plt.show()

    def _xyz_list(self, vs):
        """ Transform a nx3 array into 3 nx1 arrays."""
        
        x = [vi[0] for vi in vs]
        y = [vi[1] for vi in vs]
        z = [vi[2] for vi in vs]
        return x,y,z

ob = EFS(1,2,1)
ob.plot(False,True,False,0)

#Objectives:
#   -From points to robtargets
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
