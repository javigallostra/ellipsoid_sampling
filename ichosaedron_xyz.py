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

class EIT(ellipsoid_base):
    """ Ellipsoidal ichosaedron tesselation.

    Create an ellipsoidal ichosaedron and perform
    triangular subdivisions (1 to 4) to tesellate its surface.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        super().__init__(rx, ry, rz)
        self.points, self.faces = self._ichosaedron()

    def _ichosaedron(self):
        """ Create a 'regular' ellipsoidal ichosaedron.

        Use the radii defined in the constructor
        to shape the regular ichosaedron to the desired
        elipsoid.
        """
        
        # Create points in order (spherical coord)
        p0 = point(90, 0, None, "spherical")
        p1_5 = [point(26.57, 2*i*36, None, "spherical") for i in range(0,5)]
        p6_10 = [point(-26.57, (2*i+1)*36, None, "spherical") for i in range(0,5)]
        p11 = point(-90, 0, None, "spherical")
        # Create point list in order (cartesian coord)
        points = [self._spherical_to_cart(p0)]
        for i in range(0,5):
            points.append(self._spherical_to_cart(p1_5[i]))
            points.append(self._spherical_to_cart(p6_10[i]))
        points.append(self._spherical_to_cart(p11))
        # Create faces with the ordered vertices
        faces = []
        for i in range(0,4):
            faces.append([0,2*i+1,2*(i+1)+1])
            faces.append([11,2*(i+1),2*(i+2)])
        for i in range(0,8):
            faces.append([i+1,i+2,i+3])
        faces.append([0,1,9])
        faces.append([11,2,10])
        faces.append([9,10,1])
        faces.append([10,1,2])
        # Return
        return points, faces

    def _point_interpolation(self, p1, p2):
        """ Point interpolation function.

        Find the point equidistant to p1,p2
        lying on the ellipsoid.
        """

        # Find the middle point of p1-p2
        p = p1 + (p2 - p1)/2
        # Project the point onto the ellipsoid
        # Impossible case: 0,0,0
        if p.count(0) == 3:
            x = 0.0
            y = 0.0
            z = 0.0
        # Axial points
        elif p.count(0) == 2:
            x = self.rx * p.x / abs(p.x) if p.x else 0.0
            y = self.ry * p.y / abs(p.y) if p.y else 0.0
            z = self.rz * p.z / abs(p.z) if p.z else 0.0
        # Planar points
        elif p.count(0) == 1:
            #plane x=0
            if p.x == 0:
                x = 0.0
                #find y
                dydz = p.y / p.z
                y_num = self.ry * self.rz * abs(dydz)
                y_den = math.sqrt((self.rz * dydz)**2 + self.ry**2)
                y = y_num/y_den
                if p.y < 0: y = -y
                #compute z
                z = y / dydz
            #plane y=0
            elif p.y == 0:
                y = 0.0
                #find x
                dxdz = p.x / p.z
                x_num = self.rx * self.rz * abs(dxdz)
                x_den = math.sqrt((self.rz * dxdz)**2 + self.rx**2)
                x = x_num / x_den
                if p.x < 0: x = -x
                #compute z
                z = x / dxdz
            #plane z=0
            elif p.z == 0:
                z = 0.0
                #find x
                dxdy = p.x / p.y
                x_num = self.rx * self.ry * abs(dxdy)
                x_den = math.sqrt((self.ry * dxdy)**2 + self.rx**2)
                x = x_num / x_den
                if p.x < 0: x = -x
                #compute z
                y = x / dxdy
        # Quadrant points
        else:
            #find x
            dxdy = p.x / p.y
            dxdz = p.x / p.z
            x_num = self.rx * self.ry  *self.rz * abs(dxdy * dxdz)
            x_den = math.sqrt((dxdy * self.ry * dxdz * self.rz)**2 + (self.rx * dxdz * self.rz)**2 + (self.rx * dxdy * self.ry)**2)
            x = x_num / x_den
            if p.x < 0: x = -x
            #compute y, z
            y = x / dxdy
            z = x / dxdz
        # Return
        return point(x, y, z)

    def _subdivide(self, points, faces, times=1):
        """ Recursive triangular subdivision.

        Subdivide each surface triangle into 4
        a defined number of times to obtain a fine
        triangular tesselation.
        """

        # End case
        if times == 0:
            return (points, faces)
        # Recursive case
        else:
            new_faces = []
            new_points = []
            # Perform subdivision
            for i in range(len(faces)):
                #copy existing points
                p0 = points[faces[i][0]] * 1
                p2 = points[faces[i][1]] * 1
                p4 = points[faces[i][2]] * 1
                #interpolate new vertices
                p1 = self._point_interpolation(p0, p2)
                p3 = self._point_interpolation(p2, p4)
                p5 = self._point_interpolation(p4, p0)
                #add vertices if not repeated
                for new_point in [p0, p1, p2, p3, p4, p5]:
                    if new_point not in new_points:
                        new_points.append(new_point)
                #add faces
                new_faces.append([new_points.index(p0), new_points.index(p1), new_points.index(p5)])
                new_faces.append([new_points.index(p1), new_points.index(p2), new_points.index(p3)])
                new_faces.append([new_points.index(p3), new_points.index(p4), new_points.index(p5)])
                new_faces.append([new_points.index(p1), new_points.index(p3), new_points.index(p5)])
            # Recursion
            return self._subdivide(new_points, new_faces, times-1)

    def _plot_polygon(self, v, f, crop_z=-10, plane=True):
        """ Plot a polygon defined by vertices v and faces f."""

        # Create axes
        ax = a3.Axes3D(plt.figure())
        # Add polygon faces
        for i in range(len(f)):
            poly=[v[f[i][0]], v[f[i][1]], v[f[i][2]]]
            #add only if higher than crop distance
            if poly[0][2]>=crop_z and poly[1][2]>=crop_z and poly[2][2]>=crop_z:
                face = a3.art3d.Poly3DCollection([poly])
                face.set_alpha(1)
                face.set_color(colors.rgb2hex(sp.rand(3)))
                face.set_edgecolor('k')
                ax.add_collection3d(face)
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
        self.points, self.faces = self._subdivide(self.points, self.faces, n_times)
        self._plot_polygon(self.points, self.faces, 0)
        px, py, pz = self._xyz_list(self.points)
        self._plot_points(px, py, pz, 0)
        vs = [basis(i, self._point_normal(i)) for i in self.points]
        self._plot_bases(vs, 0.15)
        for i in range(len(vs)):
            a = quaternion()
            b = a.from_matrix(vs[i].matrix())

ob = EIT(1,1,3)
ob.tesselate(2)

#Objectives:
#   -From points to robtargets
# FIX quaternion to avoid division by zero
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
