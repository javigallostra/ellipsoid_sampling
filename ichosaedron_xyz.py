import mpl_toolkits.mplot3d as a3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy as sp
import math
from random import randint

class EIT:
    """ Elipsoidal ichosaedron tesselation.

    Create an elipsoidal ichosaedron and perform
    triangular subdivisions (1 to 4) to tesellate its surface.
    """

    def __init__(self, rx=1, ry=1, rz=1, precision=3):
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.dec_prec = precision
        self.vertices, self.faces = self._ichosaedron()

    def _ichosaedron(self):
        """ Create a 'regular' elipsoidal ichosaedron.

        Use the radii defined in the constructor
        to shape the regular ichosaedron to the desired
        elipsoid.
        """
        
        # Create vertices in order (polar coord)
        v0 = [90,0]
        v11 = [-90,0]
        v1_5 = [[26.57, 2*i*36] for i in range(0,5)]
        v6_10 = [[-26.57, (2*i+1)*36] for i in range(0,5)]
        # Create vertex list in order (cartesian coord)
        v = [self._polar_to_cart(v0)]
        for i in range(0,5):
            v.append(self._polar_to_cart(v1_5[i]))
            v.append(self._polar_to_cart(v6_10[i]))
        v.append(self._polar_to_cart(v11))
        # Create faces with the ordered vertices
        f = []
        for i in range(0,4):
            f.append([0,2*i+1,2*(i+1)+1])
            f.append([11,2*(i+1),2*(i+2)])
        for i in range(0,8):
            f.append([i+1,i+2,i+3])
        f.append([0,1,9])
        f.append([11,2,10])
        f.append([9,10,1])
        f.append([10,1,2])
        # Return
        return v,f

    def _polar_to_cart(self, v):
        """ Polar to cartesian point mapping."""
        
        z = self.rz*math.sin(math.radians(v[0]))
        xy = math.cos(math.radians(v[0]))
        x = self.rx*xy*math.cos(math.radians(v[1]))
        y = self.ry*xy*math.sin(math.radians(v[1]))
        return [round(x, self.dec_prec), round(y, self.dec_prec), round(z, self.dec_prec)]

    def _cart_to_polar(self, v):
        """ Cartesian to polar point mapping."""

        x, y, z = v
        xy = math.sqrt(x**2 + y**2)
        merid = math.degrees(math.atan2(y, x))
        paral = math.degrees(math.atan2(z, xy))
        return [paral, merid]

    def _vertex_interpolation(self, v1, v2):
        """ Vertices interpolation function.

        Find the point equidistant to v1,v2
        lying on the elipsoid.
        """

        # Find the middle point of v1-v2
        vx = v1[0] + (v2[0]-v1[0])/2
        vy = v1[1] + (v2[1]-v1[1])/2
        vz = v1[2] + (v2[2]-v1[2])/2
        # Project the point onto the elipsoid
        # Axial points
        if [vx,vy,vz].count(0) == 2:
            x = self.rx*vx/abs(vx) if vx else 0
            y = self.ry*vy/abs(vy) if vy else 0
            z = self.rz*vz/abs(vz) if vz else 0
        # Planar points
        elif [vx,vy,vz].count(0) == 1:
            #plane x=0
            if vx == 0:
                x = 0
                #find y
                dydz = vy/vz
                y_num = self.ry*self.rz*abs(dydz)
                y_den = math.sqrt((self.rz*dydz)**2 + self.ry**2)
                y = y_num/y_den
                if vy < 0: y = -y
                #compute z
                z = y/dydz
            #plane y=0
            elif vy == 0:
                y = 0
                #find x
                dxdz = vx/vz
                x_num = self.rx*self.rz*abs(dxdz)
                x_den = math.sqrt((self.rz*dxdz)**2 + self.rx**2)
                x = x_num/x_den
                if vx < 0: x = -x
                #compute z
                z = x/dxdz
            #plane z=0
            elif vz == 0:
                z = 0
                #find x
                dxdy = vx/vy
                x_num = self.rx*self.ry*abs(dxdy)
                x_den = math.sqrt((self.ry*dxdy)**2 + self.rx**2)
                x = x_num/x_den
                if vx < 0: x = -x
                #compute z
                y = x/dxdy
        # Quadrant points
        elif vx!=0 and vy!=0 and vz!=0:
            #find x
            dxdy = vx/vy
            dxdz = vx/vz
            x_num = self.rx*self.ry*self.rz*abs(dxdy*dxdz)
            x_den = math.sqrt((dxdy*self.ry*dxdz*self.rz)**2 + (self.rx*dxdz*self.rz)**2 + (self.rx*dxdy*self.ry)**2)
            x = x_num/x_den
            if vx < 0: x = -x
            #compute y, z
            y = x/dxdy
            z = x/dxdz
        # Return
        return [round(x, self.dec_prec), round(y, self.dec_prec), round(z, self.dec_prec)]

    def _subdivide(self, v, f, times=1):
        """ Recursive triangular subdivision.

        Subdivide each surface triangle into 4
        a defined number of times to obtain a fine
        triangular tesselation.
        """

        # End case
        if times == 0:
            return (v, f)
        # Recursive case
        else:
            newf = []
            newv = []
            # Perform subdivision
            for i in range(len(f)):
                #copy existing vertices
                v0 = list(v[f[i][0]])
                v2 = list(v[f[i][1]])
                v4 = list(v[f[i][2]])
                #interpolate new vertices
                v1 = self._vertex_interpolation(v0,v2)
                v3 = self._vertex_interpolation(v2,v4)
                v5 = self._vertex_interpolation(v4,v0)
                #add vertices if not repeated
                if v0 not in newv: newv.append(v0)
                if v1 not in newv: newv.append(v1)
                if v2 not in newv: newv.append(v2)
                if v3 not in newv: newv.append(v3)
                if v4 not in newv: newv.append(v4)
                if v5 not in newv: newv.append(v5)
                #add faces
                newf.append([newv.index(v0),newv.index(v1),newv.index(v5)])
                newf.append([newv.index(v1),newv.index(v2),newv.index(v3)])
                newf.append([newv.index(v3),newv.index(v4),newv.index(v5)])
                newf.append([newv.index(v1),newv.index(v3),newv.index(v5)])
            # Recursion
            return self._subdivide(newv, newf, times-1)

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

    def _plot_points_bases(self, x, y, z, vx, vy, vz, crop_z=-10, plane=True, vec_scale=5):
        """ Plot a set of 3D points with coordinate bases."""

        # Create figure with axes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Filter points lower than crop distance
        xc = [x[i] for i in range(len(x)) if z[i] >=crop_z]
        yc = [y[i] for i in range(len(y)) if z[i] >=crop_z]
        zc = [z[i] for i in range(len(z)) if z[i] >=crop_z]
        vxc = [vx[i] for i in range(len(vx)) if z[i] >=crop_z]
        vyc = [vy[i] for i in range(len(vy)) if z[i] >=crop_z]
        vzc = [vz[i] for i in range(len(vz)) if z[i] >=crop_z]
        print("# points: " + str(len(xc)))
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
        
        x = [i[0] for i in vs]
        y = [i[1] for i in vs]
        z = [i[2] for i in vs]
        return x,y,z

    def _xyz_to_normal(self, p):
        """ Compute the normal of an ellipsoid surface point.

        The normal is unitary and facing inwards the ellipsoid.
        """
        
        px, py, pz = p
        v = [2*px/(self.rx**2), 2*py/(self.ry**2), 2*pz/(self.rz**2)]
        modulo = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        normal = self._times(v, 1/modulo)
        return normal

    def _oriented_base(self, p, z_vector, flip_threshold=0.6):
        """ Compute the oriented base of the z-vector of the base.

        The assumption that the z component of the basis x-vector is zero
        is made in order to compute the orientation base.
        Flip_threshold is used to change the camera orientation based
        on the height of the measurement points.
        """
        
        x_modulo = math.sqrt(z_vector[0]**2 + z_vector[1]**2)
        # Check for vertical vectors
        if x_modulo != 0:
                x_vector = [-z_vector[1]/x_modulo, z_vector[0]/x_modulo, 0]
        else:
            x_vector = [1, 0, 0]
        y_vector = self._cross(z_vector, x_vector)

        return self._base_rotate_around_z(x_vector, y_vector, z_vector, randint(-4, 4)*90/4)
##        if p[2] >= flip_threshold:
##            return self._base_rotate_around_z(x_vector, y_vector, z_vector, 180)
##        else:
##            return [x_vector, y_vector, z_vector]

    def _base_rotate_around_z(self, vx, vy, vz, alfa):
        """ Rotate a vector basis around its z vector."""

        rad = math.radians(alfa)
        # Rotate x vector using angle-axis formula
        new1 = self._times(vz, self._dot(vz, vx))
        new2 = self._times(self._cross(self._cross(vz, vx), vz), math.cos(rad))
        new3 = self._times(self._cross(vz, vx), math.sin(rad))
        new_vx = [new1[0]+new2[0]+new3[0], new1[1]+new2[1]+new3[1], new1[2]+new2[2]+new3[2]]
        # Compute the new y vector using the cross product
        new_vy = self._cross(vz, new_vx)
        
        return [new_vx, new_vy, vz]

    def _cross(self, u, v):
        """ Compute the cross product u^v."""

        return [u[1]*v[2]-v[1]*u[2], -u[0]*v[2]+v[0]*u[2], u[0]*v[1]-v[0]*u[1]]

    def _dot(self, u, v):
        """ Compute the dot product u·v."""

        return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

    def _times(self, u, k):
        """ Compute the scalar multiplication k*u."""

        return [u[0]*k, u[1]*k, u[2]*k]

    def _rotation_to_quaternion(self, vx, vy, vz):
        """ Compute the quaternion (x,y,z,w) of a rotation matrix.

        If the rotation matrix is orthonormal, the resulting quaternion
        will be unitary.
        """

        # Using max to avoid numerical errors of close to zero negative numbers
        w = 0.5*math.sqrt(max(0,1+vx[0]+vy[1]+vz[2]))
        x = (vy[2]-vz[1])/(4*w)
        y = (vz[0]-vx[2])/(4*w)
        z = (vx[1]-vy[0])/(4*w)

        return [[x,y,z,w], math.sqrt(x**2+y**2+z**2+w**2)]
        

    def tesselate(self, n_times=1):
        self.vertices, self.faces = self._subdivide(self.vertices, self.faces, n_times)
        #self._plot_polygon(self.vertices, self.faces, 0)
        px, py, pz = self._xyz_list(self.vertices)
        #self._plot_points(px, py, pz, 0)
        vs = [self._oriented_base(i, self._xyz_to_normal(i)) for i in self.vertices]
        self._plot_points_bases(px, py, pz, [vs[i][0] for i in range(len(vs))], [vs[i][1] for i in range(len(vs))], [vs[i][2] for i in range(len(vs))], 0.15)
        #for i in range(len(vs)):
        #    print(self._rotation_to_quaternion(vs[i][0], vs[i][1], vs[i][2]))

ob = EIT(1,1,1)
ob.tesselate(2)

#Objectives:
#   -From points to robtargets
# FIX quaternion to avoid division by zero
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
