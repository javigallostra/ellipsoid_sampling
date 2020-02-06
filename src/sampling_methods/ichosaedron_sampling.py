import math
from base_classes.point_base import point
from base_classes.vector_base import vector
from base_classes.basis_base import basis
from base_classes.quaternion_base import quaternion
from base_classes.ellipsoid_sampling_base import ellipsoid_sampling_base

class EIS(ellipsoid_sampling_base):
    """ Ellipsoidal Ichosaedron Subdivision.

    Class to sample an ellipsoid using succesive Ichosaedron
    Subdivisions. The ichosaedron is created on the ellipsoid
    using spherical coordinates. Its triangles are then subdivided
    into four smaller triangles an arbitrary number of times. The
    vertices of said triangles are the sampled points.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters and create the regular ichosaedron.

        After getting the ichosaedron points, their orthonormal basis as
        well as their gathering into faces are computed.
        """
        
        super().__init__(rx, ry, rz, "ichosaedron")
        self.points, self.faces = self._ichosaedron()
        self.basis = [basis(i, self._point_normal(i)) for i in self.points]

    def _ichosaedron(self):
        """ Create a regular ichosaedron backprojected to the ellipsoid.

        Compute the ichosaedron points with spherical coordinates and
        backproject them to the ellipsoida surface in the cartesian space.
        After that, group the ichosaedron points into faces.

        Return a tuple with the list of points and the list of faces.
        """

        # @ todo: check if we can use bakcprojection instead of cartesian to spherical
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
        return (points, faces)

    def _point_interpolation(self, p1, p2):
        """ Find the middle ellipsoidal point between two points.

        Compute the middle point of vector p1-p2 and backproject
        it to the ellipsoidal surface.

        Return the resulting point.
        """

        # Find the middle point of p1-p2
        p = p1 + (p2 - p1)/2
        # Project the point onto the ellipsoid
        p = self._point_ellipsoid_projection(p)
        # Return
        return p

    def _subdivide(self, points, faces, times=1):
        """ Recursively subdivide the triangles.

        Subdivide each surface triangle into 4
        a defined number of times to obtain a fine
        triangular tesselation.

        Return a tuple with the list of points and the list of faces.
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

    def sample(self, n_times=1):
        """ Sample the ellipsoid by succesive Ichosaedron Subdivisions.

        Get the points and compute their orthogonal basis.

        Return nothing.
        """
        self.points, self.faces = self._subdivide(self.points, self.faces, n_times)
        self.basis = [basis(i, self._point_normal(i)) for i in self.points]
        return
