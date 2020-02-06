import math
from base_classes.point_base import point
from base_classes.vector_base import vector
from base_classes.basis_base import basis
from base_classes.quaternion_base import quaternion
from base_classes.ellipsoid_sampling_base import ellipsoid_sampling_base

class EFS(ellipsoid_sampling_base):
    """ Ellipsoidal Fibonacci Spiral.

    Class to sample an ellipsoid with an arbitrary
    number of points. The points are distributed using
    the Fibonacci Spiral algorithm directly on the
    ellipsoidal surface. Alternatively they can be distributed
    on a unit sphere and then backprojected to the ellipsoidal surface.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters.

        Return nothing.
        """
        
        super().__init__(rx,ry,rz,"fibonacci")
        return

    def _sphere_lowest_z(self, ellipsoid_z):
        """ Compute the lowest z value of a sphere point whose projection
        onto the ellipsoid has the given ellipsoid_z value.

        It works by taking the largest of x/y axis, finding the vector on
        the corresponding x=0/y=0 plane whose tip lies on the ellipsoid surface
        at z = ellipsoid_z. Its unitary vector lies on the unitary sphere surface
        and gives the lowest z value.

        Return the lowest z value.
        """

        longest_axis_sq = max(self.rx, self.ry)**2
        rz_sq = self.rz**2
        norm = math.sqrt(longest_axis_sq + ((rz_sq - longest_axis_sq) / rz_sq) * ellipsoid_z**2)
        lowest_z = ellipsoid_z / norm

        return lowest_z

    def _spiral_sphere(self, n_points, min_z=None):
        """ Distribute n_points on a Fibonacci Spiral.

        The spiral is created on a unit sphere and then
        the points are backprojected to the ellipsoidal
        surface.

        An optional min_z parameter is available to decide
        the z coordinate of the lowest point.

        Return the list of points.
        """

        # If no limit is given, take the full sphere
        if min_z == None:
            min_z = -1.0
        else:
            min_z = self._sphere_lowest_z(min_z)
        # Create spiral on unit sphere and project onto ellipsoid
        dlong = math.pi*(3 - math.sqrt(5))
        dz = (1.0 - min_z) / n_points
        long = 0
        z = 1 - dz/2
        points = []
        for i in range(n_points):
            r = math.sqrt(1-z*z)
            p = point(math.cos(long)*r, math.sin(long)*r, z)
            p = self._point_ellipsoid_projection(p)
            points.append(p)
            z -=  dz
            long += dlong
        # Return
        return points

    def _spiral(self, n_points, min_z=None):
        """ Distribute n_points on a Fibonacci Spiral.

        The spiral is created on the ellipsoidal
        surface.

        An optional min_z parameter is available to decide
        the z coordinate of the lowest point.

        Return the list of points.
        """
        
        # If no limit is given, take the full ellipsoid
        if min_z == None:
            min_z = -self.rz
        # Create spiral on ellipsoid
        dlong = math.pi*(3 - math.sqrt(5))
        dz = (self.rz - min_z) / n_points
        long = 0
        z = self.rz - dz/2
        points = []
        for i in range(n_points):
            # Solve ellipsoid equation
            num = self.rx**2 * self.ry**2 * (1 - (z**2) / self.rz**2)
            den = self.ry**2 + self.rx**2 * math.tan(long)**2
            x = math.sqrt(num / den)
            # Solve x's sign
            angle = math.degrees(long) % 360
            if angle > 90 and angle < 270:
                x = -x
            # Compute y
            y = x * math.tan(long)
            # Generate point
            p = point(x, y, z)
            points.append(p)
            # Next point
            z -=  dz
            long += dlong
        # Return
        return points

    def sample(self, n_points, min_z=None):
        """ Sample the ellipsoid using the Fibonacci Spiral algorithm.

        Get the points and compute their orthogonal basis.
        An optional minimum z value for the lowest desired point can be given.
        Otherwise the Fibonacci Spiral algorithm will sample n_points along
        the whole ellipsoid height.

        Return nothing.
        """
        
        self.points = self._spiral(n_points, min_z)
        self.basis = [basis(p, self._point_normal(p)) for p in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh
        return
