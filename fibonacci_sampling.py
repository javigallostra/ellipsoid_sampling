import math
from point_base import point
from vector_base import vector
from basis_base import basis
from quaternion_base import quaternion
from ellipsoid_sampling import ellipsoid_sampling

class EFS(ellipsoid_sampling):
    """ Ellipsoidal Fibonacci Spiral.

    Class to sample an ellipsoid with an arbitrary
    number of points. The points are distributed using
    the Fibonacci Spiral algorithm on a unit sphere
    and then backprojected to the ellipsoidal surface.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters.

        Return nothing.
        """
        
        super().__init__(rx,ry,rz,"fibonacci")
        return

    def _spiral(self, n_points):
        """ Distribute n_points on a Fibonacci Spiral.

        The spiral is created on a unit sphere and then
        the points are backprojected to the ellipsoidal
        surface.

        Return the list of points.
        """

        # Create spiral on unit sphere and project onto ellipsoid
        dlong = math.pi*(3 - math.sqrt(5))
        dz = 2.0 / n_points
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

    def sample(self, n_points):
        """ Sample the ellipsoid using the Fibonacci Spiral algorithm.

        Get the points and compute their orthogonal basis.

        Return nothing.
        """
        
        self.points = self._spiral(n_points)
        self.basis = [basis(p, self._point_normal(p)) for p in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh
        return
