import math
from point_base import point
from vector_base import vector
from basis_base import basis
from quaternion_base import quaternion
from ellipsoid_sampling import ellipsoid_sampling

class EFS(ellipsoid_sampling):
    """ Ellipsoidal fibonacci spiral.

    Create an ellipsoidal distribution of N points
    using the Fibonacci Sprial.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        super().__init__(rx,ry,rz)

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

    def sample(self, n_points):
        self.points = self._spiral(n_points)
        self.basis = [basis(p, self._point_normal(p)) for p in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh

ob = EFS(1,2,1)
ob.sample(100)
ob.plot(False,True,False,0)

#Objectives:
#   -From points to robtargets
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
