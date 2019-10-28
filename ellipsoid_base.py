from math import sin, cos, sqrt, degrees, radians, atan2
from point_base import point
from vector_base import vector

class ellipsoid_base:

    def __init__(self, rx=1, ry=1, rz=1):
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.points = []

    def _point_normal(self, p):
        """ Compute the normal of an ellipsoid surface point.

        The normal is unitary and facing inwards the ellipsoid.
        """
        
        vx = 2 * p.x / (self.rx**2)
        vy = 2 * p.y / (self.ry**2)
        vz = 2 * p.z / (self.rz**2)
        normal = vector(vx, vy, vz).normalized()
        
        return normal

    def _spherical_to_cart(self, p):
        """ Spherical to cartesian point mapping."""
        
        z = self.rz * sin(radians(p.lat))
        xy = cos(radians(p.lat))
        x = self.rx * xy * cos(radians(p.long))
        y = self.ry * xy * sin(radians(p.long))
        
        return point(x, y, z)

    def _cart_to_spherical(self, p):
        """ Cartesian to spherical point mapping."""

        xy = sqrt(p.x**2 + p.y**2)
        merid = degrees(atan2(p.y, p.x))
        paral = degrees(atan2(p.z, xy))
        
        return point(paral, merid, 0, "spherical")
