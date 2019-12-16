from math import sin, cos, sqrt, degrees, radians, atan2
from point_base import point
from vector_base import vector

class ellipsoid_base:
    """ Basic ellipsoid class.

    Provides methods to compute the normal vector
    of a surface points, to switch a point between cartesian
    and spherical coordinates, and to project any point in
    space onto the ellipsoidal surface.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters.

        Create an empty list of points.

        Return nothing.
        """
        
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.points = []
        return

    def _point_normal(self, p):
        """ Compute the normal of an ellipsoid surface point.

        The normal is unitary and facing inwards the ellipsoid.

        Return the resulting normal vector.
        """
        
        # Compute vector and normalize
        vx = 2 * p.x / (self.rx**2)
        vy = 2 * p.y / (self.ry**2)
        vz = 2 * p.z / (self.rz**2)
        normal = vector(vx, vy, vz).normalized()
        # return
        return normal

    def _spherical_to_cart(self, p):
        """ Spherical to cartesian point mapping.

        Return the resulting point.
        """

        # Compute coordinates
        z = self.rz * sin(radians(p.lat))
        xy = cos(radians(p.lat))
        x = self.rx * xy * cos(radians(p.long))
        y = self.ry * xy * sin(radians(p.long))
        # Return
        return point(x, y, z)

    def _cart_to_spherical(self, p):
        """ Cartesian to spherical point mapping.

        Return the resulting point.
        """

        # Compute coordinates
        # @todo: check that lat is lat and long is long...
        xy = sqrt(p.x**2 + p.y**2)
        merid = degrees(atan2(p.y, p.x))
        paral = degrees(atan2(p.z, xy))
        # Return
        return point(paral, merid, 0, "spherical")

    def _point_ellipsoid_projection(self, p):
        """ Project a point onto the ellipsoid surface.

        Return the resulting point.
        """
        
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
                y_den = sqrt((self.rz * dydz)**2 + self.ry**2)
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
                x_den = sqrt((self.rz * dxdz)**2 + self.rx**2)
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
                x_den = sqrt((self.ry * dxdy)**2 + self.rx**2)
                x = x_num / x_den
                if p.x < 0: x = -x
                #compute z
                y = x / dxdy
        # Quadrant points
        else:
            #find x
            dxdy = p.x / p.y
            dxdz = p.x / p.z
            x_num = self.rx * self.ry  * self.rz * abs(dxdy * dxdz)
            x_den = sqrt((dxdy * self.ry * dxdz * self.rz)**2 + (self.rx * dxdz * self.rz)**2 + (self.rx * dxdy * self.ry)**2)
            x = x_num / x_den
            if p.x < 0: x = -x
            #compute y, z
            y = x / dxdy
            z = x / dxdz
        # Return
        return point(x, y, z)
