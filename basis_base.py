from vector_base import vector
from math import radians, sin, cos

class basis:
    """ Basis class."""

    def __init__(self, origin = 0, vz = 0):
        if not origin:
            self.origin = vector(0,0,0)
        else:
            self.origin = vector(origin)
        if not vz:
            self.vx = vector(1,0,0)
            self.vy = vector(0,1,0)
            self.vz = vector(0,0,1)
        else:
            self.set_vz(vz)

    def set_vz(self, vz):
        """ Generate a basis based on a vz vector.

        The assumption that the z compoment of the basis
        vx vector is zero is made in order to compute the
        orientation base.
        """
        
        self.vz = vector(vz).normalized()
        # Check for vertical vz
        self.vx = vector(-vz[1], vz[0], 0)
        if self.vx.norm() == 0:
            self.vx = vector(1,0,0)
        else:
            self.vx.normalize()
        self.vy = self.vz.cross(self.vx).normalized()
        

    def rotate_z(self, alfa):
        """ Rotate the basis around its z-vector by alfa degrees."""
        
        rad = radians(alfa)
        # Rotation using angle-axis formula
        self.vx = self.vz*(self.vx*self.vz) + self.vz.cross(self.vx).cross(self.vz)*cos(rad) + self.vz.cross(self.vx)*sin(rad)
        self.vy = self.vz.cross(self.vx)

    def rotate_y(self, alfa):
        """ Rotate the basis around its y-vector by alfa degrees."""
        
        rad = radians(alfa)
        # Rotation using angle-axis formula
        self.vz = self.vy*(self.vz*self.vy) + self.vy.cross(self.vz).cross(self.vy)*cos(rad) + self.vy.cross(self.vz)*sin(rad)
        self.vx = self.vy.cross(self.vz)

    def translate(self, d, axis='z'):
        """ Translate the basis along one of its axis by d millimeters."""

        if axis == 'x':
            displacement = self.vx * d
        elif axis == 'y':
            displacement = self.vy * d
        elif axis == 'z':
            displacement = self.vz * d
        else:
            raise ValueError("Unknown axis name: " + str(axis))
        self.origin += displacement

    def matrix(self, transpose=True):
        """ Get the matrix representation of the basis orientation.

        The matrix is returned as a list of vectors.
        """

        m = [vector(self.vx), vector(self.vy), vector(self.vz)]
        if transpose:
            mt = [[0,0,0],[0,0,0],[0,0,0]]
            mt[0][0] = m[0][0]
            mt[1][1] = m[1][1]
            mt[2][2] = m[2][2]
            mt[0][1] = m[1][0]
            mt[0][2] = m[2][0]
            mt[1][0] = m[0][1]
            mt[1][2] = m[2][1]
            mt[2][0] = m[0][2]
            mt[2][1] = m[1][2]
        else:
            mt = m
        return mt
