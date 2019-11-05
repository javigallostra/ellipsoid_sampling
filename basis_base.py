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
        

    def rotate(self, alfa, axis='z'):
        """ Rotate the basis around one of its axis by alfa degrees."""
        
        rad = radians(alfa)
        # Prepare rotation
        if axis == 'x':
            a = self.vx
            b = self.vy
        elif axis == 'y':
            a = self.vy
            b = self.vz
        elif axis == 'z':
            a = self.vz
            b = self.vx
        else:
            raise ValueError("Unknown axis: " + str(axis))
        # Rotation using angle-axis formula
        new_b = a * (b * a) + a.cross(b).cross(a) * cos(rad) + a.cross(b) * sin(rad)
        new_c = a.cross(new_b)
        # Save new values
        if axis == 'x':
            self.vy = new_b
            self.vz = new_c
        elif axis == 'y':
            self.vx = new_c
            self.vz = new_b
        elif axis == 'z':
            self.vx = new_b
            self.vy = new_c
        else:
            raise ValueError("Unknown axis: " + str(axis))

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
