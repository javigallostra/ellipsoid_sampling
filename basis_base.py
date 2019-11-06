from vector_base import vector
from math import radians, sin, cos

class basis:
    """ Basis class.

    Provides methods for creating, rotating and
    translating a basis as well as a method for
    computing the matrix representation of the
    basis.
    """

    def __init__(self, origin = 0, vz = 0):
        """ Set the origin and the basis x, y, z vectors.

        If an initial z vector is given, the basis is
        constructed orthonormal and with the z-component of
        the x vector equal to zero.

        Return nothing.
        """

        # Set origin
        if not origin:
            self.origin = vector(0,0,0)
        else:
            self.origin = vector(origin)
        # Set vectors
        if not vz:
            self.vx = vector(1,0,0)
            self.vy = vector(0,1,0)
            self.vz = vector(0,0,1)
        else:
            self.set_vz(vz)
        # Return
        return

    def set_vz(self, vz):
        """ Generate an orthonormal basis based on a vz vector.

        The assumption that the z compoment of the basis
        vx vector is zero is made in order to compute the
        orientation base.

        Return nothing.
        """

        # Normalize vz
        self.vz = vector(vz).normalized()
        # Check for vertical vz
        self.vx = vector(-vz[1], vz[0], 0)
        if self.vx.norm() == 0:
            self.vx = vector(1,0,0)
        else:
            self.vx.normalize()
        self.vy = self.vz.cross(self.vx).normalized()
        # Return
        return
        

    def rotate(self, alfa, axis='z'):
        """ Rotate the basis around one of its axis by alfa degrees.

        Rotate the basis using the axis-angle formula.

        Return nothing.
        """
        
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
        # Return
        return

    def translate(self, d, axis='z'):
        """ Translate the basis along one of its axis by d millimeters.

        Return nothing.
        """

        # Compute displacement vector
        if axis == 'x':
            displacement = self.vx * d
        elif axis == 'y':
            displacement = self.vy * d
        elif axis == 'z':
            displacement = self.vz * d
        else:
            raise ValueError("Unknown axis name: " + str(axis))
        # Move origin
        self.origin += displacement
        # Return
        return

    def matrix(self, transpose=True):
        """ Get the matrix representation of the basis orientation.

        The x, y, z basis vectors are stored as column vectors in
        the matrix. The matrix can also be constructed with x, y, z
        basis vectors as row vectors by enabling the 'transpose' option.

        Return the matrix as a list of vectors.
        """

        # Column vector matrix
        m = [vector(self.vx), vector(self.vy), vector(self.vz)]
        # Transpose if needed
        if transpose:
            mt = [m[j][i] for j in range(len(m[0])) for i in range(len(m))]
        else:
            mt = m
        # Return
        return mt
