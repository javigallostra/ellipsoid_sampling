from vector_base import vector
from math import radians, sin, cos

class basis:
    """ Basis class."""

    def __init__(self):
        self.origin = vector(0,0,0)
        self.vx = vector(1,0,0)
        self.vy = vector(0,1,0)
        self.vz = vector(0,0,1)

    def set_vz(self, vz):
        """ Generate a basis based on a vz vector.

        The assumption that the z compoment of the basis
        vx vector is zero is made in order to compute the
        orientation base.
        """
        
        self.vz = vector(vz).normalized()
        # Check for vertical vz
        if self.vz[2] != 1:
            self.vx = vector(-vz[1], vz[0], 0).normalized()
        else:
            self.vx = vector(1,0,0)
        self.vy = vz.cross(vx).normalized()
        

    def rotate_z(self, alfa):
        """ Rotate the basis around its z-vector by alfa degrees."""
        
        rad = radians(alfa)
        # Rotation using angle-axis formula
        self.vx = self.vz*(self.vx*self.vz) + self.vz.cross(self.vx).cross(self.vz)*cos(rad) + self.vz.cross(self.vx)*sin(rad)
        self.vy = self.vz.cross(self.vx)
