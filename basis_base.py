from vector_base import vector
from math import radians, sin, cos

class basis:
    """ Basis class."""

    def __init__(self):
        self.origin = vector(0,0,0)
        self.vx = vector(1,0,0)
        self.vy = vector(0,1,0)
        self.vz = vector(0,0,1)

    def rotate_z(self, alfa):
        rad = radians(alfa)
        self.vx = self.vz*(self.vx*self.vz) + self.vz.cross(self.vx).cross(self.vz)*cos(rad) + self.vz.cross(self.vx)*sin(rad)
        self.vy = self.vz.cross(self.vx)
