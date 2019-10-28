import math
from random import uniform, choice
from point_base import point
from vector_base import vector
from basis_base import basis
from quaternion_base import quaternion
from ellipsoid_sampling import ellipsoid_sampling

class ETS(ellipsoid_sampling):
    """ Ellipsoidal Thomson Simulation.

    """

    def __init__(self, rx=1, ry=1, rz=1):
        super().__init__(rx, ry, rz)

    def _random_distribution(self, n_points):
        """
        """
        points = []
        for n in range(n_points):
            x = uniform(-self.rx, self.rx)
            y_max = math.sqrt((1 - x**2 / self.rx**2) * self.ry**2)
            y = uniform(-y_max, y_max)
            z_opt = math.sqrt((1 - x**2 / self.rx**2 - y**2 / self.ry**2) * self.rz**2)
            z = choice([-z_opt, z_opt])
            points.append(point(x, y, z))

        # Return
        return points

    def _point_ellipsoid_projection(self, p):
        # Project the point onto the ellipsoid
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
                y_den = math.sqrt((self.rz * dydz)**2 + self.ry**2)
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
                x_den = math.sqrt((self.rz * dxdz)**2 + self.rx**2)
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
                x_den = math.sqrt((self.ry * dxdy)**2 + self.rx**2)
                x = x_num / x_den
                if p.x < 0: x = -x
                #compute z
                y = x / dxdy
        # Quadrant points
        else:
            #find x
            dxdy = p.x / p.y
            dxdz = p.x / p.z
            x_num = self.rx * self.ry  *self.rz * abs(dxdy * dxdz)
            x_den = math.sqrt((dxdy * self.ry * dxdz * self.rz)**2 + (self.rx * dxdz * self.rz)**2 + (self.rx * dxdy * self.ry)**2)
            x = x_num / x_den
            if p.x < 0: x = -x
            #compute y, z
            y = x / dxdy
            z = x / dxdz
        # Return
        return point(x, y, z)

    def iteration(self, n_iter, pot_i):
                ###########################
        if n_iter == 0:
            print("EOI")
            return
        else:
            pot_i1 = 0
            k = 10
            vectors = []
            for p1 in self.points:
                v = vector(0, 0, 0)
                for p2 in self.points:
                    if p1 != p2:
                        diff = p1 - p2
                        v_diff = vector(diff.x, diff.y, diff.z)
                        v_normal = self._point_normal(p1)
                        v_plane = self._vector_plane_projection(v_diff, v_normal)
                        v_final = v_plane * (1 / v_diff.norm())
                        pot_i1 += 1 / v_diff.norm()
                        v = v + v_final
                vectors.append(v)
            print( pot_i, pot_i1)
            # threshood for potential change
            if abs(pot_i - pot_i1) < 0.01:
                print( pot_i, pot_i1)
                print("THSTOP: " + str(n_iter))
                return
            
            for i in range(len(self.points)):
                v = vectors[i]
                p = self.points[i] + point(v[0], v[1], v[2])
                p2 = self._point_ellipsoid_projection(p)
                dpsph = self._cart_to_spherical(p2) - self._cart_to_spherical(self.points[i])
                dvsph = vector(dpsph.lat, dpsph.long, 0).normalized()
                dvsph = dvsph * dvsph.norm() * k
                pf = self._cart_to_spherical(self.points[i]) + point(dvsph[0], dvsph[1], 0, 'spherical')
                pf = self._spherical_to_cart(pf)
                self.points[i] = pf

            return self.iteration(n_iter - 1, pot_i1)

        ############################

    def _vector_plane_projection(self, v, n):
        """ Project vector v onto plane with normal n. """
        
        v_proj = v + n * (v * n / n.norm())
        return v_proj

    def sample(self, n_points=50):
        self.points = self._random_distribution(n_points)
        self.basis = [basis(i, self._point_normal(i)) for i in self.points]

ob = ETS(1,1,1)
ob.sample(100)
ob.plot(True,False,False)
ob.iteration(50, 100000)
ob.plot(True,False,False)

#Objectives:
#   -From points to robtargets
# ¿Qué herramienta y qué workboject se usa? Pör ahora los puntos están centrados en la mesa.
