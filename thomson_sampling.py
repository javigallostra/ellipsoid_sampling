import math
from random import uniform, choice
from point_base import point
from vector_base import vector
from basis_base import basis
from quaternion_base import quaternion
from fibonacci_sampling import EFS

class ETS(EFS):
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

    def _iterate(self, n_iter, pot_i, stop_threshold):
        """
        """

        # Ending contition 1: max cycles reached
        if n_iter == 0:
            print("Simulation reached maximum number of steps")
            return
        # Max cycles not reached: compute potential
        else:
            pot_i1 = 0
            k = 1 # @todo: not fixed?
            vectors = []
            # Sum potential of each point w.r.t all the others
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
            print(pot_i1)
            # Ending condition 2: potential change below threshold
            if abs(pot_i - pot_i1) < stop_threshold:
                print("Simulation converged with " + str(n_iter) + " steps left")
                return
            # Continue simulation: move points and repeat
            for i in range(len(self.points)):
                # Compute point + vector
                v = vectors[i]
                p = self.points[i] + point(v[0], v[1], v[2])
                # Project onto ellipsoid
                p2 = self._point_ellipsoid_projection(p)
                # Find unitary displacement vector in spherical coordinates
                dpsph = self._cart_to_spherical(p2) - self._cart_to_spherical(self.points[i])
                dvsph = vector(dpsph.lat, dpsph.long, 0).normalized()
                # Scale vector
                dvsph = dvsph * dvsph.norm() * k
                # Move point in spherical coordinates
                pf = self._cart_to_spherical(self.points[i]) + point(dvsph[0], dvsph[1], 0, 'spherical')
                # Save changed point in cartesian coordinates
                pf = self._spherical_to_cart(pf)
                self.points[i] = pf
            # Recursive iteration
            return self._iterate(n_iter - 1, pot_i1, stop_threshold)

    def _vector_plane_projection(self, v, n):
        """ Project vector v onto plane with normal n. """
        
        v_proj = v + n * (v * n / n.norm())
        return v_proj

    def sample(self, n_points=50, n_iterations=100, stop_threshold=0.01):
        # 1 - unit sphere
        ellipsoid_radii = [self.rx, self.ry, self.rz]
        self.rx = self.ry = self.rz = 1
        # 2 - iterate
        self.points = self._spiral(n_points)
        self._iterate(n_iterations, n_points*100, stop_threshold)
        # 3 - backproject
        self.rx, self.ry, self.rz = ellipsoid_radii
        for i in range(len(self.points)):
            self.points[i] = self._point_ellipsoid_projection(self.points[i])
        self.basis = [basis(i, self._point_normal(i)) for i in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh

#Objectives:
#   -From points to robtargets
# ¿Qué herramienta y qué workboject se usa? Por ahora los puntos están centrados en la mesa.

# @todo: fix initial orientation to avoid randomness between runs with same number of points
# @todo: implement a stopping criterion to detect stable oscillation
