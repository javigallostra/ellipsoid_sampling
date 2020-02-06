import math
from random import uniform, choice
from base_classes.point_base import point
from base_classes.vector_base import vector
from base_classes.basis_base import basis
from base_classes.quaternion_base import quaternion
from sampling_methods.fibonacci_sampling import EFS

class ETS(EFS):
    """ Ellipsoidal Thomson Simulation.

    Class used to sample an ellipsoid with an arbitrary
    number of points. The points are distributed by running
    a Thomson Simulation on a unit sphere and then backprojecting
    the final points onto the ellipsoidal surface.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters.

        Return nothing.
        """
        
        super().__init__(rx, ry, rz, "thomson")
        self.stop_threshold = 0.00001
        self.max_iterations = 100
        self.max_oscillations = 6
        self.max_growing = 3
        return

    def _random_distribution(self, n_points):
        """ Randomly distribute n_points over the ellipsoid surface.

        Return the list of points.
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

    def _iterate(self, n_iter, pot_i, min_pot, oscillations, growing):
        """ Recursively iterate over Thomson simulation steps.

        At each step, compute the points' potential, check all the stopping
        criteria and if none is met it then proceed to update the position
        of the points.

        Return nothing.
        """
        
        ##
        # Pre-step 1: compute current potential
        ##
        pot_i1, force_vectors = self._compute_potential()
        print("[SAMPLER-THOMSON] - Iteration: %d Potential: %f" %(n_iter, pot_i1))
        ##
        # Pre-step 2: Update stopping criteria variables
        ##
        pot_incr = pot_i - pot_i1
        # oscillations
        if growing and (pot_incr > 0):
            oscillations += 1
        elif (not growing) and (pot_incr < 0):
            oscillations += 1
        else:
            oscillations = 0
        # divergence
        if pot_incr < 0:
            growing += 1
        else:
            growing = 0
        ##
        # Step 1: Check stopping criteria
        ##
        # Ending contition 1: max cycles reached
        if n_iter == self.max_iterations:
            print("[SAMPLER-THOMSON] - Simulation reached maximum number of steps (%d)." %(n_iter))
            return
        # Ending condition 2: potential increment below threshold
        elif abs(pot_incr) < self.stop_threshold:
            print("[SAMPLER-THOMSON] - Simulation converged at iteration %d." %(n_iter))
            return
        # Ending condition 3: oscillatory state
        elif oscillations == self.max_oscillations:
            print("[SAMPLER-THOMSON] - Simulation oscillates at iteration %d." %(n_iter))
            return
        # Ending condition 4: divergence
        elif growing == self.max_growing:
            print("[SAMPLER-THOMSON] - Simulation diverges. Use the Fibonacci approximation.")
            return
        ##
        # Step 2: Continue iterations
        ##
        # No stopping: move points and repeat
        else:     
            # Continue simulation: move points and repeat
            for i in range(len(self.points)):
                # Compute point + vector
                force_vector = force_vectors[i]
                p_forced = self.points[i] + point(force_vector[0], force_vector[1], force_vector[2])
                # Project onto ellipsoid
                p_forced_proj = self._point_ellipsoid_projection(p_forced)
                # Find unitary displacement vector in spherical coordinates
                dp_sph = self._cart_to_spherical(p_forced_proj) - self._cart_to_spherical(self.points[i])
                dv_sph = vector(dp_sph.lat, dp_sph.long, 0).normalized()
                # Move point in spherical coordinates
                p_final = self._cart_to_spherical(self.points[i]) + point(dv_sph[0], dv_sph[1], 0, 'spherical')
                # Save changed point in cartesian coordinates
                p_final = self._spherical_to_cart(p_final)
                self.points[i] = p_final
            # Update minimum potential
            if pot_i1 < min_pot:
                min_pot = pot_i1
            # Update iterations
            n_iter += 1
        # Recursive iteration
        return self._iterate(n_iter, pot_i1, min_pot, oscillations, growing)

    def _compute_potential(self):
        """ Compute the potential of a set of points.

        The potential is the sum of the inverse distances
        between all possible point pairs. The list of force vectors
        for each point is also computed. A force vector is the planar
        projection of the a distance vector multiplied by the inverse
        of the distance from one point to another. All force vectors
        are computed for one point with all its possible point pairs
        and the resulting vector sum is projected to the
        plane tangent at the surface at the considered point.

        Return a tuple with the total poential and the list of force vectors.
        """

        potential = 0
        force_vectors = []
        n_points = len(self.points)
        # Matrix to store relative vectors to speed up the process
        # @ todo: test if it really helps
        potential_matrix = [[-1 for j in range(n_points)] for i in range(n_points)]
        # Sum potential of each point w.r.t all the others
        for i in range(n_points):
            p1 = self.points[i]
            force_vector = vector(0, 0, 0)
            v_normal = self._point_normal(p1)
            for j in range(n_points):
                # Case 1: same point
                if i == j:
                    continue
                # Case 2: point not paired yet
                elif potential_matrix[i][j] == -1:
                    p2 = self.points[j]
                    diff = p1 - p2
                    v_diff = vector(diff.x, diff.y, diff.z)
                    root_potential = 1 / v_diff.norm()
                    potential_matrix[i][j] = (v_diff, root_potential)
                    potential_matrix[j][i] = (v_diff * -1, root_potential)
                # Case 3: point already paired
                else:
                    v_diff, root_potential = potential_matrix[i][j]
                # Project v_diff to plane tangent to ellipsoid at p1
                v_plane = self._vector_plane_projection(v_diff, v_normal)
                # Scale projected vector with 1/|r|
                v_final = v_plane * root_potential
                # Add force vector
                force_vector += v_final
                # Add potential
                if j > i:
                    potential += root_potential
            # Looped through all j for point i: append resultant force vector
            force_vectors.append(force_vector)
        # Return
        return (potential, force_vectors)

    def _vector_plane_projection(self, v, n):
        """ Project vector v onto plane with normal n.

        Return the projected vector.
        """
        
        v_proj = v + n * (v * n / n.norm())
        return v_proj

    def sample(self, n_points=50):
        """ Sample n_points on the ellipsoid using the Thomson simulation.

        Begin by distributing n_points over the unit sphere using the
        Fibonacci Spiral algorihtm. Run a Thomson simulation on the given
        points and backproject them to the ellipsoidal surface.

        After getting the points, compute their orthonormal basis.

        Return nothing.
        """
        
        # Unit sphere Fibonacci distribution
        ellipsoid_radii = [self.rx, self.ry, self.rz]
        self.rx = self.ry = self.rz = 1
        # Thomson Simulation iterations
        self.points = self._spiral(n_points)
        self._iterate(0, n_points*100, 0, 0, 0)
        # Backprojection
        self.rx, self.ry, self.rz = ellipsoid_radii
        self.points = [self._point_ellipsoid_projection(p) for p in self.points]
        self.basis = [basis(p, self._point_normal(p)) for p in self.points]
        self.faces = [] # @todo: maybe implement a method to make a mesh
        return
