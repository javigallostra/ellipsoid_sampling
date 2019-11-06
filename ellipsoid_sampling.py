import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import rand
from ellipsoid_base import ellipsoid_base

class ellipsoid_sampling(ellipsoid_base):
    """ Ellipsoid sampling class.

    Enhances the base class with some attributes
    and useful methods.
    """

    def __init__(self, rx=1, ry=1, rz=1):
        """ Set the initial parameters.

        Create empty lists of basis and faces.

        Return nothing.
        """
        
        super().__init__(rx,ry,rz)
        self.basis = []
        self.faces = []
        return

    def plot(self, points=True, basis=False, faces=False, crop_z=None):
        """ Plot the ellipsoid.

        Plot the points, basis and/or faces according
        to the selected options. If crop_z is given, draw
        a plane at the same height and do not plot the
        elements below the plane.

        Return nothing.
        """

        # Create axes
        ax = a3.Axes3D(plt.figure())
        # Plot points if desired
        if points:
            # Crop
            if crop_z != None:
                points_crop = [p for p in self.points if p[2] >=crop_z]
            else:
                points_crop = self.points
            # Prepare
            px, py, pz = self._xyz_list(points_crop)
            # Plot
            ax.scatter(px, py, pz, c='r', marker='o')
        # Plot basis if desired
        if basis:
            # Crop
            if crop_z != None:
                basis_crop = [b for b in self.basis if b.origin[2] >= crop_z]
            else:
                basis_crop = self.basis
            # Prepare
            px = [b.origin[0] for b in basis_crop]
            py = [b.origin[1] for b in basis_crop]
            pz = [b.origin[2] for b in basis_crop]
            vx = [b.vx for b in basis_crop]
            vy = [b.vy for b in basis_crop]
            vz = [b.vz for b in basis_crop]
            # Plot
            scale = 5/max([self.rx, self.ry, self.rz])
            ax.quiver(px, py, pz, [v[0]/scale for v in vx], [v[1]/scale for v in vx], [v[2]/scale for v in vx], color="red", arrow_length_ratio=0.5)
            ax.quiver(px, py, pz, [v[0]/scale for v in vy], [v[1]/scale for v in vy], [v[2]/scale for v in vy], color="green", arrow_length_ratio=0.5)
            ax.quiver(px, py, pz, [v[0]/scale for v in vz], [v[1]/scale for v in vz], [v[2]/scale for v in vz], color="blue", arrow_length_ratio=0.5)
        # Plot faces if desired
        if faces:
            # Plot each face
            for f in self.faces:
                poly=[self.points[f[0]], self.points[f[1]], self.points[f[2]]]
                # Crop
                if crop_z != None:
                    if not(poly[0][2]>=crop_z and poly[1][2]>=crop_z and poly[2][2]>=crop_z):
                        continue
                # Plot
                face = a3.art3d.Poly3DCollection([poly])
                face.set_alpha(1)
                face.set_color(colors.rgb2hex(rand(3)))
                face.set_edgecolor('k')
                ax.add_collection3d(face)
        # Add visual plane
        if crop_z != None:
            poly = [[self.rx,self.ry,crop_z],[-self.rx,self.ry,crop_z],[-self.rx,-self.ry,crop_z],[self.rx,-self.ry,crop_z]]
            face = a3.art3d.Poly3DCollection([poly])
            face.set_alpha(0.5)
            face.set_color(colors.rgb2hex(rand(3)))
            face.set_edgecolor('k')
            ax.add_collection3d(face)
        # Add labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # Scale
        if crop_z != None:
            d = max([self.rx, self.ry, (self.rz - crop_z)/2])
            z_c = crop_z + (self.rz - crop_z)/2
            ax.auto_scale_xyz([-d, d], [-d, d], [z_c - d, z_c + d])
        else:
            d = max([self.rx, self.ry, self.rz])
            ax.auto_scale_xyz([-d, d], [-d, d], [-d, d])
        # Call matplotlib plot()
        plt.show()
        return

    def _xyz_list(self, vs):
        """ Transform a list of 3D points into three X, Y, Z lists.

        Return a tuple with the X, Y, Z, lists.
        """
        
        x = [vi[0] for vi in vs]
        y = [vi[1] for vi in vs]
        z = [vi[2] for vi in vs]
        return (x,y,z)

    def crop_z(self, crop_z):
        """ Remove the elements below the given z.

        Remove points, faces and basis.

        Return nothing.
        """

        new_faces = []
        for f in self.faces:
            poly=[self.points[f[0]], self.points[f[1]], self.points[f[2]]]
            if (poly[0][2]>=crop_z and poly[1][2]>=crop_z and poly[2][2]>=crop_z):
                new_faces.append(f)
        self.faces = new_faces
        self.points = [p for p in self.points if p[2] >= crop_z]
        self.basis = [b for b in self.basis if b.origin[2] >= crop_z]
        return
