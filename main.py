from fibonacci_sampling import EFS
from ichosaedron_sampling import EIS
from thomson_sampling import ETS
from vector_base import vector
from point_base import point

from basis_to_rt import *

# configuration
SAVE = True
part_size = point(700, 700, 1400)
d_photo = 500
method = "fibonacci"
n_points = 50
ichos_tesselation = 2
crop_z = 600

# compute
# takes part_size + d_ph as the sphere radius
radius = part_size + point(1, 1, 1) * d_photo
print("[MAIN] - Sampling ellipsoid...")
if method == "fibonacci":
    ellipsoid = EFS(radius[0], radius[1], radius[2])
    ellipsoid.sample(n_points * 2)
elif method == "ichosaedron":
    ellipsoid = EIS(radius[0], radius[1], radius[2])
    ellipsoid.sample(ichos_tesselation)
elif method == "thomson":
    ellipsoid = ETS(radius[0], radius[1], radius[2])
    ellipsoid.sample(n_points * 2)
else:
    raise ValueError("Unknown sampling method: " + method)
print("[MAIN] - Sampled %d points." %(len(ellipsoid.basis)))

# plot
print("[MAIN] - Plotting sampled points...")
ellipsoid.plot(False, True, False, crop_z)

# adjust positions
print("[MAIN] - Adjusting positions to robot tool...")
# correct with tcp_tool
max_base_z = 0
min_base_z = 0
ellipsoid.crop_z(crop_z)
for b in ellipsoid.basis:
    # Get maximum and minimum base z values
    if b.origin[2] > max_base_z:
        max_base_z = b.origin[2]
    elif b.origin[2] < min_base_z:
        min_base_z = b.origin[2]

# adjust rotation around z based on z value
for b in ellipsoid.basis:
    b.rotate(-90 * (b.origin[2] - min_base_z) / (max_base_z - min_base_z), 'z')
print("[MAIN] - Adjusted %d points." %(len(ellipsoid.basis)))

# plot
print("[MAIN] - Plotting adjusted points...")
ellipsoid.plot(False, True, False, min_base_z)

# export to RAPID
if SAVE:
    print("[MAIN] - Exporting points to RAPID code...")
    create_RAPID_module(ellipsoid.basis, "PRUEBA")

print("[MAIN] - Finished.")

# -no olvidar que V-Stars recomienda girar 90 grados al menos una vez
# -emilio recomienda girar otra vez
