from sampling_methods.fibonacci_sampling import EFS
from sampling_methods.ichosaedron_sampling import EIS
from sampling_methods.thomson_sampling import ETS
from base_classes.vector_base import vector
from base_classes.point_base import point

# configuration
size = point(200, 200, 1400)
d_extra = round((1400/2)/0.6, 0)
method = "fibonacci"
n_points = 70
ichos_tesselation = 2

# compute
# takes size + d_extra as the ellipsoid radii
radii = size + point(1, 1, 1) * d_extra
print("[MAIN] - Sampling ellipsoid...")
if method == "fibonacci":
    ellipsoid = EFS(radii[0], radii[1], radii[2])
    ellipsoid.sample(n_points)
elif method == "ichosaedron":
    ellipsoid = EIS(radii[0], radii[1], radii[2])
    ellipsoid.sample(ichos_tesselation)
elif method == "thomson":
    ellipsoid = ETS(radii[0], radii[1], radii[2])
    ellipsoid.sample(n_points * 2)
else:
    raise ValueError("Unknown sampling method: " + method)
print("[MAIN] - Sampled %d points." %(len(ellipsoid.basis)))

# plot
print("[MAIN] - Plotting sampled points...")
ellipsoid.plot(True, False, False)
