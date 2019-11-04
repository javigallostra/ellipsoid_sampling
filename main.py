from fibonacci_sampling import EFS
from ichosaedron_sampling import EIS
from thomson_sampling import ETS
from vector_base import vector
from point_base import point

from basis_to_rt import *

# configuration
tcp_tool = vector(0,0,0)#(2000,0,0) #x es -z, y es x, z es -y
part_size = point(425, 425, 1400)
d_photo = 0
method = "fibonacci"
n_points = 40
ichos_tesselation = 2
crop_z = 0

# compute
# takes part_size + d_ph as the sphere radius
radius = part_size + point(1, 1, 1) * d_photo
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

# plot
ellipsoid.plot(False, True, False, None, True)

# correct with tcp_tool
ellipsoid.crop_z(crop_z)
for b in ellipsoid.basis:
    b.translate(tcp_tool[1], 'x')
    b.translate(-tcp_tool[2], 'y')
    b.translate(-tcp_tool[0], 'z')

# plot
ellipsoid.plot(False, True, False, -1, True)

# Jugamos con tDynamo:
#   -mover las bases a lo largo de su eje z la distancia tDynamo_z
#   -orientación de la herramienta:
#       -a partir de una z determinada debería girar 180 grados
#       -no olvidar que V-Stars recomienda girar 90 grados al menos una vez

# -90 y, 90 z
