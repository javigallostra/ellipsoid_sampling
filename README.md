# ellipsoid sampling

## description

This small project aims at providing several methods for uniformly sampling the surface of an arbitrary ellipsoid.

The definition of "uniformly" may vary between different methods, although the overall goal is to distribute
as equally as possible an arbitrary number of points on the surface of an arbitrary ellipsoid.

## methods

The following list shows the implemented and WIP methods.

- [x] **Ichosaedron** sample a unit sphere using the vertices of a regular ichosaedron, subdivide its edges a defined number of times, and project the resulting points onto the ellipsoid
- [x] **Fibonacci Spiral** sample a unit sphere using a Fibonacci Spiral and an arbitrary number of points, and project the resulting points onto the ellipsoid
- [ ] **Thomson Problem** sample a unit sphere by simulating the behavior of an arbitrary number of atoms lying on the surface of the sphere.

## sample results
Ellipsoid with "radii" (1,1,3) subdivided twice and tesselated.

![Go to image.][imich]

[imich]: https://github.com/javigallostra/ellipsoid_sampling/blob/master/sample_res_ich.png "Ichosaedron sampling and tesellation"

Ellipsoid with "radii" (1,1,3) sampled using the Fibonacci Spiral method with 100 points.

![Go to image.][imfib]

[imfib]: https://github.com/javigallostra/ellipsoid_sampling/blob/master/sample_res_fib.png "Fibonacci Spiral sampling using 100 points"
