# ellipsoid sampling

## description

This small project aims at providing several methods for uniformly sampling the surface of an arbitrary ellipsoid.

The definition of "uniformly" may vary between different methods, although the overall goal is to distribute
as equally as possible an arbitrary number of points on the surface of an arbitrary ellipsoid.

## methods

The following list shows the implemented and WIP methods.

- [x] **Ichosaedron** sample a unit sphere using the vertices of a regular ichosaedron, subdivide its edges a defined number of times, and project the resulting points onto the ellipsoid
- [ ] **Fibonacci Spiral** sample a unit sphere using a Fibonacci Spiral and an arbitrary number of points, and project the resulting points onto the ellipsoid

## sample results