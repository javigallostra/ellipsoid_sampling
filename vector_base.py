from math import sqrt
from point_base import point

class vector:
    """ Vector class."""

    def __init__(self, *args):
        """ Initialize a vector.

        Each initialization value corresponds to
        the vector value in the given dimension.
        Can be initialized with a list or tuple.
        """
        
        self.values = []
        self.iter_n = 0
        # Single argument: vector, scalar, list or tuple
        if len(args) == 1:
            # Invalid argument type
            if type(args[0]) not in [type(self), type(point()), type(0), type(0.0),type([]), type(())]:
                raise TypeError("Only vectors, points, scalars or a scalar list/tuple can initialize a vector")
                # @todo: delete object, initialization failed
            # Vector or point
            elif type(args[0]) == type(self) or type(args[0]) == type(point):
                for element in args[0]:
                    self.values.append(element)
            # Scalar
            elif type(args[0]) in [type(0), type(0.0)]:
                self.values.append(float(args[0]))
            # List or tuple
            else:
                for arg in args[0]:
                    # Invalid argument type
                    if type(arg) not in [type(0), type(0.0)]:
                        raise TypeError("Only a scalar list/tuple can initialize a vector")
                        # @todo: delete object, initialization failed
                    # Scalar
                    else:
                        self.values.append(float(arg))
        # Multiple arguments: scalars
        else:
            for arg in args:
                # Invalid argument type
                if type(arg) not in [type(0), type(0.0)]:
                    raise TypeError("Only scalars can initialize a vector")
                    # @todo: delete object, initialization failed
                # Scalar
                else:
                    self.values.append(float(arg))

    def __getitem__(self, i):
        """ Get a vector value.

        Only non-negative integer indices are allowed.
        If the index is higher than the vector dimension,
        the returned value will be 0 but no error will be raised.
        """
        
        # Invalid argument type
        if type(i) != type(0):
            raise TypeError("Only integer indices are allowed")
        # Invalid argument value
        elif i < 0:
            raise IndexError("vector class does not support negative indices")
        # Get vector value
        if i >= len(self.values):
            result = 0.0
        else:
            result = self.values[i]
        # Return
        return result

    def __mul__(self, b):
        """ Multiply by a scalar or by another vector (dot product)."""
        
        # Invalid argument type
        if type(b) not in [type(0), type(0.0), type(self)]:
            raise TypeError("Only scalars or vectors are allowed")
        # Scalar multiplication
        if type(b) != type(self):
            result = vector([i * b for i in self.values])
        # Dot product
        else:
            result = sum([self.values[i] * b[i] for i in range(min(len(self.values), len(b.values)))])
        # Return
        return result

    def __imul__(self,b):
        """ Multiply itself by a scalar."""
        
        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Scalar multiplication
        for i in range(len(self.values)):
            self.values[i] *= b

    def __add__(self, b):
        """ Add a vector to another."""

        # Invalid argument type
        if type(b) != type(self):
            raise TypeError("Only vectors are allowed")
        # Return
        return vector([self.values[i] + b[i] for i in range(max(len(self.values), len(b.values)))])

    def __iadd__(self, b):
        """ Add another vector to itself."""

        # Invalid argument type
        if type(b) != type(self):
            raise TypeError("Only vectors are allowed")
        # Compute addition
        # @todo que pasa cuando dimself != dimb...
        for i in range(max(len(self.values), len(b.values))):
            self.values[i] += b[i]

    def __len__(self):

        return len(self.values)

    def __str__(self):
        
        return str(self.values)

    def __repr__(self):
        
        return str(self.values)

    def __iter__(self):
        
        self.iter_n = 0
        return self

    def __next__(self):
        
        if self.iter_n < len(self.values):
            element = self.values[self.iter_n]
            self.iter_n += 1
            return element
        else:
            raise StopIteration

    def norm(self):
        """ Compute the vector 2-norm."""

        return sqrt(sum([i**2 for i in self.values]))

    def normalize(self):
        """ Normalize the vector."""

        n = self.norm()
        for i in range(len(self.values)):
            self.values[i] /= n

    def normalized(self):
        """ Compute the normalized vector."""

        n = self.norm()
        return vector([i/n for i in self.values])

    def cross(self, b):
        """ Compute the cross product a^b.

        Only works for vectors of dimension 3.
        """

        # Invalid argument type
        if type(b) != type(self):
            raise TypeError("Only vectors are allowed")
        # Invalid dimensions
        elif len(self.values) != 3 or len(b.values) != 3:
            raise TypeError("Only 3 dimensional vectors are allowed")
        # Cross product
        c1 = self.values[1]*b.values[2] - b.values[1]*self.values[2]
        c2 = -self.values[0]*b.values[2] + b.values[0]*self.values[2]
        c3 = self.values[0]*b.values[1] - b.values[0]*self.values[1]
        # Return
        return vector(c1, c2, c3)

    # queremos: alguna manera de decir que es un vector unitario
