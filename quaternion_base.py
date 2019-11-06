from math import sqrt

class quaternion:
    """ Basic quaternion class.

    Implements representation methods, length, getter,
    iterators, normalization methods and a method to
    create a quaternion from a rotation matrix.
    """

    def __init__(self, *args):
        """ Set the initial parameters.

        Each initialization value corresponds to
        the vector value in the given dimension.
        Can be initialized with a list or tuple.

        Return nothing.
        """
        
        self.values = []
        self.iter_n = 0
        # No arguments: default values
        if len(args) == 0:
                self.values = [0.0, 0.0, 0.0, 1.0]
        # Single argument: quaternion, list or tuple
        elif len(args) == 1:
            # Invalid argument type
            if type(args[0]) not in [type(self), type([]), type(())]:
                raise TypeError("Only quaternions or a scalar list/tuple can initialize a quaternion")
                # @todo: delete object, initialization failed
            # Quaternion
            elif type(args[0]) == type(self):
                for element in args[0]:
                    self.values.append(element)
            # List or tuple
            else:
                # Invalid number of arguments
                if len(args[0]) != 4:
                    raise IndexError("Four scalars are required to initialise a quaternion")
                    # @todo: delete object, initialization failed
                else:
                    for arg in args[0]:
                        # Invalid argument type
                        if type(arg) not in [type(0), type(0.0)]:
                            raise TypeError("Only a scalar list/tuple can initialize a quaternion")
                            # @todo: delete object, initialization failed
                        # Invalid scalar
                        elif abs(arg) > 1:
                            raise ValueError("Only scalars between -1 and 1 can initialize a quaternion")
                            # @todo: delete object, initialization failed
                        # Valid scalar
                        else:
                            self.values.append(float(arg))
        # Multiple arguments: scalars
        else:
            # Invalid number of arguments
            if len(args) != 4:
                raise IndexError("Four scalars are required to initialise a quaternion")
                # @todo: delete object, initialization failed
            else:
                for arg in args:
                    # Invalid argument type
                    if type(arg) not in [type(0), type(0.0)]:
                        raise TypeError("Only scalars can initialize a quaternion")
                        # @todo: delete object, initialization failed
                    # Invalid scalar
                    elif abs(arg) > 1:
                        raise ValueError("Only scalars between -1 and 1 can initialize a quaternion")
                        # @todo: delete object, initialization failed
                    # Valid scalar
                    else:
                        self.values.append(float(arg))
        # Return
        return
                        
    def __getitem__(self, i):
        """ Get a quaternion value.

        Only non-negative integer indices are allowed.
        If the index is higher than the 3,
        an error will be raised.

        Return the requested value.
        """
        
        # Invalid argument type
        if type(i) != type(0):
            raise TypeError("Only integer indices are allowed")
        # Invalid argument value 1
        elif i < 0:
            raise IndexError("Quaternion class does not support negative indices")
        # Invalid argument value 2
        elif i > 3:
            raise IndexError("Quaternion class does not support indices greater than 3")
        # Get quaternion value
        result = self.values[i]
        # Return
        return result

    def __len__(self):
        """ Return the length of the quaternion."""

        return len(self.values)

    def __str__(self):
        """ Return the string representation of the quaternion."""
        
        return str(self.values)

    def __repr__(self):
        """ Return the representation of the quaternion."""
        
        return str(self.values)

    def __iter__(self):
        """ Return itself as an iterator."""
        
        self.iter_n = 0
        return self

    def __next__(self):
        """ Return the next element of the quaternion iterator."""
        
        if self.iter_n < len(self.values):
            element = self.values[self.iter_n]
            self.iter_n += 1
            return element
        else:
            raise StopIteration

    def norm(self):
        """ Return the quaternion 2-norm."""

        return sqrt(sum([i**2 for i in self.values]))

    def normalize(self):
        """ Normalize the quaternion.

        Return nothing.
        """

        n = self.norm()
        for i in range(len(self.values)):
            self.values[i] /= n
        return

    def normalized(self):
        """ Compute the normalized quaternion.

        Return a new quaternion.
        """

        n = self.norm()
        return quaternion([i/n for i in self.values])

    def from_matrix(self, matrix):
        """ Compute a quaternion from an orthonormal rotation matrix.

        Return the new quaternion.
        """
        
        m = matrix
        # Compute trace
        trace = m[0][0] + m[1][1] + m[2][2]
        # Compute x,y,z,w
        if trace > 0: 
            S = sqrt(trace + 1) * 2
            w = 0.25 * S
            x = (m[2][1] - m[1][2]) / S
            y = (m[0][2] - m[2][0]) / S
            z = (m[1][0] - m[0][1]) / S
        elif m[0][0] > m[1][1] and m[0][0] > m[2][2]: 
            S = sqrt(1 + m[0][0] - m[1][1] - m[2][2]) * 2 
            w = (m[2][1] - m[1][2]) / S
            x = 0.25 * S
            y = (m[0][1] + m[1][0]) / S
            z = (m[0][2] + m[2][0]) / S 
        elif m[1][1] > m[2][2]: 
            S = sqrt(1 + m[1][1] - m[0][0] - m[2][2]) * 2
            w = (m[0][2] - m[2][0]) / S
            x = (m[0][1] + m[1][0]) / S 
            y = 0.25 * S
            z = (m[1][2] + m[2][1]) / S
        else:
            S = sqrt(1 + m[2][2] - m[0][0] - m[1][1]) * 2
            w = (m[1][0] - m[0][1]) / S
            x = (m[0][2] + m[2][0]) / S
            y = (m[1][2] + m[2][1]) / S
            z = 0.25 * S
        # Get result
        result = quaternion(x,y,z,w)
        # Return
        return result
