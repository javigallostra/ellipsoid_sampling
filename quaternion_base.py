from math import sqrt

class quaternion:

    def __init__(self, *args):
        """ Initialize a quaternion.

        Each initialization value corresponds to
        the vector value in the given dimension.
        Can be initialized with a list or tuple.
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
                        elif arg > 1 or arg < 0:
                            raise ValueError("Only scalars between 0 and 1 can initialize a quaternion")
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
                    elif arg > 1 or arg < 0:
                        raise ValueError("Only scalars between 0 and 1 can initialize a quaternion")
                        # @todo: delete object, initialization failed
                    # Valid scalar
                    else:
                        self.values.append(float(arg))
                        
    def __getitem__(self, i):
        """ Get a quaternion value.

        Only non-negative integer indices are allowed.
        If the index is higher than the 3,
        an error will be raised.
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
        """ Compute the quaternion 2-norm."""

        return sqrt(sum([i**2 for i in self.values]))

    def normalize(self):
        """ Normalize the quaternion."""

        n = self.norm()
        for i in range(len(self.values)):
            self.values[i] /= n

    def normalized(self):
        """ Compute the normalized quaternion."""

        n = self.norm()
        return quaternion([i/n for i in self.values])
