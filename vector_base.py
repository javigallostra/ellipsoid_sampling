class vector:
    """ Vector class."""

    def __init__(self, *args):
        """ Initialize a vector.

        Each initialization value corresponds to
        the vector value in the given dimension.
        Can be initialized with a list or tuple.
        """
        
        self.values = []
        if len(args) == 1:
            if type(args[0]) not in [type(0), type(0.0),type([]), type(())]:
                raise TypeError("Only scalars or a scalar list/tuple can initialize a vector")
                # @todo: delete object, initialization failed
            elif type(args[0]) in [type(0), type(0.0)]:
                self.values.append(float(args[0]))
            else:
                for arg in args[0]:
                    if type(arg) not in [type(0), type(0.0)]:
                        raise TypeError("Only scalars or a scalar list/tuple can initialize a vector")
                        # @todo: delete object, initialization failed
                    else:
                        self.values.append(float(arg))
        else:
            for arg in args:
                if type(arg) not in [type(0), type(0.0)]:
                    raise TypeError("Only scalars or a scalar list/tuple can initialize a vector")
                    # @todo: delete object, initialization failed
                else:
                    self.values.append(float(arg))

    def __getitem__(self, i):
        """ Get a vector value.

        Only non-negative integer indices are allowed.
        If the index is higher than the vector dimension,
        the returned value will be 0 but no error will be raised.
        """
        
        # 1 - Check errors
        # TypeError: only integers
        if type(i) != type(0):
            raise TypeError("Only integer indices are allowed")
        # IndexError: no negative indices
        elif i < 0:
            raise IndexError("vector class does not support negative indices")
        # 2 - Compute value
        if i >= len(self.values):
            result = 0.0
        else:
            result = self.values[i]
        # 3 - Return
        return result

    def __mul__(self, b):
        """ Multiply by a scalar or another vector (dot product)."""
        
        # 1 - Check errors
        # TypeErrors: integers, floats and other vectors
        if type(b) not in [type(0), type(0.0), type(self)]:
            raise TypeError("Only scalars or vectors are allowed")
        # 2 - Compute value
        # Scalar multiplication
        if type(b) != type(self):
            new_values = []
            for i in range(len(self.values)):
                new_values.append(self.values[i] * b)
            result = vector(new_values)
        # Dot product
        else:
            result = 0.0
            for i in range(min(len(self.values), len(b.values))):
                result += self.values[i] * b[i]
        # 3 - Return
        return result

    def __imult__(self,b):
        """ Multiply itself by a scalar."""
        
        # 1 - Check errors
        # TypeErrors: integers, floats and other vectors
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # 2 - Compute value
        # Scalar multiplication
        for i in range(len(self.values)):
            self.values[i] *= b

    def __len__(self):

        return len(self.values)

    def __str__(self):
        
        return str(self.values)

    def __repr__(self):
        
        return str(self.values)
            
    # queremos: producto por escalar a*b
    # queremos: producto por vector (dot product) a*b
    # queremos: producto vectorial (cross product) a.cross(b)
    # queremos: alguna manera de decir que es un vector unitario
