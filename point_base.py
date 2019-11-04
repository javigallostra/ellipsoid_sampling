class point:

    def __init__(self, x=0, y=0, z=0, cs="cartesian"):
        self.cs = cs
        if self.cs == "cartesian":
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        elif self.cs == "spherical":
            self.lat = float(x)
            self.long = float(y)
        else:
            return
            # @todo: implement errors

    def __getitem__(self, i):
        """ Get a point value.

        Only non-negative integer indices are allowed.
        If the index is higher than the 2,
        an error will be raised.
        """
        
        # Invalid argument type
        if type(i) != type(0):
            raise TypeError("Only integer indices are allowed")
        # Invalid argument value 1
        elif i < 0:
            raise IndexError("Point class does not support negative indices")
        # Invalid argument value 2
        elif i > 2:
            raise IndexError("Point class does not support indices greater than 2")
        # Get point value
        values = [self.x, self.y, self.z]
        result = values[i]
        # Return
        return result

    def __add__(self, p2):

        if type(p2) != type(self):
            raise TypeError("Unsupported operand types: %s and %s" %(type(self), type(p2)))
        elif p2.cs != self.cs:
            raise TypeError("Cannot operate with points from different coordinate systems")
        elif self.cs == "spherical":
            result = point(self.lat + p2.lat, self.long + p2.long, None, "spherical")
        else:
            result = point(self.x + p2.x, self.y + p2.y, self.z + p2.z)
        return result

    def __iadd__(self, p2):

        if type(p2) != type(self):
            raise TypeError("Unsupported operand types: %s and %s" %(type(self), type(p2)))
        elif p2.cs != self.cs:
            raise TypeError("Cannot operate with points from different coordinate systems")
        elif self.cs == "spherical":
            self.lat += p2.lat
            self.long += p2.long
        else:
            self.x += p2.x
            self.y += p2.y
            self.z += p2.z
        return self

    def __sub__(self, p2):

        if type(p2) != type(self):
            raise TypeError("Unsupported operand types: %s and %s" %(type(self), type(p2)))
        elif p2.cs != self.cs:
            raise TypeError("Cannot operate with points from different coordinate systems")
        elif self.cs == "spherical":
            result = point(self.lat - p2.lat, self.long - p2.long, None, "spherical")
        else:
            result = point(self.x - p2.x, self.y - p2.y, self.z - p2.z)
        return result

    def __isub__(self, p2):

        if type(p2) != type(self):
            raise TypeError("Unsupported operand types: %s and %s" %(type(self), type(p2)))
        elif p2.cs != self.cs:
            raise TypeError("Cannot operate with points from different coordinate systems")
        elif self.cs == "spherical":
            self.lat -= p2.lat
            self.long -= p2.long
        else:
            self.x -= p2.x
            self.y -= p2.y
            self.z -= p2.z
        return self

    def __mul__(self, b):
        """ Multiply by a scalar."""
                 
        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Scalar multiplication
        elif self.cs == "spherical":
            result = point(self.lat * b, self.long * b, None, "spherical")
        else:
            result = point(self.x * b, self.y * b, self.z * b)
        # Return
        return result

    def __imul__(self,b):
        """ Multiply itself by a scalar."""
        
        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Scalar multiplication
        elif self.cs == "spherical":
            self.lat *= b
            self.long *= b
        else:
            self.x *= b
            self.y *= b
            self.z *= b
        # Return
        return self

    def __truediv__(self, b):
        """ Divide by a scalar."""
        
        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Invalid value
        elif b == 0:
            raise ValueError("Cannot divide by zero")
        # Scalar division
        elif self.cs == "spherical":
            result = point(self.lat / b, self.long / b, None, "spherical")
        else:
            result = point(self.x / b, self.y / b, self.z / b)
        # Return
        return result

    def __idiv__(self, b):
        """ Divide itself by a scalar."""
        
        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Invalid value
        elif b == 0:
            raise ValueError("Cannot divide by zero")
        # Scalar division
        elif self.cs == "spherical":
            self.lat /= b
            self.long /= b
        else:
            self.x /= b
            self.y /= b
            self.z /= b
        # Return
        return self

    def __eq__(self, p2):
        """ Compare for equality."""

        if type(p2) != type(self):
            return False
        elif p2.cs != self.cs:
            return False
        elif self.cs == "spherical":
            if self.lat == p2.lat and self.long == p2.long:
                return True
            else:
                return False
        else:
            if self.x == p2.x and self.y == p2.y and self.z == p2.z:
                return True
            else:
                return False
            
    def __str__(self):

        if self.cs == "spherical":
            return str([self.lat, self.long])
        else:
            return str([self.x, self.y, self.z])

    def __repr__(self):
        
        return self.__str__()

    def count(self, b):

        # Invalid argument type
        if type(b) not in [type(0), type(0.0)]:
            raise TypeError("Only scalars are allowed")
        # Count occurrences
        elif self.cs == "spherical":
            count = int(self.lat == b) + int(self.long == b)
        else:
            count = int(self.x == b) + int(self.y == b) + (self.z == b)
        # Return
        return count
