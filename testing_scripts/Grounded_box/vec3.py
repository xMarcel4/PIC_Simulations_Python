# Example usage of the class with operations
import numpy as np
# from numba import jit
class Vec3:
    def __init__(self, *args):
        if len(args) == 3:
            self.d = list(args)
        elif len(args) == 1 and isinstance(args[0], list) and len(args[0]) == 3:
            self.d = args[0]
        else:
            self.d = [0, 0, 0]
        

    def __getitem__(self, i):
        return self.d[i]

    def __setitem__(self, i, value):
        self.d[i] = value

    def __call__(self, i):
        return self.d[i]

    def __eq__(self, other):
        return self.d == other.d

    def __iadd__(self, other):
        if isinstance(other, Vec3):
            self.d += other.d
        elif isinstance(other, np.ndarray):
            self.d += other
        else:
            raise TypeError("Unsupported type for addition")
        return self
    
    def __isub__(self, other):
        if isinstance(other, Vec3):
            self.d -= other.d
        elif isinstance(other, (list, tuple, np.ndarray)) and len(other) == 3:
            self.d -= other
        else:
            raise TypeError("Unsupported type for subtraction")
        return self

    def __imul__(self, other):
        if isinstance(other, Vec3):
            for i in range(3):
                self.d[i] *= other(i)
        else:  # Assume it's a scalar
            for i in range(3):
                self.d[i] *= other
        return self

    def __add__(self, other):
        return Vec3([self.d[i] + other(i) for i in range(3)])

    def __sub__(self, other):
        return Vec3([self.d[i] - other(i) for i in range(3)])

    def __mul__(self, other):
        if isinstance(other, Vec3):
            return Vec3([self.d[i] * other(i) for i in range(3)])
        else:  # Assume it's a scalar
            return Vec3([self.d[i] * other for i in range(3)])

    def __rmul__(self, other):
        # This handles cases like 2 * Vec3
        return self.__mul__(other)

    def __truediv__(self, other):
        return Vec3([self.d[i] / other(i) for i in range(3)])

    def __repr__(self):
        return f"{self.d[0]} {self.d[1]} {self.d[2]}"

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return all(self.d[i] == other[i] for i in range(3))
    
    def copy(self):
        return Vec3(self.d[:])
    
# Aliases for Vec3 with specific types
#Double3 = Vec3
#Int3 = Vec3
    
    



