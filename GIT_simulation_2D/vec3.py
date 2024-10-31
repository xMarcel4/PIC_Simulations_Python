import numpy as np




class Vec3:
    def __init__(self, *args):
        if len(args) == 3:
            self.d = list(args)
        elif len(args) == 1 and isinstance(args[0], (list, tuple, np.ndarray)) and len(args[0]) == 3:
            self.d = list(args[0])
        else:
            self.d = [0.0, 0.0, 0.0]

    def __getitem__(self, i):
        return self.d[i]

    def __setitem__(self, i, value):
        self.d[i] = value

    def __call__(self, i):
        return self.d[i]

    def __iadd__(self, other):
        if isinstance(other, Vec3):
            self.d = [self.d[i] + other.d[i] for i in range(3)]
        elif isinstance(other, (list, tuple, np.ndarray)):
            self.d = [self.d[i] + other[i] for i in range(3)]
        else:
            raise TypeError("Unsupported type for addition")
        return self

    def __add__(self, other):
        if isinstance(other, Vec3):
            return Vec3([self.d[i] + other.d[i] for i in range(3)])
        elif isinstance(other, (list, tuple, np.ndarray)):
            return Vec3([self.d[i] + other[i] for i in range(3)])
        else:
            raise TypeError("Unsupported type for addition")

    def __sub__(self, other):
        if isinstance(other, Vec3):
            return Vec3([self.d[i] - other.d[i] for i in range(3)])
        elif isinstance(other, (list, tuple, np.ndarray)):
            return Vec3([self.d[i] - other[i] for i in range(3)])
        else:
            raise TypeError("Unsupported type for subtraction")

    def __mul__(self, other):
        if isinstance(other, Vec3):
            return Vec3([self.d[i] * other.d[i] for i in range(3)])
        elif isinstance(other, (int, float)):  # Scalar multiplication
            return Vec3([self.d[i] * other for i in range(3)])
        elif isinstance(other, (list, tuple, np.ndarray)):  # Element-wise multiplication with array
            return Vec3([self.d[i] * other[i] for i in range(3)])
        else:
            raise TypeError("Unsupported type for multiplication")

    def __rmul__(self, other):
        # This handles cases like 2 * Vec3
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Vec3):
            return Vec3([self.d[i] / other.d[i] for i in range(3)])
        elif isinstance(other, (int, float)):  # Scalar division
            return Vec3([self.d[i] / other for i in range(3)])
        elif isinstance(other, (list, tuple, np.ndarray)):  # Element-wise division with array
            return Vec3([self.d[i] / other[i] for i in range(3)])
        else:
            raise TypeError("Unsupported type for division")

    def dot(self, other):
        if isinstance(other, Vec3):
            return sum(self.d[i] * other.d[i] for i in range(3))
        elif isinstance(other, (list, tuple, np.ndarray)):
            return sum(self.d[i] * other[i] for i in range(3))
        else:
            raise TypeError("Unsupported type for dot product")


    def __repr__(self):
        return f"Vec3({self.d[0]}, {self.d[1]}, {self.d[2]})"

    def __str__(self):
        return f"{self.d[0]} {self.d[1]} {self.d[2]}"
