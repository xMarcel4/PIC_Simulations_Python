# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:19:57 2024

@author: marce
"""
import numpy as np
class Vec2:
    def __init__(self, *args):
        if len(args) == 2:
            self.d = list(args)
        elif len(args) == 1 and isinstance(args[0], (list, tuple, np.ndarray)) and len(args[0]) == 2:
            self.d = list(args[0])
        else:
            self.d = [0.0, 0.0]

    def __getitem__(self, i):
        return self.d[i]

    def __setitem__(self, i, value):
        self.d[i] = value

    def __call__(self, i):
        return self.d[i]

    def __iadd__(self, other):
        if isinstance(other, Vec2):
            self.d = [self.d[i] + other.d[i] for i in range(2)]
        elif isinstance(other, (list, tuple, np.ndarray)):
            self.d = [self.d[i] + other[i] for i in range(2)]
        else:
            raise TypeError("Unsupported type for addition")
        return self

    def __add__(self, other):
        if isinstance(other, Vec2):
            return Vec2([self.d[i] + other.d[i] for i in range(2)])
        elif isinstance(other, (list, tuple, np.ndarray)):
            return Vec2([self.d[i] + other[i] for i in range(2)])
        else:
            raise TypeError("Unsupported type for addition")

    def __sub__(self, other):
        if isinstance(other, Vec2):
            return Vec2([self.d[i] - other.d[i] for i in range(2)])
        elif isinstance(other, (list, tuple, np.ndarray)):
            return Vec2([self.d[i] - other[i] for i in range(2)])
        else:
            raise TypeError("Unsupported type for subtraction")

    def __mul__(self, other):
        if isinstance(other, Vec2):
            return Vec2([self.d[i] * other.d[i] for i in range(2)])
        elif isinstance(other, (int, float)):  # Scalar multiplication
            return Vec2([self.d[i] * other for i in range(2)])
        elif isinstance(other, (list, tuple, np.ndarray)):  # Element-wise multiplication with array
            return Vec2([self.d[i] * other[i] for i in range(2)])
        else:
            raise TypeError("Unsupported type for multiplication")

    def __rmul__(self, other):
        # This handles cases like 2 * Vec2
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Vec2):
            return Vec2([self.d[i] / other.d[i] for i in range(2)])
        elif isinstance(other, (int, float)):  # Scalar division
            return Vec2([self.d[i] / other for i in range(2)])
        elif isinstance(other, (list, tuple, np.ndarray)):  # Element-wise division with array
            return Vec2([self.d[i] / other[i] for i in range(2)])
        else:
            raise TypeError("Unsupported type for division")

    def dot(self, other):
        if isinstance(other, Vec2):
            return sum(self.d[i] * other.d[i] for i in range(2))
        elif isinstance(other, (list, tuple, np.ndarray)):
            return sum(self.d[i] * other[i] for i in range(2))
        else:
            raise TypeError("Unsupported type for dot product")

    def __repr__(self):
        return f"Vec2({self.d[0]}, {self.d[1]})"

    def __str__(self):
        return f"{self.d[0]} {self.d[1]}"