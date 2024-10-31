# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:17:39 2024

@author: marce
"""

import numpy as np

class Field2D:
    def __init__(self, ni, nj, components=1):
        #print("Field2D Constructor!")
        self.ni, self.nj = ni, nj
        self.components = components
        if self.components == 1:
            self.data = np.zeros((ni, nj))
        else:
            self.data = np.zeros((ni, nj, components))

    def __del__(self):
        #print("Field2D Destructor!")
        del self.data

    def __copy__(self):
        print("Field2D Copy Constructor!")
        new_field = Field2D(self.ni, self.nj, self.components)
        np.copyto(new_field.data, self.data)
        return new_field

    def __deepcopy__(self, memodict={}):
        print("Field2D Deep Copy Constructor!")
        new_field = Field2D(self.ni, self.nj, self.components)
        new_field.data = np.copy(self.data)
        return new_field

    def move(self):
        print("Field2D Move Constructor!")
        new_field = Field2D(self.ni, self.nj, self.components)
        new_field.data, self.data = self.data, None
        return new_field

    def r_at(self, i, j, c=None):
        if self.components == 1 or c is None:
            return self.data[i, j]
        else:
            return self.data[i, j, c]

    def w_at(self, i, j, value, component=None):
        if self.components == 1:
            self.data[i, j] = value
        else:
            self.data[i, j, component] = value

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __call__(self, i, j, component=None):
        if self.components == 1:
            return self.data[i, j]
        else:
            return self.data[i, j, component]

    def __repr__(self):
        return str(self.data)

    def __str__(self):
        return self.__repr__()

    def __iadd__(self, other):
        if isinstance(other, Field2D):
            self.data += other.data  # Element-wise addition
        else:
            raise ValueError("Cannot add non-Field2D instance to Field2D")
        return self

    def __imul__(self, scalar):
        self.data *= scalar
        return self

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __itruediv__(self, other):
        if not isinstance(other, Field2D):
            raise ValueError("Cannot divide by a non-Field2D instance")

        if self.data.shape != other.data.shape:
            raise ValueError("Field2D shapes do not match for division")

        try:
            self.data = np.where(other.data != 0, self.data / other.data, 0)
        except Exception as e:
            print(f"Error during division: {e}")
            raise

        return self

    def scatter(self, lc, value):
        if lc[0] < 0 or lc[0] > (self.ni - 1) or lc[1] < 0 or lc[1] > (self.nj - 1):
            return
        
        i = int(lc[0])
        di = lc[0] - i
        j = int(lc[1])
        dj = lc[1] - j

        if self.components == 1:
            self.data[i, j] += value * (1 - di) * (1 - dj)
            self.data[i + 1, j] += value * di * (1 - dj)
            self.data[i, j + 1] += value * (1 - di) * dj
            self.data[i + 1, j + 1] += value * di * dj
        elif self.components > 1:
            for comp in range(self.components):
                self.data[i, j, comp] += value[comp] * (1 - di) * (1 - dj)
                self.data[i + 1, j, comp] += value[comp] * di * (1 - dj)
                self.data[i, j + 1, comp] += value[comp] * (1 - di) * dj
                self.data[i + 1, j + 1, comp] += value[comp] * di * dj
        else:
            raise NotImplementedError("Scatter operation is not implemented for this field type.")

    def gather(self, lc):
        i = int(lc[0])
        di = lc[0] - i
        j = int(lc[1])
        dj = lc[1] - j
    
        # Clamp i and j to valid indices (including boundary cases)
        if i < 0:
            i = 0
            di = 0
        elif i >= self.ni - 1:
            i = self.ni - 2
            di = 1
    
        if j < 0:
            j = 0
            dj = 0
        elif j >= self.nj - 1:
            j = self.nj - 2
            dj = 1
    
        if self.components == 1:
            # Scalar field gather operation
            val = (
                self.r_at(i, j) * (1 - di) * (1 - dj) +
                self.r_at(i + 1, j) * di * (1 - dj) +
                self.r_at(i, j + 1) * (1 - di) * dj +
                self.r_at(i + 1, j + 1) * di * dj
            )
            return val
    
        elif self.components > 1:
            # Vector field gather operation
            interpolated = np.zeros(self.components)
            for comp in range(self.components):
                interpolated[comp] = (
                    self.r_at(i, j, comp) * (1 - di) * (1 - dj) +
                    self.r_at(i + 1, j, comp) * di * (1 - dj) +
                    self.r_at(i, j + 1, comp) * (1 - di) * dj +
                    self.r_at(i + 1, j + 1, comp) * di * dj
                )
            return interpolated
    
        else:
            raise NotImplementedError("Gather operation is not implemented for this field type.")
    
    

    def __mul__(self, scalar):
        if isinstance(scalar, (int, float)):
            result = Field2D(*self.data.shape)
            result.data = self.data * scalar
            return result
        else:
            raise TypeError("Scalar multiplication only supports int or float types")

# Overloading the << operator to mimic the stream output in C++
def field_output(field):
    output = "\n"
    for j in range(field.nj):
        for i in range(field.ni):
            if field.components == 1:
                output += f"{field(i, j)} "
            else:
                output += f"({field(i, j, 0)}, {field(i, j, 1)}) "
        output += "\n"
    return output
