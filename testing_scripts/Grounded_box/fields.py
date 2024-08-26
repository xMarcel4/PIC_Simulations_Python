# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:50:59 2024

@author: marce
"""
"""
In the tutorial/book the focus is on creating multi dimensional arrays in c++ on the heap

1. 

"""
import numpy as np
from vec3 import Vec3

class Field:
    def __init__(self, ni, nj, nk):
        self.ni, self.nj, self.nk = ni, nj, nk
        # Allocate a 3D numpy array for data
        self.data = np.zeros((ni, nj, nk))

    def __getitem__(self, index):
        return self.data[index]

    def __call__(self, i, j, k):
        return self.data[i, j, k]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __iadd__(self, other):
        self.data += other.data
        return self

    def __imul__(self, scalar):
        self.data *= scalar
        return self

    def __truediv__(self, other):
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.true_divide(self.data, other.data)
            result[~np.isfinite(result)] = 0  # -inf inf NaN will be set to 0
        self.data = result
        return self

    def __repr__(self):
        return str(self.data)

    def clear(self):
        self.data.fill(0)

    def scatter(self, lc, value):
        i, di = int(lc[0]), lc[0] - int(lc[0])
        j, dj = int(lc[1]), lc[1] - int(lc[1])
        k, dk = int(lc[2]), lc[2] - int(lc[2])

        self.data[i, j, k] += value * (1 - di) * (1 - dj) * (1 - dk)
        self.data[i + 1, j, k] += value * di * (1 - dj) * (1 - dk)
        self.data[i + 1, j + 1, k] += value * di * dj * (1 - dk)
        self.data[i, j + 1, k] += value * (1 - di) * dj * (1 - dk)
        self.data[i, j, k + 1] += value * (1 - di) * (1 - dj) * dk
        self.data[i + 1, j, k + 1] += value * di * (1 - dj) * dk
        self.data[i + 1, j + 1, k + 1] += value * di * dj * dk
        self.data[i, j + 1, k + 1] += value * (1 - di) * dj * dk

    def gather(self, lc):
        i, di = int(lc[0]), lc[0] - int(lc[0])
        j, dj = int(lc[1]), lc[1] - int(lc[1])
        k, dk = int(lc[2]), lc[2] - int(lc[2])

        val = (self.data[i, j, k] * (1 - di) * (1 - dj) * (1 - dk) +
               self.data[i + 1, j, k] * di * (1 - dj) * (1 - dk) +
               self.data[i + 1, j + 1, k] * di * dj * (1 - dk) +
               self.data[i, j + 1, k] * (1 - di) * dj * (1 - dk) +
               self.data[i, j, k + 1] * (1 - di) * (1 - dj) * dk +
               self.data[i + 1, j, k + 1] * di * (1 - dj) * dk +
               self.data[i + 1, j + 1, k + 1] * di * dj * dk +
               self.data[i, j + 1, k + 1] * (1 - di) * dj * dk)
        
        return val

    def __str__(self):
        return str(self.data)


