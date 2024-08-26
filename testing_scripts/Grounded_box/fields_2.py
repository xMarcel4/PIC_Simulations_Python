import numpy as np

class Field:
    def __init__(self, ni, nj, nk):
        print("Field Constructor!")
        self.ni, self.nj, self.nk = ni, nj, nk
        self.data = np.zeros((ni, nj, nk))  # allocate memory and initialize with zeros

    def __del__(self):
        print("Field Destructor!")
        del self.data

    def __copy__(self):
        print("Field Copy Constructor!")
        new_field = Field(self.ni, self.nj, self.nk)
        np.copyto(new_field.data, self.data)
        return new_field

    def __deepcopy__(self, memodict={}):
        print("Field Copy Constructor!")
        new_field = Field(self.ni, self.nj, self.nk)
        new_field.data = np.copy(self.data)
        return new_field

    def move(self):
        print("Field Move Constructor!")
        new_field = Field(self.ni, self.nj, self.nk)
        new_field.data, self.data = self.data, None
        return new_field

    def r_at(self, i, j, k):
        return self.data[i, j, k]

    def w_at(self, i, j, k, value):
        self.data[i, j, k] = value

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __call__(self, i, j, k):
        return self.data[i, j, k]

    def __repr__(self):
        return str(self.data)

    def __str__(self):
        return self.__repr__()

    def __iadd__(self, other):
        self.data += other.data
        return self

    def __imul__(self, scalar):
        self.data *= scalar
        return self

    def __itruediv__(self, other):
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.true_divide(self.data, other.data)
            result[~np.isfinite(result)] = 0  # Set -inf, inf, NaN to 0
        self.data = result
        return self

    def scatter(self, lc, value):
        # Ensure the logical coordinate is within bounds
        if lc[0] < 0 or lc[0] > (self.ni - 1) or \
           lc[1] < 0 or lc[1] > (self.nj - 1) or \
           lc[2] < 0 or lc[2] > (self.nk - 1):
            return

        # Compute the cell index and fractional distances
        i, di = int(lc[0]), lc[0] - int(lc[0])
        j, dj = int(lc[1]), lc[1] - int(lc[1])
        k, dk = int(lc[2]), lc[2] - int(lc[2])

        # Deposit fractional values to the 8 surrounding nodes
        self.data[i, j, k] += value * (1 - di) * (1 - dj) * (1 - dk)
        self.data[i + 1, j, k] += value * di * (1 - dj) * (1 - dk)
        self.data[i, j + 1, k] += value * (1 - di) * dj * (1 - dk)
        self.data[i + 1, j + 1, k] += value * di * dj * (1 - dk)
        self.data[i, j, k + 1] += value * (1 - di) * (1 - dj) * dk
        self.data[i + 1, j, k + 1] += value * di * (1 - dj) * dk
        self.data[i, j + 1, k + 1] += value * (1 - di) * dj * dk
        self.data[i + 1, j + 1, k + 1] += value * di * dj * dk

    def gather(self, lc):
        # Compute the cell index and fractional distances
        i, di = int(lc[0]), lc[0] - int(lc[0])
        j, dj = int(lc[1]), lc[1] - int(lc[1])
        k, dk = int(lc[2]), lc[2] - int(lc[2])

        # Interpolate data onto the logical coordinate
        val = self.r_at(i, j, k) * (1 - di) * (1 - dj) * (1 - dk) + \
              self.r_at(i + 1, j, k) * di * (1 - dj) * (1 - dk) + \
              self.r_at(i, j + 1, k) * (1 - di) * dj * (1 - dk) + \
              self.r_at(i + 1, j + 1, k) * di * dj * (1 - dk) + \
              self.r_at(i, j, k + 1) * (1 - di) * (1 - dj) * dk + \
              self.r_at(i + 1, j, k + 1) * di * (1 - dj) * dk + \
              self.r_at(i, j + 1, k + 1) * (1 - di) * dj * dk + \
              self.r_at(i + 1, j + 1, k + 1) * di * dj * dk

        return val

    def __deepcopy__(self, memodict={}):
        new_field = Field(self.ni, self.nj, self.nk)
        new_field.data = np.copy(self.data)
        return new_field

# Overloading the * operator for scalar multiplication
def multiply_field(scalar, field):
    new_field = Field(field.ni, field.nj, field.nk)
    new_field.data = field.data * scalar
    return new_field

# Overloading the << operator to mimic the stream output in C++
def field_output(field):
    output = "\n"
    for k in range(field.nk):
        for j in range(field.nj):
            for i in range(field.ni):
                output += f"{field(i, j, k)} "
        output += "\n"
    return output
