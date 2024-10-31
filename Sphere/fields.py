import numpy as np

class Field:        
    def __init__(self, ni, nj, nk, components=1):
        #print("Field Constructor!")
        self.ni, self.nj, self.nk = ni, nj, nk
        self.components = components
        if self.components == 1:
            self.data = np.zeros((ni, nj, nk))
        else:
            self.data = np.zeros((ni, nj, nk, components))

    

    def __del__(self):
        #print("Field Destructor!")
        del self.data

    def __copy__(self):
        print("Field Copy Constructor!")
        new_field = Field(self.ni, self.nj, self.nk, self.components)
        np.copyto(new_field.data, self.data)
        return new_field

    def __deepcopy__(self, memodict={}):
        print("Field Copy Constructor!")
        new_field = Field(self.ni, self.nj, self.nk, self.components)
        new_field.data = np.copy(self.data)
        return new_field

    def move(self):
        print("Field Move Constructor!")
        new_field = Field(self.ni, self.nj, self.nk, self.components)
        new_field.data, self.data = self.data, None
        return new_field

    def r_at(self, i, j, k, c=None):
        if self.components == 1 or c is None:
            return self.data[i, j, k]
        else:
            return self.data[i, j, k, c]

    def w_at(self, i, j, k, value, component=None):
        if self.components == 1:
            self.data[i, j, k] = value
        else:
            self.data[i, j, k, component] = value

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __call__(self, i, j, k, component=None):
        if self.components == 1:
            return self.data[i, j, k]
        else:
            return self.data[i, j, k, component]

    def __repr__(self):
        return str(self.data)

    def __str__(self):
        return self.__repr__()

    def __iadd__(self, other):
        if isinstance(other, Field):
            self.data += other.data  # Element-wise addition
        else:
            raise ValueError("Cannot add non-Field instance to Field")
        return self

    def __imul__(self, scalar):
        self.data *= scalar
        return self
    
    def __rmul__(self, scalar):
        return self.__mul__(scalar)
    

    def __itruediv__(self, other):
     # print(f"Dividing: self.type={type(self)}, other.type={type(other)}")
     # print(f"self.data.shape: {self.data.shape}, other.data.shape: {other.data.shape}")
     
     # Check if the other instance is a Field and has matching dimensions
     if not isinstance(other, Field):
         print(f"Error: Expected Field, but got {type(other)} instead.")
         raise ValueError("Cannot divide by a non-Field instance")

     if self.data.shape != other.data.shape:
         print(f"Error: Shape mismatch. self.data.shape={self.data.shape}, other.data.shape={other.data.shape}")
         raise ValueError("Field shapes do not match for division")

     try:
         # Perform element-wise division
         self.data = np.where(other.data != 0, self.data / other.data, 0)
     except Exception as e:
         print(f"Error during division: {e}")
         raise
     
     return self
    
    
    # @nijt
    def scatter(self, lc, value):
        # Ensure the logical coordinate is within bounds
        if lc[0] < 0 or lc[0] > (self.ni - 1) or \
           lc[1] < 0 or lc[1] > (self.nj - 1) or \
           lc[2] < 0 or lc[2] > (self.nk - 1):
            return
    
       # Compute the cell index and the fractional distances
        i = int(lc[0])
        di = lc[0] - i
        j = int(lc[1])
        dj = lc[1] - j
        k = int(lc[2])
        dk = lc[2] - k
    
        if self.components == 1:
            # Deposit fractional values to the 8 surrounding nodes for scalar fields
            self.data[i, j, k] += value * (1 - di) * (1 - dj) * (1 - dk)
            self.data[i + 1, j, k] += value * di * (1 - dj) * (1 - dk)
            self.data[i, j + 1, k] += value * (1 - di) * dj * (1 - dk)
            self.data[i + 1, j + 1, k] += value * di * dj * (1 - dk)
            self.data[i, j, k + 1] += value * (1 - di) * (1 - dj) * dk
            self.data[i + 1, j, k + 1] += value * di * (1 - dj) * dk
            self.data[i, j + 1, k + 1] += value * (1 - di) * dj * dk
            self.data[i + 1, j + 1, k + 1] += value * di * dj * dk
    
        elif self.components > 1:
            print("yes")
            # Deposit fractional values to the 8 surrounding nodes for vector fields
            for comp in range(self.components):
                self.data[i, j, k, comp] += value[comp] * (1 - di) * (1 - dj) * (1 - dk)
                self.data[i + 1, j, k, comp] += value[comp] * di * (1 - dj) * (1 - dk)
                self.data[i, j + 1, k, comp] += value[comp] * (1 - di) * dj * (1 - dk)
                self.data[i + 1, j + 1, k, comp] += value[comp] * di * dj * (1 - dk)
                self.data[i, j, k + 1, comp] += value[comp] * (1 - di) * (1 - dj) * dk
                self.data[i + 1, j, k + 1, comp] += value[comp] * di * (1 - dj) * dk
                self.data[i, j + 1, k + 1, comp] += value[comp] * (1 - di) * dj * dk
                self.data[i + 1, j + 1, k + 1, comp] += value[comp] * di * dj * dk
    
        else:
            raise NotImplementedError("Scatter operation is not implemented for this field type.")

    # @nijt
    def gather(self, lc):
        
        i = int(lc[0])
        di = lc[0] - i
        j = int(lc[1])
        dj = lc[1] - j
        k = int(lc[2])
        dk = lc[2] - k
        if self.components == 1:
            # Scalar field gather operation
            val = (
                self.r_at(i, j, k) * (1 - di) * (1 - dj) * (1 - dk) +
                self.r_at(i + 1, j, k) * di * (1 - dj) * (1 - dk) +
                self.r_at(i, j + 1, k) * (1 - di) * dj * (1 - dk) +
                self.r_at(i + 1, j + 1, k) * di * dj * (1 - dk) +
                self.r_at(i, j, k + 1) * (1 - di) * (1 - dj) * dk +
                self.r_at(i + 1, j, k + 1) * di * (1 - dj) * dk +
                self.r_at(i, j + 1, k + 1) * (1 - di) * dj * dk +
                self.r_at(i + 1, j + 1, k + 1) * di * dj * dk
            )
            return val

        elif self.components > 1:
            # print("Executed components bigger 1")
            # Vector field gather operation
            interpolated = np.zeros(self.components)
            for comp in range(self.components):
                interpolated[comp] = (
                    self.r_at(i, j, k, comp) * (1 - di) * (1 - dj) * (1 - dk) +
                    self.r_at(i + 1, j, k, comp) * di * (1 - dj) * (1 - dk) +
                    self.r_at(i, j + 1, k, comp) * (1 - di) * dj * (1 - dk) +
                    self.r_at(i + 1, j + 1, k, comp) * di * dj * (1 - dk) +
                    self.r_at(i, j, k + 1, comp) * (1 - di) * (1 - dj) * dk +
                    self.r_at(i + 1, j, k + 1, comp) * di * (1 - dj) * dk +
                    self.r_at(i, j + 1, k + 1, comp) * (1 - di) * dj * dk +
                    self.r_at(i + 1, j + 1, k + 1, comp) * di * dj * dk
                )
            return interpolated

        else:
            raise NotImplementedError("Gather operation is not implemented for this field type.")

    def __deepcopy__(self, memodict={}):
        new_field = Field(self.ni, self.nj, self.nk, self.components)
        new_field.data = np.copy(self.data)
        return new_field

    def __mul__(self, scalar):
        if isinstance(scalar, (int, float)):
            result = Field(*self.data.shape)
            result.data = self.data * scalar
            return result
        else:
            raise TypeError("Scalar multiplication only supports int or float types")
    
# Overloading the << operator to mimic the stream output in C++
def field_output(field):
    output = "\n"
    for k in range(field.nk):
        for j in range(field.nj):
            for i in range(field.ni):
                if field.components == 1:
                    output += f"{field(i, j, k)} "
                else:
                    output += f"({field(i, j, k, 0)}, {field(i, j, k, 1)}, {field(i, j, k, 2)}) "
        output += "\n"
    return output
