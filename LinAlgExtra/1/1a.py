class ComplexNumber:
    def __init__(self, real_part, imaginary_part):
        # Constructor 
        self.real = real_part
        self.imaginary = imaginary_part

    # Addition 
    def __add__(self, other):
        real_sum = self.real + other.real
        imag_sum = self.imaginary + other.imaginary
        return ComplexNumber(real_sum, imag_sum)

    # Multiplication
    def __mul__(self, other):
        real_product = self.real * other.real - self.imaginary * other.imaginary
        imag_product = self.real * other.imaginary + self.imaginary * other.real
        return ComplexNumber(real_product, imag_product)

    # Division
    def __truediv__(self, other):
        if other.real == 0 and other.imaginary == 0:
            raise ZeroDivisionError("Cannot divide by zero!")
        divisor = other.real**2 + other.imaginary**2
        real_quotient = (self.real * other.real + self.imaginary * other.imaginary) / divisor
        imag_quotient = (self.imaginary * other.real - self.real * other.imaginary) / divisor
        return ComplexNumber(real_quotient, imag_quotient)

    # Magnitude (Absolute Value)
    def magnitude(self):
        return (self.real**2 + self.imaginary**2) ** 0.5

    # Conjugate 
    def conjugate(self):
        return ComplexNumber(self.real, -self.imaginary)

    def __str__(self):
        return f"{self.real} + {self.imaginary}i"