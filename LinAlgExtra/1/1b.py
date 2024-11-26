class Vector:
    def __init__(self, dimension, field_type="real"):
        # Constructor
        self.length = dimension
        self.field_type = field_type
        self.elements = []
        
        if field_type == "real":
            # Input Real Numbers
            self.elements = [float(input(f"Enter component {i+1}: ")) for i in range(dimension)]
        elif field_type == "complex":
            # Input Complex Numbers
            self.elements = [
                ComplexNumber(
                    float(input(f"Real part of component {i+1}: ")),
                    float(input(f"Imaginary part of component {i+1}: "))
                )
                for i in range(dimension)
            ]

    # Vector Addition
    def add(self, another_vector):
        if len(self.elements) != len(another_vector.elements):
            raise ValueError("Vectors must be of the same length for addition.")
        
        summed_elements = [
            self.elements[i] + another_vector.elements[i] for i in range(self.length)
        ]
        
        return Vector(self.length, self.field_type)

    # Scalar Multiplication
    def scale(self, scalar):
        if self.field_type == "real":
            scaled_elements = [comp * scalar for comp in self.elements]
        elif self.field_type == "complex":
            scaled_elements = [comp * ComplexNumber(scalar, 0) for comp in self.elements]
        
        return Vector(self.length, self.field_type)

    def __str__(self):
        return f"Vector({self.elements})"