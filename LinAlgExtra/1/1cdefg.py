class Matrix:
    def __init__(self, rows, cols, entry_type="real"):
        # Constructor
        self.row_count = rows
        self.column_count = cols
        self.entry_type = entry_type
        
        if entry_type == "real":
            self.values = [
                [float(input(f"Enter value for entry ({i+1},{j+1}): ")) for j in range(cols)]
                for i in range(rows)
            ]
        elif entry_type == "complex":
            self.values = [
                [
                    ComplexNumber(
                        float(input(f"Real part ({i+1},{j+1}): ")),
                        float(input(f"Imaginary part ({i+1},{j+1}): "))
                    )
                    for j in range(cols)
                ]
                for i in range(rows)
            ]

    # Matrix Addition
    def add(self, another_matrix):
        if self.row_count != another_matrix.row_count or self.column_count != another_matrix.column_count:
            raise ValueError("Matrix dimensions must match for addition.")
        
        summed_values = [
            [self.values[i][j] + another_matrix.values[i][j] for j in range(self.column_count)]
            for i in range(self.row_count)
        ]
        
        return Matrix(self.row_count, self.column_count, self.entry_type)

    # Matrix Multiplication
    def multiply(self, another_matrix):
        if self.column_count != another_matrix.row_count:
            raise ValueError("Matrix multiplication requires compatible dimensions.")
        
        product_values = [[ComplexNumber(0, 0) for _ in range(another_matrix.column_count)] for _ in range(self.row_count)]
        
        for i in range(self.row_count):
            for j in range(another_matrix.column_count):
                product_values[i][j] = sum(
                    self.values[i][k] * another_matrix.values[k][j] for k in range(self.column_count)
                )
        
        return Matrix(self.row_count, another_matrix.column_count, self.entry_type)

    def get_row(self, index):
        return self.values[index]

    def get_column(self, index):
        return [self.values[i][index] for i in range(self.row_count)]

    def transpose(self):
        transposed_values = [[self.values[j][i] for j in range(self.row_count)] for i in range(self.column_count)]
        
        return Matrix(self.column_count, self.row_count, self.entry_type)

    def conjugate(self):
        if self.entry_type != "complex":
            raise ValueError("Conjugate is only valid for complex matrices.")
        
        conjugated_values = [[entry.conjugate() for entry in row] for row in self.values]
        
        return Matrix(self.row_count, self.column_count, self.entry_type)

    def hermitian(self):
        return self.transpose().conjugate()

    def __str__(self):
        return "\n".join([" ".join(map(str, row)) for row in self.values])

class MatrixFromVectors(Matrix):
    def __init__(self, vectors):
        if not vectors:
            raise ValueError("At least one vector is required to create a matrix.")

        col_count = len(vectors)
        row_count = len(vectors[0].elements)
        
        for vector in vectors:
            if len(vector.elements) != row_count:
                raise ValueError("All vectors must have the same length.")

        super().__init__(row_count, col_count, vectors[0].field_type)
        
        # Populate values from vectors
        self.values = [[vector.elements[row] for vector in vectors] for row in range(row_count)]

    def __str__(self):
        return "\n".join([" ".join(map(str, row)) for row in self.values])