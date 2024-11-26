class MatrixProperties:
    def __init__(self, num_rows, num_cols, elements=None, matrix_type="real"):
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.matrix_type = matrix_type
        if elements:
            self.data = elements
        else:
            self.data = []
            print(f"Enter {num_rows * num_cols} values for the {num_rows}x{num_cols} matrix:")
            if matrix_type == "real":
                self.data = [
                    [float(input(f"Element at ({i+1},{j+1}): ")) for j in range(num_cols)]
                    for i in range(num_rows)
                ]
            elif matrix_type == "complex":
                self.data = [
                    [
                        ComplexNumber(
                            float(input(f"Real part at ({i+1},{j+1}): ")),
                            float(input(f"Imaginary part at ({i+1},{j+1}): "))
                        )
                        for j in range(num_cols)
                    ]
                    for i in range(num_rows)
                ]

    def __str__(self):
        return "\n".join([" ".join(map(str, row)) for row in self.data])

    def is_zero_matrix(self):
        return all(all(cell == 0 for cell in row) for row in self.data)

    def is_symmetric(self):
        if not self.is_square():
            return False
        return all(
            self.data[i][j] == self.data[j][i]
            for i in range(self.num_rows) for j in range(self.num_cols)
        )

    def is_hermitian(self):
        if self.matrix_type != "complex" or not self.is_square():
            return False
        return all(
            self.data[i][j] == self.data[j][i].conjugate()
            for i in range(self.num_rows) for j in range(self.num_cols)
        )

    def is_square(self):
        return self.num_rows == self.num_cols

    def is_orthogonal(self):
        if not self.is_square():
            return False
        transpose_matrix = self.transpose()
        product_matrix = self.multiply(transpose_matrix)
        identity_matrix = MatrixProperties.identity_matrix(self.num_rows)
        return product_matrix.data == identity_matrix.data

    def is_unitary(self):
        if not self.is_square():
            return False
        hermitian_matrix = self.hermitian()
        product_matrix = self.multiply(hermitian_matrix)
        identity_matrix = MatrixProperties.identity_matrix(self.num_rows)
        return product_matrix.data == identity_matrix.data

    def is_scalar(self):
        if not self.is_square():
            return False
        diag_value = self.data[0][0]
        return all(
            (self.data[i][j] == diag_value if i == j else self.data[i][j] == 0)
            for i in range(self.num_rows) for j in range(self.num_cols)
        )

    def is_singular(self):
        return self.determinant() == 0

    def is_invertible(self):
        return not self.is_singular()

    def is_identity(self):
        if not self.is_square():
            return False
        return all(
            (self.data[i][j] == 1 if i == j else self.data[i][j] == 0)
            for i in range(self.num_rows) for j in range(self.num_cols)
        )

    def is_nilpotent(self):
        if not self.is_square():
            return False
        power_matrix = self
        for _ in range(1, self.num_rows + 1):
            power_matrix = power_matrix.multiply(self)
            if power_matrix.is_zero_matrix():
                return True
        return False

    def is_diagonalizable(self):
        eigenvalues = self.get_eigenvalues()
        eigenspaces = [self.get_eigenspace(value) for value in eigenvalues]
        total_space_dimension = sum(len(space) for space in eigenspaces)
        return total_space_dimension == self.num_rows

    def has_lu_decomposition(self):
        if not self.is_square():
            return False
        return not self.is_singular()

    def transpose(self):
        transposed_data = [[self.data[j][i] for j in range(self.num_rows)] for i in range(self.num_cols)]
        return MatrixProperties(
            num_rows=self.num_cols,
            num_cols=self.num_rows,
            elements=transposed_data,
            matrix_type=self.matrix_type
        )

    def multiply(self, other):
        if self.num_cols != other.num_rows:
            raise ValueError("Invalid dimensions for multiplication.")
        
        result_data = [
            [
                sum(
                    (self.data[i][k] * other.data[k][j] if isinstance(self.data[i][k], (int, float)) 
                     else (self.data[i][k] * other.data[k][j]) ) 
                    for k in range(self.num_cols)
                )
                for j in range(other.num_cols)
            ]
            for i in range(self.num_rows)
        ]
        
        return MatrixProperties(num_rows=self.num_rows, num_cols=other.num_cols, elements=result_data, matrix_type=self.matrix_type)

    @staticmethod
    def identity_matrix(size):
        identity_elements = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
        return MatrixProperties(size, size, elements=identity_elements)

    # Method to calculate the determinant using cofactor expansion
    def determinant(self):
        if not self.is_square():
            raise ValueError("Determinant can only be calculated for square matrices.")
        
        if self.num_rows == 1:
            return self.data[0][0]
        
        if self.num_rows == 2:
            # For a 2x2 matrix: |a b| -> ad - bc |c d|
            a, b = self.data[0]
            c, d = self.data[1]
            return a * d - b * c
        
        det_value = 0
        for c in range(self.num_cols):
            # Minor matrix after removing first row and current column
            minor = [row[:c] + row[c+1:] for row in (self.data[1:])]
            det_value += ((-1) ** c) * self.data[0][c] * MatrixProperties(len(minor), len(minor[0]), elements=minor).determinant()
        
        return det_value

    # Method to get eigenvalues using characteristic polynomial (simplified approach)
    def get_eigenvalues(self):
        # This is a placeholder; actual implementation would require solving the characteristic polynomial.
        
       # For demonstration purposes, returning a fixed set of eigenvalues.
       # Replace this with an actual computation based on the characteristic polynomial.
       eigenvalues_placeholder = [1, -1]  
       return eigenvalues_placeholder 

    # Method to get eigenspace corresponding to an eigenvalue (placeholder implementation)
    def get_eigenspace(self, eigenvalue):
       # This is a placeholder; actual implementation would require solving (A - λI)x = 0.
       
       # For demonstration purposes, returning a fixed eigenspace.
       # Replace this with an actual computation based on the eigenvalue.
       eigenspace_placeholder = []  
       return eigenspace_placeholder 