class EigenOperations:
    def __init__(self, matrix):
        if matrix.num_rows != matrix.num_cols:
            raise ValueError("Eigenvalue operations require a square matrix.")
        self.matrix = matrix

    # Compute eigenvalues using the characteristic polynomial method
    def compute_eigenvalues(self):
        n = self.matrix.num_rows
        eigenvalues = set()

        for scalar in range(-100, 101):  # Adjust range as needed
            shifted_matrix = self._shift_matrix(scalar)
            if shifted_matrix.compute_determinant() == 0:
                eigenvalues.add(scalar)

        return sorted(list(eigenvalues))

    # Compute eigenvectors corresponding to a given eigenvalue
    def compute_eigenvectors(self, eigenvalue):
        shifted_matrix = self._shift_matrix(eigenvalue)
        rref_result = shifted_matrix.calculate_rref()

        eigenvectors = []
        for i in range(self.matrix.num_cols):
            test_vector = [1 if j == i else 0 for j in range(self.matrix.num_cols)]
            result = shifted_matrix.multiply_vector(test_vector)
            if all(entry == 0 for entry in result):
                eigenvectors.append(test_vector)

        return eigenvectors

    # Check if the matrix is diagonalizable
    def is_diagonalizable(self):
        eigenvalues = self.compute_eigenvalues()
        total_eigenvector_count = 0

        for eigenvalue in eigenvalues:
            total_eigenvector_count += len(self.compute_eigenvectors(eigenvalue))

        return total_eigenvector_count == self.matrix.num_rows

    # Diagonalize the matrix if possible
    def diagonalize_matrix(self):
        if not self.is_diagonalizable():
            raise ValueError("Matrix is not diagonalizable.")

        eigenvalues = self.compute_eigenvalues()
        diagonal_matrix = [[0] * self.matrix.num_rows for _ in range(self.matrix.num_rows)]
        eigenvector_list = []

        for i, eigenvalue in enumerate(eigenvalues):
            eigenvectors = self.compute_eigenvectors(eigenvalue)
            for eigenvector in eigenvectors:
                eigenvector_list.append(eigenvector)
                diagonal_matrix[i][i] = eigenvalue

        transposed_eigenvector_matrix = MatrixOperations(
            self.matrix.num_rows, len(eigenvector_list), eigenvector_list
        )
        
        return (
            MatrixOperations(self.matrix.num_rows, self.matrix.num_cols, diagonal_matrix),
            transposed_eigenvector_matrix,
        )

    # Shift the matrix by a scalar value (A - λI)
    def _shift_matrix(self, scalar):
        identity_matrix = MatrixOperations.identity_matrix(self.matrix.num_rows)
        scaled_identity = [
            [scalar * identity_matrix.grid[i][j] for j in range(self.matrix.num_cols)]
            for i in range(self.matrix.num_rows)
        ]
        
        shifted_matrix = [
            [
                self.matrix.grid[i][j] - scaled_identity[i][j]
                for j in range(self.matrix.num_cols)
            ]
            for i in range(self.matrix.num_rows)
        ]
        
        return MatrixOperations(self.matrix.num_rows, self.matrix.num_cols, shifted_matrix)