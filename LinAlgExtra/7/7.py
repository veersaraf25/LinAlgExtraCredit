class MatrixDeterminant:
    def __init__(self, matrix):
        if matrix.num_rows != matrix.num_cols:
            raise ValueError("Determinant can only be calculated for square matrices.")
        self.matrix = matrix

    # Compute determinant using cofactor expansion
    def compute_via_cofactor(self):
        return self._determinant_recursive(self.matrix)

    # Compute determinant using PLU decomposition
    def compute_via_plu(self):
        perm_matrix, lower_matrix, upper_matrix = self.matrix.perform_plu_decomposition()
        perm_sign = self._permutation_sign(perm_matrix)
        upper_det = self._upper_triangular_determinant(upper_matrix)
        return perm_sign * upper_det

    # Compute determinant using RREF
    def compute_via_rref(self):
        rref_matrix = self.matrix.calculate_rref()
        diagonal_product = 1
        for i in range(self.matrix.num_rows):
            diagonal_product *= rref_matrix.grid[i][i]
        return diagonal_product

    # Recursive method to calculate determinant
    def _determinant_recursive(self, matrix):
        if matrix.num_rows == 1:
            return matrix.grid[0][0]

        det_value = 0
        for col_idx in range(matrix.num_cols):
            minor_matrix = self._extract_minor(matrix, 0, col_idx)
            cofactor = ((-1) ** col_idx) * matrix.grid[0][col_idx] * self._determinant_recursive(minor_matrix)
            det_value += cofactor
        return det_value

    # Extract minor matrix by excluding a specific row and column
    def _extract_minor(self, matrix, exclude_row, exclude_col):
        minor = [
            [
                matrix.grid[row_idx][col_idx]
                for col_idx in range(matrix.num_cols) if col_idx != exclude_col
            ]
            for row_idx in range(matrix.num_rows) if row_idx != exclude_row
        ]
        return MatrixOperations(matrix.num_rows - 1, matrix.num_cols - 1, minor)

    # Calculate the sign of the permutation from the permutation matrix
    def _permutation_sign(self, perm_matrix):
        perm_vector = [perm_matrix.grid[i].index(1) for i in range(self.matrix.num_rows)]
        sign = 1
        for i in range(len(perm_vector)):
            for j in range(i + 1, len(perm_vector)):
                if perm_vector[i] > perm_vector[j]:
                    sign *= -1
        return sign

    # Calculate determinant of an upper triangular matrix
    def _upper_triangular_determinant(self, upper_matrix):
        determinant = 1
        for i in range(self.matrix.num_rows):
            determinant *= upper_matrix.grid[i][i]
        return determinant