class MatrixInversion:
    def __init__(self, square_matrix):
        if square_matrix.num_rows != square_matrix.num_cols:
            raise ValueError("Matrix must be square for inversion operations.")
        self.square_matrix = square_matrix

    # Check if the matrix is invertible
    def check_invertibility(self):
        det_value = self.calculate_determinant()
        return det_value != 0

    # Inverse using row reduction method
    def inverse_via_row_reduction(self):
        if not self.check_invertibility():
            raise ValueError("Matrix is not invertible.")
        
        n = self.square_matrix.num_rows
        augmented_matrix = [
            self.square_matrix.grid[i] + ([1 if i == j else 0 for j in range(n)])
            for i in range(n)
        ]
        
        rref = MatrixOperations(n, 2 * n, augmented_matrix).calculate_rref()

        inverse_matrix = [
            row[n:] for row in rref.grid
        ]
        return MatrixOperations(n, n, inverse_matrix)

    # Inverse using adjoint method
    def inverse_via_adjoint(self):
        if not self.check_invertibility():
            raise ValueError("Matrix is not invertible.")
        
        det_value = self.calculate_determinant()
        cofactor_matrix = self.generate_cofactor_matrix()
        adjoint_matrix = cofactor_matrix.get_transpose()

        inverse_matrix = [
            [entry / det_value for entry in row]
            for row in adjoint_matrix.grid
        ]
        return MatrixOperations(self.square_matrix.num_rows, self.square_matrix.num_cols, inverse_matrix)

    # Calculate determinant using cofactor expansion
    def calculate_determinant(self):
        if self.square_matrix.num_rows == 1:
            return self.square_matrix.grid[0][0]
        
        determinant_value = 0
        for col_idx in range(self.square_matrix.num_cols):
            minor_matrix = self.extract_minor(0, col_idx)
            cofactor = ((-1) ** col_idx) * self.square_matrix.grid[0][col_idx] * minor_matrix.calculate_determinant()
            determinant_value += cofactor
        return determinant_value

    # Generate the cofactor matrix
    def generate_cofactor_matrix(self):
        cofactor_grid = []
        for row_idx in range(self.square_matrix.num_rows):
            cofactor_row = []
            for col_idx in range(self.square_matrix.num_cols):
                minor_matrix = self.extract_minor(row_idx, col_idx)
                cofactor_value = ((-1) ** (row_idx + col_idx)) * minor_matrix.calculate_determinant()
                cofactor_row.append(cofactor_value)
            cofactor_grid.append(cofactor_row)
        return MatrixOperations(self.square_matrix.num_rows, self.square_matrix.num_cols, cofactor_grid)

    # Extract the minor matrix by excluding a specific row and column
    def extract_minor(self, exclude_row, exclude_col):
        minor = [
            [
                self.square_matrix.grid[row_idx][col_idx]
                for col_idx in range(self.square_matrix.num_cols) if col_idx != exclude_col
            ]
            for row_idx in range(self.square_matrix.num_rows) if row_idx != exclude_row
        ]
        return MatrixOperations(self.square_matrix.num_rows - 1, self.square_matrix.num_cols - 1, minor)