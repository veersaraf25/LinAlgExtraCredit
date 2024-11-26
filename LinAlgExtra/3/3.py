class MatrixOperations:
    def __init__(self, num_rows, num_cols, matrix_entries=None, entry_type="real"):
        # Constructor 
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.entry_type = entry_type
        self.grid = matrix_entries if matrix_entries else self._get_matrix_from_user()

    def __str__(self):
        return "\n".join([" ".join(map(str, row)) for row in self.grid])

    def _get_matrix_from_user(self):
        print(f"Enter {self.num_rows * self.num_cols} values for the {self.num_rows}x{self.num_cols} matrix:")
        if self.entry_type == "real":
            return [
                [float(input(f"Value at ({i+1},{j+1}): ")) for j in range(self.num_cols)]
                for i in range(self.num_rows)
            ]
        elif self.entry_type == "complex":
            return [
                [
                    ComplexNumber(
                        float(input(f"Real part ({i+1},{j+1}): ")),
                        float(input(f"Imaginary part ({i+1},{j+1}): "))
                    )
                    for j in range(self.num_cols)
                ]
                for i in range(self.num_rows)
            ]

    @staticmethod
    def compute_vector_length(vector):
        return sum(component ** 2 for component in vector) ** 0.5

    def get_matrix_size(self):
        return self.num_rows, self.num_cols

    def compute_rank(self):
        rref_matrix = self.calculate_rref()
        rank_count = sum(1 for row in rref_matrix.grid if any(val != 0 for val in row))
        return rank_count

    def compute_nullity(self):
        return self.num_cols - self.compute_rank()

    # Calculate Reduced Row Echelon Form (RREF)
    def calculate_rref(self):
        matrix_copy = [row[:] for row in self.grid]
        pivot_row_idx = 0

        for col_idx in range(self.num_cols):
            # Find the first non-zero pivot
            pivot_idx = None
            for row_idx in range(pivot_row_idx, self.num_rows):
                if matrix_copy[row_idx][col_idx] != 0:
                    pivot_idx = row_idx
                    break

            if pivot_idx is None:
                continue

            # Swap rows to move pivot into position
            matrix_copy[pivot_row_idx], matrix_copy[pivot_idx] = (
                matrix_copy[pivot_idx],
                matrix_copy[pivot_row_idx],
            )

            pivot_element = matrix_copy[pivot_row_idx][col_idx]
            # Normalize pivot row
            matrix_copy[pivot_row_idx] = [value / pivot_element for value in matrix_copy[pivot_row_idx]]

            # Eliminate column entries below the pivot
            for row_idx in range(self.num_rows):
                if row_idx != pivot_row_idx:
                    multiplier = matrix_copy[row_idx][col_idx]
                    matrix_copy[row_idx] = [
                        current_value - multiplier * matrix_copy[pivot_row_idx][col_idx]
                        for col_idx in range(self.num_cols)
                    ]

            pivot_row_idx += 1

        return MatrixOperations(self.num_rows, self.num_cols, matrix_copy)

    # Check linear dependence of vectors
    @staticmethod
    def check_linear_dependence(vectors):
        vector_matrix = MatrixOperations(len(vectors), len(vectors[0]), vectors)
        rank_of_matrix = vector_matrix.compute_rank()
        return rank_of_matrix < len(vectors)

    # Find properties of the span of given vectors
    @staticmethod
    def find_subspace_properties(vectors):
        rref_representation = MatrixOperations(len(vectors), len(vectors[0]), vectors).calculate_rref()
        basis = [
            vector
            for vector, row in zip(vectors, rref_representation.grid)
            if any(entry != 0 for entry in row)
        ]
        dimension = len(basis)
        return dimension, basis

    # Compute rank factorization of a given matrix
    def compute_rank_factorization(self):
        rref_matrix = self.calculate_rref()
        row_space = [row for row in rref_matrix.grid if any(cell != 0 for cell in row)]
        
        column_space = [
            [self.grid[row_index][col_index] for row_index in range(self.num_rows) if self.grid[row_index][col_index] != 0]
            for col_index in range(self.num_cols)
        ]
        
        row_space_matrix = MatrixOperations(len(row_space), self.num_cols, row_space)
        column_space_matrix = MatrixOperations(self.num_rows, len(column_space), column_space)
        
        return row_space_matrix, column_space_matrix

    # Perform LU decomposition on the matrix
    def perform_lu_decomposition(self):
        if self.num_rows != self.num_cols:
            raise ValueError("LU Decomposition is applicable only to square matrices.")
        
        n = self.num_rows
        lower_matrix = [[0] * n for _ in range(n)]
        upper_matrix = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i, n):
                upper_matrix[i][j] = self.grid[i][j] - sum(
                    lower_matrix[i][k] * upper_matrix[k][j] for k in range(i)
                )

            for j in range(i, n):
                if i == j:
                    lower_matrix[i][i] = 1
                else:
                    lower_matrix[j][i] = (self.grid[j][i] - sum(
                        lower_matrix[j][k] * upper_matrix[k][i] for k in range(i)
                    )) / upper_matrix[i][i]

        return MatrixOperations(n, n, lower_matrix), MatrixOperations(n, n, upper_matrix)

    # Perform PLU decomposition on the matrix
    def perform_plu_decomposition(self):
        n = self.num_rows
        permutation_vector = list(range(n))
        
        lower_matrix = [[0] * n for _ in range(n)]
        upper_matrix = [row[:] for row in self.grid]

        for i in range(n):
            pivot_row = max(range(i, n), key=lambda r: abs(upper_matrix[r][i]))
            
            if i != pivot_row:
                upper_matrix[i], upper_matrix[pivot_row] = upper_matrix[pivot_row], upper_matrix[i]
                permutation_vector[i], permutation_vector[pivot_row] = permutation_vector[pivot_row], permutation_vector[i]

            for j in range(i + 1, n):
                factor = upper_matrix[j][i] / upper_matrix[i][i]
                lower_matrix[j][i] = factor
                
                for k in range(i, n):
                    upper_matrix[j][k] -= factor * upper_matrix[i][k]

        for i in range(n):
            lower_matrix[i][i] = 1

        permutation_matrix = [
            [1 if j == permutation_vector[i] else 0 for j in range(n)]
            for i in range(n)
        ]
        
        return (
            MatrixOperations(n, n, permutation_matrix),
            MatrixOperations(n, n, lower_matrix),
            MatrixOperations(n, n, upper_matrix),
        )