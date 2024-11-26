class LinearEquationSolver:
    def __init__(self, coefficient_matrix, rhs_vector):
        if coefficient_matrix.num_cols != len(rhs_vector):
            raise ValueError("Number of columns in A must match the size of b.")
        self.matrix_A = coefficient_matrix
        self.vector_b = rhs_vector    

    # Check if dimensions are compatible
    def check_dimensions(self):
        return self.matrix_A.num_rows == len(self.vector_b)

    # Check if the system is consistent
    def is_consistent(self):
        augmented_matrix = [
            self.matrix_A.grid[row_idx] + [self.vector_b[row_idx]]
            for row_idx in range(self.matrix_A.num_rows)
        ]
        rref_augmented = MatrixOperations(
            len(augmented_matrix), len(augmented_matrix[0]), augmented_matrix
        ).calculate_rref()
        
        for row in rref_augmented.grid:
            if all(val == 0 for val in row[:-1]) and row[-1] != 0:
                return False  # Inconsistent system
        return True

    # Solve the system using Gaussian elimination
    def solve_via_gaussian(self):
        if not self.is_consistent():
            raise ValueError("The system is inconsistent and cannot be solved.")
        
        augmented_matrix = [
            self.matrix_A.grid[row_idx] + [self.vector_b[row_idx]]
            for row_idx in range(self.matrix_A.num_rows)
        ]
        
        rref_solution = MatrixOperations(
            len(augmented_matrix), len(augmented_matrix[0]), augmented_matrix
        ).calculate_rref()
        
        solution_vector = [row[-1] for row in rref_solution.grid]
        return solution_vector

    # Check if one subspace is contained within another
    @staticmethod
    def is_subspace_contained(subspace1, subspace2):
        combined_space = subspace1 + subspace2
        combined_matrix = MatrixOperations(
            len(combined_space), len(combined_space[0]), combined_space
        )
        return combined_matrix.compute_rank() == len(subspace2)

    # Express solutions in terms of free variables
    def express_solution_free_vars(self):
        augmented_matrix = [
            self.matrix_A.grid[row_idx] + [self.vector_b[row_idx]]
            for row_idx in range(self.matrix_A.num_rows)
        ]
        
        rref_matrix = MatrixOperations(
            len(augmented_matrix), len(augmented_matrix[0]), augmented_matrix
        ).calculate_rref()
        
        total_vars = len(rref_matrix.grid[0]) - 1
        basic_vars = []
        free_vars = []

        for col in range(total_vars):
            if any(rref_matrix.grid[row][col] != 0 for row in range(len(rref_matrix.grid))):
                basic_vars.append(col)
            else:
                free_vars.append(col)

        general_solution = ["0"] * total_vars
        for row_idx, col_idx in enumerate(basic_vars):
            general_solution[col_idx] = str(rref_matrix.grid[row_idx][-1]) + " + " + " + ".join(
                f"{-rref_matrix.grid[row_idx][free_idx]}*x{free_idx}"
                for free_idx in free_vars if rref_matrix.grid[row_idx][free_idx] != 0
            )

        return general_solution, free_vars

    # Solve using PLU decomposition method
    def solve_with_plu(self):
        if not self.is_consistent():
            raise ValueError("The system is inconsistent and cannot be solved.")
        
        permutation_matrix, lower_matrix, upper_matrix = self.matrix_A.perform_plu_decomposition()
        
        permuted_rhs = [self.vector_b[permutation_matrix.grid[i].index(1)] for i in range(self.matrix_A.num_rows)]

        y_vector = [0] * self.matrix_A.num_rows
        for i in range(self.matrix_A.num_rows):
            y_vector[i] = permuted_rhs[i] - sum(
                lower_matrix.grid[i][j] * y_vector[j] for j in range(i)
            )

        x_vector = [0] * self.matrix_A.num_rows
        for i in range(self.matrix_A.num_rows - 1, -1, -1):
            x_vector[i] = (y_vector[i] - sum(
                upper_matrix.grid[i][j] * x_vector[j] for j in range(i + 1, self.matrix_A.num_cols)
            )) / upper_matrix.grid[i][i]

        return x_vector