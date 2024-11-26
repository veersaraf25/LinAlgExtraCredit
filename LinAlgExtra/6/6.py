class BasisOperations:
    def __init__(self, basis_vectors):
        self.basis = basis_vectors
        self.dimension = len(basis_vectors)

    # Check if a vector is within the span of a set of vectors
    @staticmethod
    def is_within_span(target_vector, spanning_vectors):
        augmented_set = spanning_vectors + [target_vector]
        span_matrix = MatrixOperations(len(augmented_set), len(target_vector), augmented_set)
        return span_matrix.compute_rank() == len(spanning_vectors)

    # Represent a target vector as a linear combination of basis vectors
    @staticmethod
    def represent_as_combination(target_vector, basis_set):
        augmented = basis_set + [target_vector]
        augmented_matrix = MatrixOperations(len(augmented), len(target_vector), augmented)
        rref_result = augmented_matrix.calculate_rref()

        coefficients = []
        for row in rref_result.grid:
            if any(row[:-1]):
                coefficients.append(row[-1])

        return coefficients if len(coefficients) == len(basis_set) else None

    # Check if two sets of vectors span the same subspace
    @staticmethod
    def span_equality_check(set_a, set_b):
        combined_set = set_a + set_b
        combined_rank = MatrixOperations(len(combined_set), len(combined_set[0]), combined_set).compute_rank()
        return combined_rank == len(set_a) == len(set_b)

    # Calculate coordinates of a target vector in terms of an ordered basis
    @staticmethod
    def calculate_coordinates(target_vector, ordered_basis):
        representation = BasisOperations.represent_as_combination(target_vector, ordered_basis)
        if representation is None:
            raise ValueError("The vector does not belong to the span of the given basis.")
        return representation

    # Derive the change of basis matrix from one basis to another
    @staticmethod
    def derive_basis_change_matrix(starting_basis, target_basis):
        starting_inverse = MatrixInversion(MatrixOperations(len(starting_basis), len(starting_basis), starting_basis)).inverse_via_row_reduction()
        target_matrix = MatrixOperations(len(target_basis), len(target_basis[0]), target_basis)
        
        basis_change_matrix = starting_inverse.multiply(target_matrix)
        return basis_change_matrix

    # Apply a change of basis transformation to vector coordinates
    @staticmethod
    def apply_basis_change(vector_coordinates, transformation_matrix):
        new_coordinates = [
            sum(vector_coordinates[j] * transformation_matrix.grid[j][i] for j in range(len(vector_coordinates)))
            for i in range(len(transformation_matrix.grid[0]))
        ]
        return new_coordinates