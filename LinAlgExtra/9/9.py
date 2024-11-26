class QRDecomposition:
    def __init__(self, matrix):
        if matrix.num_cols > matrix.num_rows:
            raise ValueError("Matrix must have rows >= columns for QR decomposition.")
        self.matrix = matrix

    # Perform Gram-Schmidt orthogonalization
    def gram_schmidt(self):
        vectors = [list(col) for col in zip(*self.matrix.grid)]  # Columns as vectors
        orthogonal_set = []

        for vector in vectors:
            projection = [0] * len(vector)
            for basis_vector in orthogonal_set:
                proj_coefficient = self._dot_product(vector, basis_vector) / self._dot_product(basis_vector, basis_vector)
                projection = [proj + proj_coefficient * basis_comp for proj, basis_comp in zip(projection, basis_vector)]
            orthogonal_vector = [vec - proj for vec, proj in zip(vector, projection)]
            orthogonal_set.append(orthogonal_vector)

        orthogonal_matrix = list(zip(*orthogonal_set)) 
        return MatrixOperations(len(orthogonal_matrix), len(orthogonal_matrix[0]), orthogonal_matrix)

    # Normalize the orthogonal basis to produce an orthonormal basis
    def normalize_basis(self, orthogonal_matrix):
        normalized_vectors = []

        for vector in zip(*orthogonal_matrix.grid):  # Treat columns as vectors
            norm = self._vector_norm(vector)
            normalized_vector = [component / norm for component in vector]
            normalized_vectors.append(normalized_vector)

        normalized_matrix = list(zip(*normalized_vectors))  # Convert back to matrix form
        return MatrixOperations(len(normalized_matrix), len(normalized_matrix[0]), normalized_matrix)

    # Perform QR decomposition
    def qr_decomposition(self):
        orthogonal_matrix = self.gram_schmidt()
        q_matrix = self.normalize_basis(orthogonal_matrix)

        q_transpose = q_matrix.get_transpose()
        r_matrix = MatrixOperations.multiply_matrices(q_transpose, self.matrix)

        return q_matrix, r_matrix

    @staticmethod
    def _dot_product(vector1, vector2):
        return sum(v1 * v2 for v1, v2 in zip(vector1, vector2))

    @staticmethod
    def _vector_norm(vector):
        return sum(component ** 2 for component in vector) ** 0.5