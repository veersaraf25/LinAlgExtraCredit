class SingularValueDecomposition:
    def __init__(self, matrix):
        # Constructor for SVD operations
        self.matrix = matrix

    # Compute singular values using the characteristic polynomial method
    def compute_singular_values(self):
        transpose_matrix = self.matrix.get_transpose()
        ata_matrix = MatrixOperations.multiply_matrices(transpose_matrix, self.matrix)

        eigen_operations = EigenOperations(ata_matrix)
        eigenvalues = eigen_operations.compute_eigenvalues()

        singular_values = [abs(eigenval) ** 0.5 for eigenval in eigenvalues if eigenval >= 0]
        return sorted(singular_values, reverse=True)

    # Compute U matrix from singular values
    def compute_u_matrix(self, singular_values):
        m, n = self.matrix.num_rows, self.matrix.num_cols
        u_matrix = []

        for singular_value in singular_values:
            if singular_value == 0:
                continue
            left_vector = self._compute_left_singular_vector(singular_value)
            u_matrix.append(left_vector)

        # Ensure U is m x m by appending orthogonal vectors if needed
        while len(u_matrix) < m:
            u_matrix.append([0] * m)

        return MatrixOperations(m, m, u_matrix)

    # Compute V matrix from the matrix
    def compute_v_matrix(self):
        transpose_matrix = self.matrix.get_transpose()
        eigen_operations = EigenOperations(MatrixOperations.multiply_matrices(transpose_matrix, self.matrix))
        eigenvalues = eigen_operations.compute_eigenvalues()

        v_matrix = []
        for eigenvalue in eigenvalues:
            eigenvectors = eigen_operations.compute_eigenvectors(eigenvalue)
            v_matrix.extend(eigenvectors)

        return MatrixOperations(self.matrix.num_cols, self.matrix.num_cols, v_matrix)

    # Construct the diagonal S matrix from singular values
    def compute_s_matrix(self, singular_values):
        s_matrix = [[0] * self.matrix.num_cols for _ in range(self.matrix.num_rows)]

        for i in range(min(self.matrix.num_rows, self.matrix.num_cols)):
            if i < len(singular_values):
                s_matrix[i][i] = singular_values[i]

        return MatrixOperations(self.matrix.num_rows, self.matrix.num_cols, s_matrix)

    # Perform SVD and return U, S, and V matrices
    def perform_svd(self):
        singular_values = self.compute_singular_values()
        u_matrix = self.compute_u_matrix(singular_values)
        v_matrix = self.compute_v_matrix()
        s_matrix = self.compute_s_matrix(singular_values)
        
        return u_matrix, s_matrix, v_matrix

    # Compute left singular vector corresponding to a given singular value
    def _compute_left_singular_vector(self, singular_value):
        transpose_matrix = self.matrix.get_transpose()
        ata_matrix = MatrixOperations.multiply_matrices(transpose_matrix, self.matrix)
        
        eigen_operations = EigenOperations(ata_matrix)
        eigenvectors = eigen_operations.compute_eigenvectors(singular_value ** 2)

        left_vector = eigenvectors[0]
        norm = sum(x ** 2 for x in left_vector) ** 0.5
        return [x / norm for x in left_vector]