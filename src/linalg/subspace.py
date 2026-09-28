import sympy as sp

from .vectors import vector, zero_vector, is_zero


class Subspace:
    def __init__(self, generators, ambient_dimension=None):
        """
        generators:
            Vectors whose span defines the subspace.

        ambient_dimension:
            Dimension of the surrounding vector space F^n.
            Usually inferred from the generators.
        """

        self.generators = [vector(v) for v in generators]

        if not self.generators:
            if ambient_dimension is None:
                raise ValueError(
                    "ambient_dimension is required for the zero subspace."
                )

            self.ambient_dimension = ambient_dimension
            return

        self.ambient_dimension = len(self.generators[0])

        for v in self.generators:
            if len(v) != self.ambient_dimension:
                raise ValueError(
                    "All generators must have the same dimension."
                )

    # ---------------------------------------------------------
    # Basic helpers
    # ---------------------------------------------------------

    def matrix(self):
        """Matrix whose columns are the generators."""

        if not self.generators:
            return sp.zeros(self.ambient_dimension, 0)

        return sp.Matrix.hstack(*self.generators)

    def basis(self):
        """
        Return an independent set of vectors spanning the subspace.
        """

        return self.matrix().columnspace()

    @classmethod
    def zero(cls, ambient_dimension):
        return cls([], ambient_dimension=ambient_dimension)

    def is_zero_subspace(self):
        return len(self.basis()) == 0

    # ---------------------------------------------------------
    # Membership
    # ---------------------------------------------------------

    def contains(self, v):
        """
        Check whether v belongs to this subspace.
        """

        v = vector(v)

        if len(v) != self.ambient_dimension:
            return False

        if self.is_zero_subspace():
            return is_zero(v)

        A = self.matrix()

        return A.rank() == A.row_join(v).rank()

    def __contains__(self, v):
        return self.contains(v)

    # ---------------------------------------------------------
    # Sum of arbitrary subspaces
    # ---------------------------------------------------------

    def sum(self, *others):
        """
        Return

            self + U2 + ... + Um.

        The sum is spanned by all generators from all subspaces.
        """

        subspaces = (self,) + others

        self._check_same_ambient_spaces(*others)

        generators = []

        for subspace in subspaces:
            generators.extend(subspace.generators)

        return Subspace(
            generators,
            ambient_dimension=self.ambient_dimension,
        )

    def __add__(self, other):
        """
        Allows:

            U + W
        """

        return self.sum(other)

    # ---------------------------------------------------------
    # Intersection of arbitrary subspaces
    # ---------------------------------------------------------

    def intersection(self, *others):
        """
        Return

            self ∩ U2 ∩ ... ∩ Um.
        """

        self._check_same_ambient_spaces(*others)

        result = self

        for other in others:
            result = result._intersection_two(other)

        return result

    def _intersection_two(self, other):
        """
        Compute the intersection of exactly two subspaces.

        If

            x = A a
            x = B b

        then

            A a - B b = 0.

        Thus we find the null space of

            [A  -B].
        """

        if self.is_zero_subspace() or other.is_zero_subspace():
            return Subspace.zero(self.ambient_dimension)

        A_basis = self.basis()
        B_basis = other.basis()

        A = sp.Matrix.hstack(*A_basis)
        B = sp.Matrix.hstack(*B_basis)

        combined = A.row_join(-B)

        solutions = combined.nullspace()

        intersection_vectors = []

        for solution in solutions:
            a = solution[:A.cols, :]

            common_vector = A * a

            if not is_zero(common_vector):
                intersection_vectors.append(common_vector)

        if not intersection_vectors:
            return Subspace.zero(self.ambient_dimension)

        basis = sp.Matrix.hstack(
            *intersection_vectors
        ).columnspace()

        return Subspace(
            basis,
            ambient_dimension=self.ambient_dimension,
        )

    # ---------------------------------------------------------
    # Direct sum of arbitrary subspaces
    # ---------------------------------------------------------

    def is_direct_sum(self, *others):
        """
        Check whether

            self + U2 + ... + Um

        is a direct sum.

        The sum is direct iff

            u1 + u2 + ... + um = 0

        implies

            u1 = u2 = ... = um = 0.
        """

        subspaces = (self,) + others

        self._check_same_ambient_spaces(*others)

        basis_vectors = []

        for subspace in subspaces:
            basis_vectors.extend(subspace.basis())

        if not basis_vectors:
            return True

        A = sp.Matrix.hstack(*basis_vectors)

        return A.rank() == A.cols

    # ---------------------------------------------------------
    # Decomposition across arbitrary subspaces
    # ---------------------------------------------------------

    def decompose(self, v, *others):
        """
        Given

            V = self ⊕ U2 ⊕ ... ⊕ Um,

        find the unique vectors

            u1, u2, ..., um

        such that

            v = u1 + u2 + ... + um.
        """

        subspaces = (self,) + others

        self._check_same_ambient_spaces(*others)

        v = vector(v)

        if len(v) != self.ambient_dimension:
            raise ValueError(
                "Vector has the wrong ambient dimension."
            )

        if not self.is_direct_sum(*others):
            raise ValueError(
                "The subspaces do not form a direct sum."
            )

        total_space = self.sum(*others)

        if not total_space.contains(v):
            raise ValueError(
                "v does not belong to the sum of the subspaces."
            )

        bases = [
            subspace.basis()
            for subspace in subspaces
        ]

        all_basis_vectors = [
            b
            for basis in bases
            for b in basis
        ]

        if not all_basis_vectors:
            if is_zero(v):
                return tuple(
                    zero_vector(self.ambient_dimension)
                    for _ in subspaces
                )

            raise ValueError(
                "Nonzero vector cannot belong to the zero subspace."
            )

        A = sp.Matrix.hstack(*all_basis_vectors)

        solution, _ = A.gauss_jordan_solve(v)

        components = []

        start = 0

        for basis in bases:
            count = len(basis)

            if count == 0:
                component = zero_vector(
                    self.ambient_dimension
                )
            else:
                B = sp.Matrix.hstack(*basis)

                coefficients = solution[
                    start:start + count,
                    :
                ]

                component = B * coefficients

            components.append(component)

            start += count

        return tuple(components)

    # ---------------------------------------------------------
    # Validation
    # ---------------------------------------------------------

    def _check_same_ambient_spaces(self, *others):
        for other in others:
            if (
                self.ambient_dimension
                != other.ambient_dimension
            ):
                raise ValueError(
                    "All subspaces must belong to "
                    "the same ambient vector space."
                )