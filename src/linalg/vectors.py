import sympy as sp

Vector = sp.Matrix 

def vector(values) -> Vector: 
    """ Create a column vector """
    return sp.Matrix(values)


def zero_vector(dimension:int) -> Vector: 
    """Return the zero vector in F^n."""
    return sp.zeros(dimension, 1)

def add(v: Vector, w: Vector) -> Vector:
    """Return v + w."""
    _check_same_dimension(v, w)
    return v + w


def subtract(v: Vector, w: Vector) -> Vector:
    """Return v - w."""
    _check_same_dimension(v, w)
    return v - w


def scalar_multiply(a, v: Vector) -> Vector:
    """Return a*v."""
    return a * v


def linear_combination(coefficients, vectors) -> Vector:
    """
    Compute

        a1*v1 + ... + an*vn
    """
    if len(coefficients) != len(vectors):
        raise ValueError(
            "Number of coefficients must equal number of vectors."
        )

    if not vectors:
        raise ValueError("At least one vector is required.")

    dimension = len(vectors[0])

    for v in vectors:
        if len(v) != dimension:
            raise ValueError("All vectors must have the same dimension.")

    result = zero_vector(dimension)

    for a,v in zip(coefficients, vectors):
        result += a * v 

    return result 

def is_zero(v: Vector) -> bool:
    """Check whether v is the zero vector."""
    return v == zero_vector(len(v))


def _check_same_dimension(v: Vector, w: Vector):
    if len(v) != len(w):
        raise ValueError("Vectors must have the same dimension.")
