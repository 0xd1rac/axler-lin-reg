# axler-lin-reg

An educational linear algebra library built while working through Sheldon Axler's *Linear Algebra Done Right*.

The goal of this project is **not** to replace NumPy, SciPy, or SymPy.

The goal is to make abstract linear algebra concepts:

- executable
- inspectable
- easier to reason about
- closely connected to their mathematical definitions

The library is built chapter-by-chapter alongside *Linear Algebra Done Right*.

## Why

Textbook definitions like "the sum $U_1 + \cdots + U_m$ is direct if the only way to write $0$ as $u_1 + \cdots + u_m$ is by taking every $u_k = 0$" are precise, but they're hard to *poke at*. This library turns those definitions into code you can run on concrete examples, so you can test your intuition, check exercises, and see where a definition does real work.

Some design choices that follow from this:

- **Exact arithmetic.** Vectors are SymPy matrices, so `1/3` stays `1/3` and rank computations don't suffer from floating-point error.
- **Code follows the definitions.** Docstrings state the mathematical fact being implemented, and function names match the book's vocabulary (`linear_combination`, `is_direct_sum`, `decompose`).
- **Readability over speed.** When there's a choice between a clever algorithm and one that mirrors the proof, the proof wins.

## Installation

Requires Python 3.14+ and [SymPy](https://www.sympy.org/).

```bash
git clone <repo-url> axler-lin-reg
cd axler-lin-reg
uv sync          # or: pip install -e .
```

This installs the `linalg` package (from `src/linalg`) along with SymPy. Then:

```bash
uv run python
```

## Quick tour

### Vectors in $\mathbf{F}^n$

```python
from linalg.vectors import vector, add, scalar_multiply, linear_combination, zero_vector, is_zero

v = vector([1, 2, 3])
w = vector([4, 5, 6])

add(v, w)                                # (5, 7, 9)
scalar_multiply(2, v)                    # (2, 4, 6)
linear_combination([1, -1], [v, w])      # 1·v + (-1)·w = (-3, -3, -3)
is_zero(zero_vector(3))                  # True
```

### Subspaces, sums, and direct sums (§1C)

A `Subspace` is defined by a list of vectors whose span is the subspace.

```python
from linalg.subspace import Subspace

U = Subspace([[1, 0, 0]])     # the x-axis in F^3
W = Subspace([[0, 1, 0]])     # the y-axis
X = Subspace([[1, 1, 0]])     # the line y = x

# Membership
[3, 4, 0] in U + W            # True
[0, 0, 1] in U + W            # False

# Direct sums
U.is_direct_sum(W)            # True:  U ∩ W = {0}
U.is_direct_sum(W, X)         # False: (1,0,0) + (0,1,0) - (1,1,0) = 0

# Unique decomposition in a direct sum
U.decompose([3, 4, 0], W)     # ((3, 0, 0), (0, 4, 0))
```

Intersections:

```python
P = Subspace([[1, 0, 0], [0, 1, 0]])   # the xy-plane
Q = Subspace([[0, 1, 0], [0, 0, 1]])   # the yz-plane

P.intersection(Q).basis()             # [(0, 1, 0)]  — the y-axis
```

## What's implemented

| Book section | Concept | Where |
| --- | --- | --- |
| 1A | Vectors in $\mathbf{F}^n$, addition, scalar multiplication | `linalg/vectors.py` |
| 1A | Linear combinations, the zero vector | `linalg/vectors.py` |
| 1C | Subspaces (as spans), membership, zero subspace | `linalg/subspace.py` |
| 1C | Sums of subspaces $U_1 + \cdots + U_m$ | `Subspace.sum`, `U + W` |
| 1C | Intersections $U_1 \cap \cdots \cap U_m$ | `Subspace.intersection` |
| 1C | Direct sums and the unique-decomposition property | `Subspace.is_direct_sum`, `Subspace.decompose` |

## Roadmap

Following the chapters of *Linear Algebra Done Right* (4th edition):

- [x] **Ch. 1 — Vector Spaces:** $\mathbf{F}^n$, subspaces, sums and direct sums
- [ ] **Ch. 2 — Finite-Dimensional Vector Spaces:** span, linear independence, bases, dimension
- [ ] **Ch. 3 — Linear Maps:** null spaces, ranges, matrices, invertibility, quotient spaces, duality
- [ ] **Ch. 4 — Polynomials**
- [ ] **Ch. 5 — Eigenvalues and Eigenvectors**
- [ ] **Ch. 6 — Inner Product Spaces**
- [ ] **Ch. 7 — Operators on Inner Product Spaces**
- [ ] **Ch. 8 — Operators on Complex Vector Spaces**
- [ ] **Ch. 9 — Multilinear Algebra and Determinants**

## Project layout

```
src/linalg/
├── vectors.py     # vectors in F^n and basic operations
└── subspace.py    # subspaces, sums, intersections, direct sums
```

## Acknowledgements

Built as a study companion to Sheldon Axler, *Linear Algebra Done Right* (Springer), which is [freely available online](https://linear.axler.net/).
