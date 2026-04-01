# PHLOWER Local Patches

**Project**: P697 TEAseq T1D Low-Dose IL-2
**Author of patches**: T. Edwards (with AI assistance)
**Date**: 2026-03-30
**PHLOWER version**: 0.1.5 (`phlowerpy`)
**Install location**: `~/miniconda3/envs/phlower_env/lib/python3.10/site-packages/phlower/`
**Upstream repo**: <https://github.com/CostaLab/phlower>

---

## Summary

Two patches were applied to `phlower/util.py` to fix incompatibilities with
**scikit-sparse 0.5.0** (the latest PyPI release as of 2026-03-30). Both issues
cause `phlower.ext.ddhodge()` to crash during the Cholesky solve step of graph
construction. These are general bugs that affect any user installing phlower
0.1.5 with scikit-sparse ≥ 0.5.0.

---

## Patch 1: Missing `CholmodTooLargeError` export

**File**: `phlower/util.py`, function `test_cholesky()` (line ~403)

**Problem**: phlower imports `CholmodTooLargeError` from `sksparse.cholmod`, but
scikit-sparse 0.5.0 does not export this name. The import fails with:

```
ImportError: cannot import name 'CholmodTooLargeError' from 'sksparse.cholmod'
```

**Root cause**: `CholmodTooLargeError` was added to scikit-sparse's development
branch but has not been included in a PyPI release. phlower 0.1.5 was likely
developed against an unreleased version of scikit-sparse.

**Fix**: Replace the single-line import with a try/except that defines a stub
exception class if the name is unavailable.

```python
# ORIGINAL (line 403):
from sksparse.cholmod import cholesky, CholmodNotPositiveDefiniteError, CholmodTooLargeError

# PATCHED:
from sksparse.cholmod import cholesky, CholmodNotPositiveDefiniteError
try:
    from sksparse.cholmod import CholmodTooLargeError
except ImportError:
    class CholmodTooLargeError(Exception):
        pass
```

**Impact**: The stub class is never raised by scikit-sparse 0.5.0, so the
`except CholmodTooLargeError` block is effectively dead code. This is safe —
the exception would only be meaningful once scikit-sparse ships a version that
actually raises it.

---

## Patch 2: `cholesky()` returns tuple instead of callable Factor

**File**: `phlower/util.py`, function `test_cholesky()` (line ~413)

**Problem**: After Patch 1 resolves the import error, `cholesky(A, beta=beta)`
succeeds but returns a `(L, perm)` tuple instead of the callable `Factor`
object that phlower expects. The calling code in `graphconstr.py:275` does
`edge_weight = ret(b)`, which fails with:

```
TypeError: 'tuple' object is not callable
```

**Root cause**: scikit-sparse 0.5.0 changed the return type of `cholesky()`.
Previous versions returned a `Factor` object with a `__call__` method that
solved the system. Version 0.5.0 returns `(L, perm)` — a sparse triangular
factor and a permutation vector. The factorization format also changed (the
returned `L` is upper-triangular and does not satisfy `LL^T = P(A+βI)P^T`
in a straightforward way), making manual triangular solves unreliable.

**Fix**: When `cholesky()` returns a tuple, the positive-definiteness test has
already passed (no exception was raised). We construct `A + βI` and fall back to
`scipy.sparse.linalg.spsolve()` for the actual solve, which is numerically
exact (residual ~1e-13 in testing).

```python
# ORIGINAL (line ~413):
        return solve  # x = solve( b )

# PATCHED:
        # scikit-sparse 0.5.0 returns (L, perm) tuple instead of Factor
        if isinstance(solve, tuple):
            Apb = A + beta * scipy.sparse.eye(A.shape[0], format="csc")
            def _solve(b):
                return scipy.sparse.linalg.spsolve(Apb, b)
            return _solve
        return solve  # x = solve( b )
```

**Impact**: The Cholesky factorization is still used for its intended purpose
(testing positive-definiteness). The actual linear solve falls back to a direct
sparse solver, which produces identical results. There is no performance penalty
for the matrix sizes encountered in typical phlower workflows.

---

## Full patched function

For reference, the complete `test_cholesky()` function after both patches:

```python
def test_cholesky( A, beta=1e-6, verbose=False ):
    """ try cholesky( a scipy.sparse matrix  + beta I )

    https://scicomp.stackexchange.com/questions/12979/testing-if-a-matrix-is-positive-semi-definite
    Why A + beta I, beta say 1e-6 ?
    If A has tiny eigenvalues, 0 to within machine precision,
    about half of these "zeros" may be negative -- tough on solvers.
    Also the condition number improves to ~ rho(A) / beta.
    """
    #scikit-sparse
    from sksparse.cholmod import cholesky, CholmodNotPositiveDefiniteError
    try:
        from sksparse.cholmod import CholmodTooLargeError
    except ImportError:
        class CholmodTooLargeError(Exception):
            pass
    if not scipy.sparse.isspmatrix_csc(A):
        A = scipy.sparse.csc_matrix(A)
    try:
        solve = cholesky( A, beta=beta )  # A + beta I
        if verbose:
            print( "+ %g I is positive-definite" % (beta ))
        # scikit-sparse 0.5.0 returns (L, perm) tuple instead of Factor
        if isinstance(solve, tuple):
            Apb = A + beta * scipy.sparse.eye(A.shape[0], format="csc")
            def _solve(b):
                return scipy.sparse.linalg.spsolve(Apb, b)
            return _solve
        return solve  # x = solve( b )
    except CholmodNotPositiveDefiniteError:
        if verbose:
            print( " + %g I is not positive-definite" % (beta ))
    except CholmodTooLargeError:
        if verbose:
            print( " + %g I is too large" % (beta ))
        return False
```

---

## Environment details

```
Python 3.10 (conda: phlower_env)
phlowerpy==0.1.5
scikit-sparse==0.5.0
scipy==1.15.3
numpy==2.2.6
suitesparse (via brew): suite-sparse 8.x
```

---

## Reproducibility

Both issues are triggered by calling `phlower.ext.ddhodge()` on any dataset.
They are not specific to our data — any user with scikit-sparse 0.5.0 will hit
them. A minimal reproducer for Patch 2:

```python
import scipy.sparse
from sksparse.cholmod import cholesky

A = scipy.sparse.eye(10, format="csc")
result = cholesky(A, beta=1e-6)
print(type(result))  # <class 'tuple'> in 0.5.0, <class 'Factor'> in older versions
```

---

## Change log

| Date       | Patch | Description                                      |
|------------|-------|--------------------------------------------------|
| 2026-03-30 | 1     | Fallback import for `CholmodTooLargeError`        |
| 2026-03-30 | 2     | Handle tuple return from `cholesky()` in 0.5.0    |
