"""Array API namespace entry point.

Code that should run on any Array API compatible library asks for the namespace of its
inputs through :func:`xp_of` instead of importing NumPy directly.
"""

from array_api_compat import array_namespace

__all__ = ["xp_of"]


def xp_of(*arrays):
    """Return the Array API namespace shared by ``arrays``.

    This is ``array_api_compat.array_namespace(*arrays)``: NumPy arrays map to the
    ``array_api_compat.numpy`` wrapper, and libraries that already follow the standard
    (for example ``array_api_strict``) map to themselves. Mixing arrays from different
    libraries raises ``TypeError``.
    """
    return array_namespace(*arrays)
