"""Load a COBRA Toolbox .mat model into plain numpy/scipy arrays.

This avoids depending on cobrapy, so the same model data can be fed to a
reference implementation, to cobrapy (via an adapter), or to anything else.
"""
from __future__ import annotations

import numpy as np
import scipy.io as sio
import scipy.sparse as sp


def _col(x, n, default):
    if x is None:
        return np.full(n, default, dtype=float)
    return np.asarray(x, dtype=float).ravel()


def _strs(x):
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return [str(s) for s in np.asarray(x, dtype=object).ravel()]


def load_cobra_mat(path: str, var: str | None = None) -> dict:
    """Return a dict with S (csr), b, c, lb, ub, csense, rxns, mets, osenseStr."""
    raw = sio.loadmat(path, simplify_cells=True)
    names = [k for k in raw if not k.startswith("__")]
    if var is None:
        if len(names) != 1:
            raise ValueError(f"{path} holds {names}; pass var=")
        var = names[0]
    m = raw[var]
    S = sp.csr_matrix(m["S"], dtype=float)
    nm, nr = S.shape
    csense = m.get("csense")
    if csense is None or len(np.atleast_1d(csense)) == 0:
        csense = "E" * nm
    csense = "".join(np.atleast_1d(csense).astype(str).ravel()).replace(" ", "")
    model = {
        "S": S,
        "b": _col(m.get("b"), nm, 0.0),
        "c": _col(m.get("c"), nr, 0.0),
        "lb": _col(m.get("lb"), nr, -1000.0),
        "ub": _col(m.get("ub"), nr, 1000.0),
        "csense": csense[:nm],
        "rxns": _strs(m.get("rxns")),
        "mets": _strs(m.get("mets")),
        "osenseStr": m.get("osenseStr") or "max",
        "name": var,
        "path": path,
    }
    assert len(model["rxns"]) == nr and len(model["mets"]) == nm
    return model
