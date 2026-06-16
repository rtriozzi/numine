# showerreco/geometry.py

# Lightweight replacements for LArSoft geo::Point_t / geo::Vector_t and the
# ShowerElementHolder that ferries outputs between tools in the Pandora outerfcace.
# All coordinates are in the LArSoft convention (x = drift, y = vertical, z = beam).
# Units are centimetres throughout unless otherwise noted.

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy

Point  = numpy.ndarray  # shape (3,), dtype float64 — [x, y, z]
def point(x: float, y: float, z: float) -> Point:
    return numpy.array([x, y, z], dtype=numpy.float64)
UNSET_POINT  = point(-999., -999., -999.)

Vector = numpy.ndarray  # shape (3,), dtype float64 — [dx, dy, dz]
def vector(dx: float, dy: float, dz: float) -> Vector:
    return numpy.array([dx, dy, dz], dtype=numpy.float64)
UNSET_VECTOR = vector(-999., -999., -999.)

@dataclass
class SpacePoint:
    position: Point       # [x, y, z] in cm
    charge:   float = 0.0 # integrated ADC (arbitrary units)
    time:     float = 0.0 # drift time in µs (used for lifetime correction)

@dataclass
class PCAxis:
    ok:           bool
    n_hits:       int
    eigenvalues:  numpy.ndarray # shape (3,), descending
    eigenvectors: numpy.ndarray # shape (3, 3), rows are eigenvectors
    ave_position: Point         # charge-weighted centroid

class ShowerElementHolder:
    '''
        Mimics the Pandora outerface `ShowerElementHolder`: a typed key-value store 
        that tools write to and read from.  
        Values and their errors are stored separately.
    '''

    def __init__(self) -> None:
        self._elements:     dict[str, Any]  = {}
        self._errors:       dict[str, Any]  = {}
        self._set_flags:    dict[str, bool] = {}

    # writing
    def set_element(self, key: str, value: Any, error: Any = None) -> None:
        self._elements[key]   = value
        self._errors[key]     = error
        self._set_flags[key]  = True

    # reading 
    def check_element(self, key: str) -> bool:
        return self._set_flags.get(key, False)

    def get_element(self, key: str) -> Any:
        if not self.check_element(key):
            raise KeyError(f"ShowerElementHolder: element '{key}' has not been set.")
        return self._elements[key]

    def get_error(self, key: str) -> Any:
        if not self.check_element(key):
            raise KeyError(f"ShowerElementHolder: element '{key}' has not been set.")
        return self._errors[key]

    # convenience
    def __repr__(self) -> str:
        keys = [k for k, v in self._set_flags.items() if v]
        return f"ShowerElementHolder(set={keys})"